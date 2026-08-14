// anyloop:fsm -- the whole fast steering mirror loop as one pipeline stage.
//
// This device owns a short pipeline of its own: sensor -> controller -> clamp
// -> actuator, each an ordinary anyloop device initialized through the same
// init_device() the top level uses. A config that used to spell out
// center_of_mass, fsp, clamp and parport_dac as four pipeline entries becomes a
// camera and this.
//
// What it adds over writing those four stages out flat:
//
//  1. The substages are typechecked against each other at init, the same way
//     anyloop.c typechecks the top-level pipeline, and the composite publishes
//     the resulting type_in/type_out outward -- so it still slots into a
//     pipeline that is itself typechecked.
//  2. The mirror is PARKED when the loop exits: one zero command written
//     through the actuator before anything is closed. parport_dac_fini()
//     deliberately leaves its outputs where they are (the mirror is presumably
//     still holding a beam), and parking has so far been a separate anyloop run
//     of conf_parport_dac_park.json that run_steering.py fires before the loop.
//     Here it also happens after a crash-out or a Ctrl-C, which is when it
//     matters most.
//  3. An error out of a substage that was supposed to change the pipeline type
//     stops the loop instead of being swallowed. At the top level anyloop.c
//     exits outright in that case; inside a composite the mirror is holding a
//     beam, so we ask for a clean shutdown (AYLP_DONE) and let the park write
//     happen on the way out.
//
// It deliberately does NOT implement any beam-loss policy of its own: the
// sensors raise AYLP_BEAM_LOST/AYLP_NO_SIGNAL, and fsp's gap handling and pid's
// park_after already act on them. A second, outer policy here would fight them.

#include <gsl/gsl_vector.h>
#include <json-c/json.h>

#include "anyloop.h"
#include "logging.h"
#include "fsm.h"
#include "xalloc.h"
#include "device.h"


// Role names as they appear in the params object, and the device each role
// falls back to when the config gives params but no uri. The order of this
// array is the order the stages run in.
static const struct {
	const char *key;
	const char *def_uri;
} fsm_roles[AYLP_FSM_N_ROLES] = {
	[AYLP_FSM_SENSOR]	= { "sensor", "anyloop:center_of_mass" },
	[AYLP_FSM_CONTROLLER]	= { "controller", "anyloop:fsp" },
	[AYLP_FSM_CLAMP]	= { "clamp", "anyloop:clamp" },
	[AYLP_FSM_ACTUATOR]	= { "actuator", "anyloop:parport_dac" },
};


// A bare name like "fit_com" means the built-in device of that name; anything
// carrying a scheme ("anyloop:fit_com", "file:/opt/foo/bar.so") is passed
// through untouched. Returns an owned string.
static char *normalize_uri(const char *uri)
{
	if (strchr(uri, ':')) return xstrdup(uri);
	size_t n = strlen(uri) + sizeof("anyloop:");
	char *ret = xmalloc(n);
	snprintf(ret, n, "anyloop:%s", uri);
	return ret;
}


// Give `params` a default value for `key`, unless the key is already there (in
// which case `val` is dropped). Consumes the caller's reference to both `params`
// and `val`, and returns an owned reference to the result. Never mutates the
// caller's object: if the key has to be added we shallow-copy first, because
// `params` belongs to the config file, not to us.
static struct json_object *params_with_default(struct json_object *params,
	const char *key, struct json_object *val
){
	if (params && json_object_object_get_ex(params, key, 0)) {
		json_object_put(val);
		return params;
	}
	struct json_object *ret = json_object_new_object();
	if (params) {
		json_object_object_foreach(params, k, v) {
			json_object_object_add(ret, k, json_object_get(v));
		}
		json_object_put(params);
	}
	json_object_object_add(ret, key, val);
	return ret;
}


static struct aylp_fsm_stage *push_stage(struct aylp_fsm_data *data)
{
	data->stages = xrealloc(data->stages,
		(data->n_stages + 1) * sizeof(struct aylp_fsm_stage)
	);
	struct aylp_fsm_stage *stage = &data->stages[data->n_stages++];
	*stage = (struct aylp_fsm_stage){0};
	return stage;
}


// Add one substage from one json value. Accepts:
//   "fit_com"                        -- device name, no params
//   {"uri": ..., "params": {...}}    -- device name and its params
//   {...}                            -- params for this role's default device
static int add_stage(struct aylp_fsm_data *data, enum aylp_fsm_role role,
	struct json_object *val
){
	const char *uri = 0;
	struct json_object *params = 0;
	if (json_object_is_type(val, json_type_string)) {
		uri = json_object_get_string(val);
	} else if (json_object_is_type(val, json_type_object)) {
		struct json_object *sub;
		if (json_object_object_get_ex(val, "uri", &sub)) {
			uri = json_object_get_string(sub);
			if (json_object_object_get_ex(val, "params", &sub))
				params = json_object_get(sub);
		} else {
			// no uri key, so the whole object is this role's params
			params = json_object_get(val);
		}
	} else {
		log_error("fsm: \"%s\" must be a device name, an object, or an "
			"array of either", fsm_roles[role].key
		);
		return -1;
	}
	if (!uri || !*uri) uri = fsm_roles[role].def_uri;

	// anyloop:clamp needs to be told which type it is clamping, and inside
	// here that is always the command vector -- don't make every config
	// repeat it
	if (!strcmp(uri, "anyloop:clamp") || !strcmp(uri, "clamp")) {
		params = params_with_default(params, "type",
			json_object_new_string("vector")
		);
	}

	struct aylp_fsm_stage *stage = push_stage(data);
	stage->role = role;
	stage->dev.uri = normalize_uri(uri);
	stage->params = params;
	stage->dev.params = params;
	return 0;
}


// Add every substage configured under one role, in order.
static int add_role(struct aylp_fsm_data *data, enum aylp_fsm_role role,
	struct json_object *val
){
	if (!val || json_object_is_type(val, json_type_null)) return 0;
	if (json_object_is_type(val, json_type_array)) {
		size_t n = json_object_array_length(val);
		for (size_t idx = 0; idx < n; idx++) {
			struct json_object *el =
				json_object_array_get_idx(val, idx);
			if (json_object_is_type(el, json_type_array)) {
				log_error("fsm: \"%s\" may not nest arrays",
					fsm_roles[role].key
				);
				return -1;
			}
			if (add_stage(data, role, el)) return -1;
		}
		return 0;
	}
	return add_stage(data, role, val);
}


// Chain the substages' declared types together, exactly as anyloop.c does for
// the top-level pipeline, and work out what the composite itself accepts and
// emits. While the running type is still AYLP_T_UNCHANGED, the data in flight is
// whatever WE were handed, so those stages' requirements become our own; once a
// stage sets a concrete output type we can check the rest here and now.
static int typecheck(struct aylp_device *self, struct aylp_fsm_data *data)
{
	aylp_type type_in = AYLP_T_ANY;
	aylp_units units_in = AYLP_U_ANY;
	aylp_type type_cur = AYLP_T_UNCHANGED;
	aylp_units units_cur = AYLP_U_UNCHANGED;
	for (size_t idx = 0; idx < data->n_stages; idx++) {
		struct aylp_device *dev = &data->stages[idx].dev;
		if (!type_cur) {
			type_in &= dev->type_in;
		} else if (type_cur & ~dev->type_in) {
			log_error("fsm: %s with input type 0x%hhX is "
				"incompatible with previous type 0x%hhX",
				dev->uri, dev->type_in, type_cur
			);
			return -1;
		}
		if (!units_cur) {
			units_in &= dev->units_in;
		} else if (units_cur & ~dev->units_in) {
			log_error("fsm: %s with input units 0x%hhX is "
				"incompatible with previous units 0x%hhX",
				dev->uri, dev->units_in, units_cur
			);
			return -1;
		}
		if (dev->type_out) type_cur = dev->type_out;
		if (dev->units_out) units_cur = dev->units_out;
	}
	if (!type_in || !units_in) {
		log_error("fsm: no input type/units (0x%hhX/0x%hhX) can satisfy "
			"every stage before the first one that changes it",
			type_in, units_in
		);
		return -1;
	}
	self->type_in = type_in;
	self->units_in = units_in;
	self->type_out = type_cur;
	self->units_out = units_cur;
	return 0;
}


int fsm_init(struct aylp_device *self)
{
	self->device_data = xcalloc(1, sizeof(struct aylp_fsm_data));
	struct aylp_fsm_data *data = self->device_data;

	// attach these now so a failure partway through still gets cleaned up
	self->proc = &fsm_proc;
	self->fini = &fsm_fini;

	data->park_on_exit = true;

	if (!self->params) {
		log_error("No params object found.");
		return -1;
	}
	json_object_object_foreach(self->params, key, val) {
		// parse parameters
		if (key[0] == '_') {
			// keys starting with _ are comments
		} else if (!strcmp(key, "park_on_exit")) {
			data->park_on_exit = json_object_get_boolean(val);
			log_trace("park_on_exit = %d", data->park_on_exit);
		} else if (!strcmp(key, "park_value")) {
			data->park_value = json_object_get_double(val);
			log_trace("park_value = %g", data->park_value);
		} else {
			bool is_role = false;
			for (size_t r = 0; r < AYLP_FSM_N_ROLES; r++)
				is_role |= !strcmp(key, fsm_roles[r].key);
			// the roles themselves are handled below, in role
			// order, so that the order of the keys in the json
			// object can't silently reorder the loop
			if (!is_role)
				log_warn("Unknown parameter \"%s\"", key);
		}
	}
	for (size_t r = 0; r < AYLP_FSM_N_ROLES; r++) {
		struct json_object *val;
		if (!json_object_object_get_ex(self->params,
				fsm_roles[r].key, &val))
			continue;
		if (add_role(data, r, val)) return -1;
	}
	if (!data->n_stages) {
		log_error("fsm: no substages; give me at least one of "
			"\"sensor\", \"controller\", \"clamp\", \"actuator\""
		);
		return -1;
	}

	// the last actuator stage is the one that gets parked
	data->actuator = data->n_stages;
	for (size_t idx = 0; idx < data->n_stages; idx++)
		if (data->stages[idx].role == AYLP_FSM_ACTUATOR)
			data->actuator = idx;

	for (size_t idx = 0; idx < data->n_stages; idx++) {
		struct aylp_fsm_stage *stage = &data->stages[idx];
		log_info("fsm: initializing %s stage",
			fsm_roles[stage->role].key
		);
		if (init_device(&stage->dev)) {
			log_error("fsm: could not initialize %s",
				stage->dev.uri
			);
			return -1;
		}
		stage->inited = true;
	}

	if (typecheck(self, data)) return -1;

	if (data->park_on_exit && data->actuator == data->n_stages) {
		log_warn("fsm: park_on_exit is set but there is no actuator "
			"stage to park"
		);
		data->park_on_exit = false;
	}
	log_info("fsm: %zu stage(s), types 0x%hhX/0x%hhX -> 0x%hhX/0x%hhX, "
		"park on exit %s", data->n_stages, self->type_in,
		self->units_in, self->type_out, self->units_out,
		data->park_on_exit ? "on" : "off"
	);
	return 0;
}


int fsm_proc(struct aylp_device *self, struct aylp_state *state)
{
	struct aylp_fsm_data *data = self->device_data;
	for (size_t idx = 0; idx < data->n_stages; idx++) {
		struct aylp_device *dev = &data->stages[idx].dev;
		if (UNLIKELY(idx == data->actuator)) {
			// remember the command the actuator was handed, so we
			// can reproduce its shape when parking at exit
			data->park_header = state->header;
			data->park_size =
				state->header.type == AYLP_T_VECTOR
				? state->vector->size : 0;
		}
		if (UNLIKELY(!dev->proc)) continue;
		int err = dev->proc(dev, state);
		if (UNLIKELY(err)) {
			// same rule anyloop.c applies at the top level: an
			// error is recoverable unless the stage was supposed to
			// change the pipeline type, in which case everything
			// downstream of it -- the command we are about to put
			// on the mirror -- is garbage
			if (dev->type_out && dev->type_out != dev->type_in) {
				log_fatal("fsm: %s returned error %d but was "
					"expected to change the pipeline type; "
					"stopping the loop", dev->uri, err
				);
				aylp_set_done(&state->header);
				return err;
			}
			log_error("fsm: %s returned error %d", dev->uri, err);
		}
	}
	// If something has already asked the loop to stop, make sure the status
	// word says only that. This matters when the device that set AYLP_DONE
	// belongs to a host built before the measurement flags existed -- a
	// stock anyloop's stop_after_count, say -- because that host's exit test
	// is `status ^ AYLP_DONE`, which a sensor's AYLP_NO_SIGNAL turns into a
	// loop that will not end. Our own sensors re-raise their flag every
	// iteration, so the clearing has to happen here, after they have run:
	// the loop then exits one iteration later instead of never.
	if (UNLIKELY(state->header.status & AYLP_DONE))
		aylp_set_done(&state->header);
	return 0;
}


// Write one park command through the actuator. Called before anything is
// closed, so the actuator is still live.
static void fsm_park(struct aylp_fsm_data *data)
{
	struct aylp_device *dev = &data->stages[data->actuator].dev;
	if (!dev->proc) return;
	if (data->park_header.type != AYLP_T_VECTOR || !data->park_size) {
		log_warn("fsm: not parking %s: no command vector has reached "
			"it (last type 0x%hhX)", dev->uri,
			data->park_header.type
		);
		return;
	}
	gsl_vector *park = xcalloc_type(gsl_vector, data->park_size);
	if (data->park_value != 0.0)
		gsl_vector_set_all(park, data->park_value);
	// reuse the last header so log_dim/units/pitch match what the actuator
	// has been seeing, but clear the status flags: this is a command we are
	// issuing deliberately, not a held or lost measurement
	struct aylp_state state = {.header = data->park_header};
	state.header.status = 0;
	state.vector = park;
	log_info("fsm: parking %s at %g", dev->uri, data->park_value);
	if (dev->proc(dev, &state))
		log_error("fsm: park write to %s failed", dev->uri);
	// free our own vector, not state.vector: a substage is free to have
	// swapped in a buffer of its own
	xfree_type(gsl_vector, park);
}


int fsm_fini(struct aylp_device *self)
{
	struct aylp_fsm_data *data = self->device_data;
	if (!data) return 0;
	// the inited check matters when init failed part way through: the
	// actuator index is only meaningful once every stage has been parsed
	if (data->park_on_exit && data->actuator < data->n_stages
	&& data->stages[data->actuator].inited)
		fsm_park(data);
	for (size_t idx = 0; idx < data->n_stages; idx++) {
		struct aylp_fsm_stage *stage = &data->stages[idx];
		if (stage->inited && stage->dev.fini)
			stage->dev.fini(&stage->dev);
		if (stage->params) json_object_put(stage->params);
		if (stage->dev.uri) xfree(stage->dev.uri);
	}
	if (data->stages) xfree(data->stages);
	xfree(data);
	self->device_data = 0;
	return 0;
}
