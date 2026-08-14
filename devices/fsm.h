#ifndef AYLP_DEVICES_FSM_H_
#define AYLP_DEVICES_FSM_H_

#include <stdbool.h>
#include <gsl/gsl_vector.h>

#include "anyloop.h"

/** The roles a substage can play, in the order they are run.
* This order is the pipeline order inside the device; the order the keys happen
* to appear in the json object is irrelevant.
*/
enum aylp_fsm_role {
	AYLP_FSM_SENSOR = 0,
	AYLP_FSM_CONTROLLER,
	AYLP_FSM_CLAMP,
	AYLP_FSM_ACTUATOR,
	AYLP_FSM_N_ROLES,
};

/** One substage: an ordinary aylp_device that we own and drive ourselves. */
struct aylp_fsm_stage {
	/** The substage itself, initialized through init_device() exactly as a
	* top-level pipeline device would be. */
	struct aylp_device dev;
	/** Which role this stage was configured under. */
	enum aylp_fsm_role role;
	/** Our reference to the substage's params (put back in fini). */
	struct json_object *params;
	/** Whether dev.init() succeeded, so dev.fini() is safe to call. */
	bool inited;
};

struct aylp_fsm_data {
	/** Substages, already in the order they run. */
	struct aylp_fsm_stage *stages;
	size_t n_stages;
	/** Index of the actuator stage, or n_stages if there is none. */
	size_t actuator;
	/** Whether to drive the actuator to park_value when the loop exits. */
	bool park_on_exit;
	/** The command every element is parked at (default 0, which puts each
	* actuator channel at its own configured offset/center). */
	double park_value;
	/** The header the actuator saw on the last iteration, so the park write
	* can reproduce it. Its type stays zero until the first proc(). */
	struct aylp_header park_header;
	/** Size of the last command vector, or 0 if we never saw one. */
	size_t park_size;
};

// initialize fsm device
int fsm_init(struct aylp_device *self);

// run every substage once, in role order
int fsm_proc(struct aylp_device *self, struct aylp_state *state);

// park the mirror and close every substage when the loop exits
int fsm_fini(struct aylp_device *self);

#endif
