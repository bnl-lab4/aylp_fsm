#include <time.h>
#include <math.h>
#include <stdbool.h>

#include "anyloop.h"
#include "logging.h"
#include "xalloc.h"
#include "pid.h"


/** One step of the series lead section, as a pure function of its input and the
* stored filter state. Keeping it pure lets the anti-windup logic below evaluate
* a candidate command before deciding whether to commit the filter state. */
static inline double lead_eval(bool on, double b0, double b1, double a1,
	double raw, double lead_in, double lead_out)
{
	return on ? (b0*raw + b1*lead_in - a1*lead_out) : raw;
}


/** Largest |accumulator| whose integral term -i*a still fits inside the output
* clamp. With i == 0 the accumulator does not reach the output at all, so just
* keep it bounded. */
static inline double acc_limit(double i, double clamp)
{
	return (i != 0.0) ? clamp / fabs(i) : clamp;
}


/** Parse a scalar or array param into dst[0..max); returns the count. */
static size_t parse_double_list(struct json_object *val, double *dst, size_t max)
{
	if (json_object_is_type(val, json_type_array)) {
		size_t n = json_object_array_length(val);
		if (n > max) {
			log_warn("pid: list truncated to %zu entries", max);
			n = max;
		}
		for (size_t k = 0; k < n; k++)
			dst[k] = json_object_get_double(
				json_object_array_get_idx(val, k));
		return n;
	}
	dst[0] = json_object_get_double(val);
	return 1;
}


int pid_init(struct aylp_device *self)
{
	int err;
	self->proc = &pid_proc;
	self->fini = &pid_fini;
	self->device_data = xcalloc(1, sizeof(struct aylp_pid_data));
	struct aylp_pid_data *data = self->device_data;

	// set defaults
	data->p = 1.0;
	data->i = 0.0;
	data->d = 0.0;
	data->clamp = 1.0;
	data->start_delay = 0.0;	// no startup hold by default
	data->g = 1.0;
	data->p_y = data->i_y = data->d_y = data->g_y = -1.0;
	data->dfilt = 0.0;	// derivative filter off by default
	data->dfilt_y = -1.0;	// <0 means "use base dfilt"
	data->deadband = 0.0;	// off by default
	data->deadband_y = -1.0;	// <0 means "use base deadband"
	// lead compensator off by default (0 freqs disable it), unity gain
	data->lead_fz = data->lead_fp = 0.0;
	data->lead_g = 1.0;
	data->lead_fz_y = data->lead_fp_y = 0.0;
	data->lead_g_y = 1.0;

	// line-rejection staging; assembled into data->line after the parse.
	// index [1] = x (base params), [0] = y, matching the p/p_y convention
	double lfreq[2][AYLP_PID_MAX_LINES];
	double lphase[2][AYLP_PID_MAX_LINES];
	size_t n_lphase[2] = {0, 0};
	double lgain[2][AYLP_PID_MAX_LINES];
	size_t n_lgain[2] = {0, 0};

	// parse parameters
	if (!self->params) {
		log_error("No params object found.");
		return -1;
	}
	json_object_object_foreach(self->params, key, val) {
		if (key[0] == '_') {
			// keys starting with _ are comments
		} else if (!strcmp(key, "type")) {
			const char *s = json_object_get_string(val);
			if (!strcmp(s, "vector"))
				data->type = AYLP_T_VECTOR;
			else if (!strcmp(s, "matrix"))
				data->type = AYLP_T_MATRIX;
			else log_error("Unrecognized type: %s", s);
			log_trace("type = %s (0x%hhX)", s, data->type);
		} else if (!strcmp(key, "units")) {
			const char *s = json_object_get_string(val);
			data->units = aylp_units_from_string(s);
			log_trace("units = %s (0x%hhX)", s, data->units);
		} else if (!strcmp(key, "p")) {
			data->p = json_object_get_double(val);
			log_trace("p = %G", data->p);
		} else if (!strcmp(key, "i")) {
			data->i = json_object_get_double(val);
			log_trace("i = %G", data->i);
		} else if (!strcmp(key, "d")) {
			data->d = json_object_get_double(val);
			log_trace("d = %G", data->d);
		} else if (!strcmp(key, "g")) {
			data->g = json_object_get_double(val);
			log_trace("g = %G", data->g);
		} else if (!strcmp(key, "py")) { data->p_y = json_object_get_double(val);
		} else if (!strcmp(key, "iy")) { data->i_y = json_object_get_double(val);
		} else if (!strcmp(key, "dy")) { data->d_y = json_object_get_double(val);
		} else if (!strcmp(key, "gy")) { data->g_y = json_object_get_double(val);
		} else if (!strcmp(key, "dfilt")) {
			data->dfilt = json_object_get_double(val);
			log_trace("dfilt = %G Hz", data->dfilt);
		} else if (!strcmp(key, "dfilty")) {
			data->dfilt_y = json_object_get_double(val);
			log_trace("dfilty = %G Hz", data->dfilt_y);
		} else if (!strcmp(key, "deadband")) {
			data->deadband = json_object_get_double(val);
			if (data->deadband < 0) data->deadband = -data->deadband;
			log_trace("deadband = %G", data->deadband);
		} else if (!strcmp(key, "deadbandy")) {
			data->deadband_y = json_object_get_double(val);
			log_trace("deadbandy = %G", data->deadband_y);
		} else if (!strcmp(key, "start_delay")) {
			data->start_delay = json_object_get_double(val);
			if (data->start_delay < 0) data->start_delay = 0;
			log_trace("start_delay = %G s", data->start_delay);
		} else if (!strcmp(key, "park_after")) {
			data->park_after = json_object_get_uint64(val);
			log_trace("park_after = %zu", data->park_after);
		} else if (!strcmp(key, "park_value")) {
			data->park_value = json_object_get_double(val);
			log_trace("park_value = %G", data->park_value);
		} else if (!strcmp(key, "clamp")) {
			data->clamp = json_object_get_double(val);
			if (data->clamp < 0)
				data->clamp = -data->clamp;
			log_trace("clamp = ±%G", data->clamp);
		} else if (!strcmp(key, "lead_fz")) {
			data->lead_fz = json_object_get_double(val);
			log_trace("lead_fz = %G", data->lead_fz);
		} else if (!strcmp(key, "lead_fp")) {
			data->lead_fp = json_object_get_double(val);
			log_trace("lead_fp = %G", data->lead_fp);
		} else if (!strcmp(key, "lead_g")) {
			data->lead_g = json_object_get_double(val);
			log_trace("lead_g = %G", data->lead_g);
		} else if (!strcmp(key, "lead_fzy")) {
			data->lead_fz_y = json_object_get_double(val);
			log_trace("lead_fzy = %G", data->lead_fz_y);
		} else if (!strcmp(key, "lead_fpy")) {
			data->lead_fp_y = json_object_get_double(val);
			log_trace("lead_fpy = %G", data->lead_fp_y);
		} else if (!strcmp(key, "lead_gy")) {
			data->lead_g_y = json_object_get_double(val);
			log_trace("lead_gy = %G", data->lead_g_y);
		} else if (!strcmp(key, "lines")) {
			data->n_lines[1] = parse_double_list(val, lfreq[1],
				AYLP_PID_MAX_LINES);
		} else if (!strcmp(key, "linesy")) {
			data->n_lines[0] = parse_double_list(val, lfreq[0],
				AYLP_PID_MAX_LINES);
		} else if (!strcmp(key, "line_gain")) {
			n_lgain[1] = parse_double_list(val, lgain[1],
				AYLP_PID_MAX_LINES);
		} else if (!strcmp(key, "line_gainy")) {
			n_lgain[0] = parse_double_list(val, lgain[0],
				AYLP_PID_MAX_LINES);
		} else if (!strcmp(key, "line_delay")) {
			data->line_delay = json_object_get_double(val);
			log_trace("line_delay = %G s", data->line_delay);
		} else if (!strcmp(key, "line_phase")) {
			n_lphase[1] = parse_double_list(val, lphase[1],
				AYLP_PID_MAX_LINES);
		} else if (!strcmp(key, "line_phasey")) {
			n_lphase[0] = parse_double_list(val, lphase[0],
				AYLP_PID_MAX_LINES);
		} else {
			log_warn("Unknown parameter \"%s\"", key);
		}
	}
	if (data->p_y < 0) data->p_y = data->p;
	if (data->i_y < 0) data->i_y = data->i;
	if (data->d_y < 0) data->d_y = data->d;
	if (data->g_y < 0) data->g_y = data->g;
	// assemble the line-rejection oscillators; a gain list shorter than
	// the line list extends with its last element (scalar = same for all)
	if ((data->n_lines[0] || data->n_lines[1])
			&& data->type != AYLP_T_VECTOR) {
		log_error("line rejection needs type=vector.");
		return -1;
	}
	for (int ax = 0; ax < 2; ax++) for (size_t k = 0;
			k < data->n_lines[ax]; k++) {
		struct aylp_pid_line *ln = &data->line[ax][k];
		ln->f = lfreq[ax][k];
		if (ln->f <= 0.0) {
			log_error("line frequency must be > 0; got %G", ln->f);
			return -1;
		}
		if (k < n_lgain[ax]) ln->g = lgain[ax][k];
		else if (n_lgain[ax]) ln->g = lgain[ax][n_lgain[ax]-1];
		else if (n_lgain[1]) ln->g = lgain[1][n_lgain[1]-1];
		else ln->g = 20.0;
		ln->phi = (k < n_lphase[ax])
			? lphase[ax][k] * M_PI / 180.0
			: -2.0 * M_PI * ln->f * data->line_delay;
		log_info("pid: %c line rejection at %G Hz, gain %G/s, "
			"path phase %.1f deg", ax ? 'x' : 'y', ln->f, ln->g,
			ln->phi * 180.0 / M_PI);
	}
	if (data->dfilt_y < 0) data->dfilt_y = data->dfilt;
	if (data->deadband_y < 0) data->deadband_y = data->deadband;
	// make sure we didn't miss any params
	if (!data->type) {
		log_error("You must provide valid type param.");
		return -1;
	}
	// the pole enables the lead; a zero without a pole can't be realized
	if (data->lead_fz > 0 && data->lead_fp <= 0)
		log_warn("x lead has a zero but no pole (lead_fp); lead disabled.");
	if (data->lead_fz_y > 0 && data->lead_fp_y <= 0)
		log_warn("y lead has a zero but no pole (lead_fpy); lead disabled.");
	int lead_on = data->lead_fp > 0 || data->lead_fp_y > 0;
	if (lead_on && data->type != AYLP_T_VECTOR)
		log_warn("Lead compensator only applies for type=vector; ignoring.");

	// allocate dummy vectors or matrices so we can skip checking if they
	// exist in proc()
	switch (data->type) {
	case AYLP_T_VECTOR:
		data->acc_v = xmalloc_type(gsl_vector, 0);
		data->pre_v = xmalloc_type(gsl_vector, 0);
		data->res_v = xmalloc_type(gsl_vector, 0);
		data->dfd_v = xmalloc_type(gsl_vector, 0);
		data->lead_in_v = xmalloc_type(gsl_vector, 0);
		data->lead_out_v = xmalloc_type(gsl_vector, 0);
		break;
	case AYLP_T_MATRIX:
		data->acc_m = xmalloc_type(gsl_matrix, 0, 0);
		data->pre_m = xmalloc_type(gsl_matrix, 0, 0);
		data->res_m = xmalloc_type(gsl_matrix, 0, 0);
		data->dfd_m = xmalloc_type(gsl_matrix, 0, 0);
		break;
	}

	err = clock_gettime(CLOCK_MONOTONIC, &data->tp);
	if (err) {
		log_error("Couldn't get time: %s", strerror(err));
		return -1;
	}

	// set types and units
	self->type_in = data->type;
	self->units_in = AYLP_U_ANY;
	self->type_out = AYLP_T_UNCHANGED;
	self->units_out = data->units;
	return 0;
}


/** Empty every piece of controller memory.
*
* Called when the output parks: each of these holds a correction derived from
* measurements the sensor has told us not to believe, so carrying any of it
* across the outage would put a stale correction on the first good frame. The
* line oscillators lose their weights but keep their phase, since their phase is
* a property of the disturbance and of wall-clock time, not of the loop. */
static void pid_reset_state(struct aylp_pid_data *data)
{
	if (data->type == AYLP_T_VECTOR) {
		if (data->acc_v) gsl_vector_set_zero(data->acc_v);
		if (data->pre_v) gsl_vector_set_zero(data->pre_v);
		if (data->dfd_v) gsl_vector_set_zero(data->dfd_v);
		if (data->lead_in_v) gsl_vector_set_zero(data->lead_in_v);
		if (data->lead_out_v) gsl_vector_set_zero(data->lead_out_v);
	} else {
		if (data->acc_m) gsl_matrix_set_zero(data->acc_m);
		if (data->pre_m) gsl_matrix_set_zero(data->pre_m);
		if (data->dfd_m) gsl_matrix_set_zero(data->dfd_m);
	}
	for (size_t j = 0; j < 2; j++)
		for (size_t k = 0; k < data->n_lines[j]; k++)
			data->line[j][k].a = data->line[j][k].b = 0.0;
}


/** Drive every element of the output to park_value and emit it.
*
* Writes into the same result object the normal path uses, so the pipeline sees
* an ordinary correction of the expected size, type and units -- nothing
* downstream needs to know that this one came from the park rather than from the
* controller. */
static int pid_park(struct aylp_device *self, struct aylp_state *state)
{
	struct aylp_pid_data *data = self->device_data;
	if (data->type == AYLP_T_VECTOR) {
		gsl_vector *s = state->vector;
		if (data->res_v->size != s->size) {
			xfree_type(gsl_vector, data->res_v);
			data->res_v = xcalloc_type(gsl_vector, s->size);
		}
		gsl_vector_set_all(data->res_v, data->park_value);
		state->vector = data->res_v;
	} else {
		gsl_matrix *s = state->matrix;
		if (data->res_m->size1 != s->size1
				|| data->res_m->size2 != s->size2) {
			xfree_type(gsl_matrix, data->res_m);
			data->res_m = xcalloc_type(gsl_matrix,
				s->size1, s->size2
			);
		}
		gsl_matrix_set_all(data->res_m, data->park_value);
		state->matrix = data->res_m;
	}
	// deliberately no header housekeeping here: pid's type_out is
	// AYLP_T_UNCHANGED, so writing it into the header would null the
	// pipeline type and every device after this one would fail on a park.
	// The normal path leaves the header alone for the same reason.
	return 0;
}


int pid_proc(struct aylp_device *self, struct aylp_state *state)
{
	int err;
	struct aylp_pid_data *data = self->device_data;
	struct timespec tp1;
	err = clock_gettime(CLOCK_MONOTONIC, &tp1);
	if (err) {
		log_error("Couldn't get time: %s", strerror(err));
		return -1;
	}
	double dt = tp1.tv_sec - data->tp.tv_sec;
	dt += 1E-9 * (tp1.tv_nsec - data->tp.tv_nsec);
	data->tp = tp1;
	log_trace("dt = %G s", dt);

	// --- park on a genuinely lost signal --------------------------------
	// Count consecutive iterations the sensor flagged as carrying no fresh
	// measurement. A held error is not an error the controller can act on:
	// the integrator winds at a constant rate on it, which is how a lost
	// beam becomes a large excursion rather than simply a stalled one.
	if (data->park_after) {
		if (LIKELY(!(state->header.status & AYLP_NO_SIGNAL))) {
			// log before zeroing the counter, so the message can
			// say how long the outage actually lasted
			if (UNLIKELY(data->parked)) {
				log_info("pid: measurement returned after %zu "
					"iterations parked; resuming control "
					"from rest", data->no_signal
				);
				data->parked = false;
			}
			data->no_signal = 0;
		} else if (++data->no_signal >= data->park_after) {
			if (!data->parked) {
				log_warn("pid: no fresh measurement for %zu "
					"iterations; parking output at %G and "
					"resetting controller state",
					data->no_signal, data->park_value
				);
				// empty the accumulator and the filters, so the
				// loop restarts from rest when the signal comes
				// back rather than dumping a stale integral
				// into the first good frame
				pid_reset_state(data);
				data->parked = true;
			}
			return pid_park(self, state);
		}
	}

	switch (data->type) {
	case AYLP_T_VECTOR: {
		gsl_vector *s = state->vector;
		// check if we need to (re)initialize
		if (data->acc_v->size != s->size) {
			xfree_type(gsl_vector, data->acc_v);
			data->acc_v = xcalloc_type(gsl_vector, s->size);
		}
		if (data->pre_v->size != s->size) {
			xfree_type(gsl_vector, data->pre_v);
			data->pre_v = xcalloc_type(gsl_vector, s->size);
		}
		if (data->res_v->size != s->size) {
			xfree_type(gsl_vector, data->res_v);
			data->res_v = xcalloc_type(gsl_vector, s->size);
		}
		// derivative-filter state, zeroed so it starts at rest
		if (data->dfd_v->size != s->size) {
			xfree_type(gsl_vector, data->dfd_v);
			data->dfd_v = xcalloc_type(gsl_vector, s->size);
		}
		// lead compensator state, zeroed so the filter starts at rest
		if (data->lead_in_v->size != s->size) {
			xfree_type(gsl_vector, data->lead_in_v);
			data->lead_in_v = xcalloc_type(gsl_vector, s->size);
		}
		if (data->lead_out_v->size != s->size) {
			xfree_type(gsl_vector, data->lead_out_v);
			data->lead_out_v = xcalloc_type(gsl_vector, s->size);
		}
		// Read these only AFTER the reallocs above. They free the old vectors
		// and install new ones, so any pointer captured beforehand dangles --
		// and the very first proc always reallocs, because init() sizes these
		// to zero.
		gsl_vector *a = data->acc_v;
		gsl_vector *p = data->pre_v;
		gsl_vector *r = data->res_v;

		// Startup hold: park the command at zero and keep the integrator
		// empty until start_delay has elapsed. A coarse stage takes time to
		// walk the beam near centre; without this the integrator spends that
		// time winding on a huge error and the fine actuator is already at its
		// rail the moment it is allowed to move.
		if (UNLIKELY(!data->started)) {
			double now = tp1.tv_sec + 1E-9 * tp1.tv_nsec;
			if (data->t0 == 0.0) data->t0 = now;
			if (now - data->t0 < data->start_delay) {
				for (size_t j = 0; j < s->size; j++) {
					a->data[j*a->stride] = 0.0;
					// track the error so the derivative doesn't
					// spike on the first live sample
					p->data[j*p->stride] = s->data[j*s->stride];
					data->dfd_v->data[
						j*data->dfd_v->stride] = 0.0;
					data->lead_in_v->data[
						j*data->lead_in_v->stride] = 0.0;
					data->lead_out_v->data[
						j*data->lead_out_v->stride] = 0.0;
				}
				gsl_vector_set_zero(r);
				// and the line-rejection weights empty
				for (int ax = 0; ax < 2; ax++)
				for (size_t k = 0; k < data->n_lines[ax]; k++)
					data->line[ax][k].a =
						data->line[ax][k].b = 0.0;
				state->vector = data->res_v;
				break;
			}
			data->started = true;
			if (data->start_delay > 0.0)
				log_info("pid: released after %.2f s startup hold",
					now - data->t0);
		}

		// loop over elements and apply PID control
		for (size_t j = 0; j < s->size; j++) {
			double pj = j ? data->p : data->p_y;
			double ij = j ? data->i : data->i_y;
			double dj = j ? data->d : data->d_y;
			double gj = j ? data->g : data->g_y;
			double dfj = j ? data->dfilt : data->dfilt_y;
			double dbj = j ? data->deadband : data->deadband_y;
			// deadband: ignore sub-threshold error so the loop stops
			// chattering on quantization/relay nonlinearities at the null
			double e = s->data[j*s->stride];
			if (e < dbj && e > -dbj) e = 0.0;
			// derivative term, optionally low-passed at dfj Hz to limit
			// noise gain: D = s/(1+s/wf) via backward Euler. Reduces to the
			// raw difference (e-e_prev)/dt as the cutoff -> infinity.
			// Computed before the integrator so the anti-windup check below
			// can evaluate the whole command, not just its integral part.
			double de = (e - p->data[j*p->stride]) / dt;
			if (dfj > 0.0) {
				double wf = 2.0*M_PI*dfj;
				double dprev = data->dfd_v->data[j*data->dfd_v->stride];
				de = (dprev + wf*(e - p->data[j*p->stride]))
					/ (1.0 + wf*dt);
			}
			data->dfd_v->data[j*data->dfd_v->stride] = de;
			// per-axis series lead compensator on the PID output:
			//   C(s) = g*(s + 2pi*fz) / (s + 2pi*fp)
			// bilinear-transformed with the live dt so it tracks the
			// variable loop rate. The pole (fp>0) enables the section;
			// the zero may sit at 0 Hz (fz==0) for a band-limited
			// differentiator s/(s+wp). DC gain = g*fz/fp, HF gain = g.
			double fzj = j ? data->lead_fz : data->lead_fz_y;
			double fpj = j ? data->lead_fp : data->lead_fp_y;
			double glj = j ? data->lead_g  : data->lead_g_y;
			bool lead_on = fpj > 0.0;
			double b0 = 0.0, b1 = 0.0, a1 = 0.0;
			if (lead_on) {
				double k = 2.0 / dt;            // bilinear 2/T
				double wz = 2.0*M_PI*fzj;
				double wp = 2.0*M_PI*fpj;
				double a0 = k + wp;             // normalizing term
				b0 = glj * (k + wz) / a0;
				b1 = glj * (wz - k) / a0;
				a1 = (wp - k) / a0;
			}
			double li = data->lead_in_v->data[j*data->lead_in_v->stride];
			double lo = data->lead_out_v->data[j*data->lead_out_v->stride];

			// --- integrator, with anti-windup ---------------------------
			// The accumulator reaches the command as -ij*a, but the command
			// is clamped to +-clamp. Bounding `a` at +-clamp (as this used
			// to) therefore lets the integral term range over +-ij*clamp --
			// ij times past anything the output can express. After a
			// saturating transient the accumulator has to unwind all of that
			// excess before the command even leaves the rail, so recovery
			// takes ij times longer than it should. That is windup, and at
			// ij=300 it is a factor of 300. Two guards:
			//   1. bound the accumulator so |ij*a| <= clamp;
			//   2. conditional integration -- when the command is already
			//      railed, don't integrate a sample that pushes it further
			//      into the rail (the leak still applies).
			double a_prev = a->data[j*a->stride];
			double a_leak = gj * a_prev;	// leak, before this sample's error
			double a_try = a_leak + dt * e;
			// the integral term reaches the output THROUGH the lead,
			// whose DC gain is g*fz/fp -- bound the accumulator by what
			// the output can express after that attenuation, or a lead
			// with fz<fp throttles the integrator's authority to
			// clamp*fz/fp and the loop parks off-setpoint (seen
			// 2026-07-14 with lead 8/32: x stuck ~2 px off, integral
			// path pinned at clamp/4)
			double dcj = (lead_on && fzj > 0.0) ? glj*fzj/fpj : 1.0;
			double a_max = acc_limit(ij*dcj, data->clamp);
			if (a_try > a_max) a_try = a_max;
			else if (a_try < -a_max) a_try = -a_max;

			double raw = - pj*e - ij*a_try - dj*de;
			double out = lead_eval(lead_on, b0, b1, a1, raw, li, lo);
			// how far this sample's integration moved the command
			double push = -ij * (a_try - a_leak);
			if (UNLIKELY((out > data->clamp && push > 0.0)
					|| (out < -data->clamp && push < 0.0))) {
				a_try = a_leak;
				raw = - pj*e - ij*a_try - dj*de;
				out = lead_eval(lead_on, b0, b1, a1, raw, li, lo);
			}
			a->data[j*a->stride] = a_try;
			if (lead_on) {
				data->lead_in_v->data[j*data->lead_in_v->stride] = raw;
				data->lead_out_v->data[j*data->lead_out_v->stride] = out;
			}
			// narrowband internal-model line rejection ([y, x] only):
			// each oscillator's weights integrate the error demodulated
			// at its frequency, rotated by the command->error path
			// phase (filtered-x LMS). Adding the remodulated weights to
			// the command gives infinite loop gain exactly at f, so a
			// steady line is rejected even above crossover.
			if (j < 2 && data->n_lines[j]) {
				double u = 0.0;
				// like the integrator's conditional integration:
				// don't adapt on samples the rail would distort
				bool railed = out >= data->clamp
					|| out <= -data->clamp;
				for (size_t k = 0; k < data->n_lines[j]; k++) {
					struct aylp_pid_line *ln =
						&data->line[j][k];
					ln->th += 2.0*M_PI*ln->f*dt;
					if (ln->th > 2.0*M_PI)
						ln->th -= 2.0*M_PI;
					if (LIKELY(!railed)) {
						double w = ln->g * dt * e;
						ln->a -= w*cos(ln->th + ln->phi);
						ln->b -= w*sin(ln->th + ln->phi);
						// one line can never rail the
						// command on its own
						if (ln->a > data->clamp)
							ln->a = data->clamp;
						else if (ln->a < -data->clamp)
							ln->a = -data->clamp;
						if (ln->b > data->clamp)
							ln->b = data->clamp;
						else if (ln->b < -data->clamp)
							ln->b = -data->clamp;
					}
					u += ln->a*cos(ln->th)
						+ ln->b*sin(ln->th);
				}
				out += u;
			}
			r->data[j*r->stride] = out;
			// clamp result if needed
			if (r->data[j*r->stride] > data->clamp)
				r->data[j*r->stride] = data->clamp;
			else if (r->data[j*r->stride] < -data->clamp)
				r->data[j*r->stride] = -data->clamp;
			// update previous
			p->data[j*p->stride] = e;
		}
		state->vector = data->res_v;
		break;
	}

	case AYLP_T_MATRIX: {
		gsl_matrix *s = state->matrix;
		// check if we need to (re)initialize
		if (UNLIKELY(data->acc_m->size1 != s->size1
				|| data->acc_m->size2 != s->size2)) {
			xfree_type(gsl_matrix, data->acc_m);
			data->acc_m = xcalloc_type(gsl_matrix,
				s->size1, s->size2
			);
		}
		if (UNLIKELY(data->pre_m->size1 != s->size1
				|| data->pre_m->size2 != s->size2)) {
			xfree_type(gsl_matrix, data->pre_m);
			data->pre_m = xcalloc_type(gsl_matrix,
				s->size1, s->size2
			);
		}
		if (UNLIKELY(data->res_m->size1 != s->size1
				|| data->res_m->size2 != s->size2)) {
			xfree_type(gsl_matrix, data->res_m);
			data->res_m = xmalloc_type(gsl_matrix,
				s->size1, s->size2
			);
		}
		if (UNLIKELY(data->dfd_m->size1 != s->size1
				|| data->dfd_m->size2 != s->size2)) {
			xfree_type(gsl_matrix, data->dfd_m);
			data->dfd_m = xcalloc_type(gsl_matrix,
				s->size1, s->size2
			);
		}
		// read only after the reallocs above; see the vector case
		gsl_matrix *a = data->acc_m;
		gsl_matrix *p = data->pre_m;
		gsl_matrix *r = data->res_m;

		// startup hold; see the vector case
		if (UNLIKELY(!data->started)) {
			double now = tp1.tv_sec + 1E-9 * tp1.tv_nsec;
			if (data->t0 == 0.0) data->t0 = now;
			if (now - data->t0 < data->start_delay) {
				for (size_t y = 0; y < s->size1; y++)
				for (size_t x = 0; x < s->size2; x++) {
					a->data[y*a->tda+x] = 0.0;
					p->data[y*p->tda+x] = s->data[y*s->tda+x];
					r->data[y*r->tda+x] = 0.0;
					data->dfd_m->data[y*data->dfd_m->tda+x] = 0.0;
				}
				state->matrix = data->res_m;
				break;
			}
			data->started = true;
			if (data->start_delay > 0.0)
				log_info("pid: released after %.2f s startup hold",
					now - data->t0);
		}
		// loop over elements and apply PID control
		for (size_t y = 0; y < s->size1; y++) {
			for (size_t x = 0; x < s->size2; x++) {
				// deadband: ignore sub-threshold error (see vector case)
				double e = s->data[y*s->tda+x];
				if (e < data->deadband && e > -data->deadband) e = 0.0;
				// derivative term, optionally low-passed (see vector case);
				// computed first so anti-windup sees the whole command
				double de = (e - p->data[y*p->tda+x]) / dt;
				if (data->dfilt > 0.0) {
					double wf = 2.0*M_PI*data->dfilt;
					double dprev = data->dfd_m->data[y*data->dfd_m->tda+x];
					de = (dprev + wf*(e - p->data[y*p->tda+x]))
						/ (1.0 + wf*dt);
				}
				data->dfd_m->data[y*data->dfd_m->tda+x] = de;
				// integrator with anti-windup; see the vector case for why
				// the accumulator is bounded by clamp/|i| and not by clamp
				// (note this path has never applied the `g` leak)
				double a_prev = a->data[y*a->tda+x];
				double a_try = a_prev + dt * e;
				double a_max = acc_limit(data->i, data->clamp);
				if (a_try > a_max) a_try = a_max;
				else if (a_try < -a_max) a_try = -a_max;
				double out = - data->p*e - data->i*a_try - data->d*de;
				double push = -data->i * (a_try - a_prev);
				if (UNLIKELY((out > data->clamp && push > 0.0)
						|| (out < -data->clamp && push < 0.0))) {
					a_try = a_prev;
					out = - data->p*e - data->i*a_try - data->d*de;
				}
				a->data[y*a->tda+x] = a_try;
				r->data[y*r->tda+x] = out;
				// clamp result if needed
				if (r->data[y*r->tda+x] > data->clamp)
					r->data[y*r->tda+x] = data->clamp;
				if (r->data[y*r->tda+x] < -data->clamp)
					r->data[y*r->tda+x] = -data->clamp;
				// update previous
				p->data[y*p->tda+x] = e;
			}
		}
		state->matrix = data->res_m;
		break;
	}
	}

	return 0;
}


int pid_fini(struct aylp_device *self)
{
	struct aylp_pid_data *data = self->device_data;
	switch (data->type) {
	case AYLP_T_VECTOR:
		xfree_type(gsl_vector, data->acc_v);
		xfree_type(gsl_vector, data->pre_v);
		xfree_type(gsl_vector, data->res_v);
		xfree_type(gsl_vector, data->dfd_v);
		xfree_type(gsl_vector, data->lead_in_v);
		xfree_type(gsl_vector, data->lead_out_v);
		break;
	case AYLP_T_MATRIX:
		xfree_type(gsl_matrix, data->acc_m);
		xfree_type(gsl_matrix, data->pre_m);
		xfree_type(gsl_matrix, data->res_m);
		xfree_type(gsl_matrix, data->dfd_m);
		break;
	}
	xfree(data);
	return 0;
}


