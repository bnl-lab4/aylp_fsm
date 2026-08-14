#include <time.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_vector.h>

#include "anyloop.h"
#include "logging.h"
#include "xalloc.h"
#include "kalman_filter.h"


// (re)allocate per-element state if the pipeline element count changed
static void kalman_filter_alloc_state(
	struct aylp_kalman_filter_data *data, size_t n
){
	if (data->n == n) return;
	xfree(data->hist);
	xfree(data->w);
	xfree(data->bypass_cnt);
	data->hist = xcalloc(n * data->hist_len, sizeof(double));
	data->w = xcalloc(n * data->order, sizeof(double));
	data->bypass_cnt = xcalloc(n, sizeof(size_t));
	data->n = n;
	data->head = 0;
	data->n_seen = 0;
	data->t_train = 0;
	data->engaged = false;
}

// copy the window of `order` samples ending `back` samples before the most
// recent write into xbuf, most recent first: xbuf[i] = s(t - back - i)
static inline void kalman_filter_window(
	struct aylp_kalman_filter_data *data, size_t j, size_t back,
	double *xbuf
){
	size_t L = data->hist_len;
	const double *h = data->hist + j * L;
	// data->head is the index of the most recent sample
	size_t idx = (data->head + L - back) % L;
	for (size_t i = 0; i < data->order; i++) {
		xbuf[i] = h[idx];
		idx = idx ? idx - 1 : L - 1;
	}
}


int kalman_filter_init(struct aylp_device *self)
{
	self->proc = &kalman_filter_proc;
	self->fini = &kalman_filter_fini;
	self->device_data = xcalloc(1, sizeof(struct aylp_kalman_filter_data));
	struct aylp_kalman_filter_data *data = self->device_data;

	// defaults; horizon 5 samples ~= the measured 2.1 ms loop delay at
	// the 2404 Hz 64x64 frame rate. mu 0.03 sits between the offline
	// optimum (0.02, lowest misadjustment) and faster tracking (0.05).
	data->order = 60;
	data->horizon = 5;
	data->mu = 0.03;
	data->gain = 1.0;
	data->gain_y = -1.0;	// <0 means "inherit base gain"
	data->ramp = 10.0;
	data->clamp = 1.0;
	data->start_delay = 0.0;

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
			else log_error("Unrecognized type: %s", s);
		} else if (!strcmp(key, "units")) {
			const char *s = json_object_get_string(val);
			data->units = aylp_units_from_string(s);
		} else if (!strcmp(key, "order")) {
			data->order = json_object_get_uint64(val);
			log_trace("order = %zu", data->order);
		} else if (!strcmp(key, "horizon")) {
			data->horizon = json_object_get_uint64(val);
			log_trace("horizon = %zu", data->horizon);
		} else if (!strcmp(key, "mu")) {
			data->mu = json_object_get_double(val);
			log_trace("mu = %G", data->mu);
		} else if (!strcmp(key, "gain")) {
			data->gain = json_object_get_double(val);
			log_trace("gain = %G", data->gain);
		} else if (!strcmp(key, "gainy")) {
			data->gain_y = json_object_get_double(val);
			log_trace("gainy = %G", data->gain_y);
		} else if (!strcmp(key, "ramp")) {
			data->ramp = json_object_get_double(val);
			log_trace("ramp = %G s", data->ramp);
		} else if (!strcmp(key, "clamp")) {
			data->clamp = fabs(json_object_get_double(val));
			log_trace("clamp = ±%G", data->clamp);
		} else if (!strcmp(key, "start_delay")) {
			data->start_delay = json_object_get_double(val);
			if (data->start_delay < 0) data->start_delay = 0;
			log_trace("start_delay = %G s", data->start_delay);
		} else if (!strcmp(key, "bypass")) {
			data->bypass = json_object_get_double(val);
			log_trace("bypass = %G", data->bypass);
		} else {
			log_warn("Unknown parameter \"%s\"", key);
		}
	}

	if (data->type != AYLP_T_VECTOR) {
		log_error("kalman_filter: type must be \"vector\".");
		return -1;
	}
	if (!data->order || !data->horizon) {
		log_error("kalman_filter: order and horizon must be >= 1.");
		return -1;
	}
	if (data->mu <= 0.0 || data->mu >= 2.0) {
		log_error("kalman_filter: need 0 < mu < 2 (got %G).",
			data->mu);
		return -1;
	}
	if (data->gain < 0.0 || data->gain > 1.0) {
		log_error("kalman_filter: need 0 <= gain <= 1 (got %G).",
			data->gain);
		return -1;
	}
	if (data->gain_y < 0.0) data->gain_y = data->gain;
	if (data->gain_y > 1.0) {
		log_error("kalman_filter: need 0 <= gainy <= 1 (got %G).",
			data->gain_y);
		return -1;
	}
	if (data->ramp < 0.0) data->ramp = 0.0;
	data->hist_len = data->order + data->horizon;
	data->xbuf = xcalloc(data->order, sizeof(double));

	log_info("kalman_filter: predicting %zu samples ahead, order %zu, "
		"mu %G, gain %G over a %G s ramp", data->horizon, data->order,
		data->mu, data->gain, data->ramp);

	self->type_in = AYLP_T_VECTOR;
	self->units_in = AYLP_U_ANY;
	self->type_out = AYLP_T_UNCHANGED;
	self->units_out = data->units;
	return 0;
}


int kalman_filter_proc(struct aylp_device *self, struct aylp_state *state)
{
	struct aylp_kalman_filter_data *data = self->device_data;
	gsl_vector *s = state->vector;
	size_t p = data->order;

	struct timespec tp;
	int err = clock_gettime(CLOCK_MONOTONIC, &tp);
	if (err) {
		log_error("Couldn't get time: %s", strerror(err));
		return -1;
	}
	double now = tp.tv_sec + 1E-9 * tp.tv_nsec;
	if (UNLIKELY(!data->t0)) data->t0 = now;
	// hold: pass the input through untouched and don't train, so the
	// open-loop startup transient doesn't poison the weights
	if (UNLIKELY(now - data->t0 < data->start_delay)) return 0;

	if (UNLIKELY(!data->res_v || data->res_v->size != s->size)) {
		if (data->res_v) xfree_type(gsl_vector, data->res_v);
		data->res_v = xmalloc_type(gsl_vector, s->size);
	}
	kalman_filter_alloc_state(data, s->size);
	gsl_vector *r = data->res_v;

	// advance the shared ring head once, then handle each element
	data->head = (data->head + 1) % data->hist_len;
	data->n_seen += 1;
	bool can_train = data->n_seen >= data->hist_len;
	// blend factor: ramp from 0 to gain over `ramp` seconds of training,
	// so the loop stays closed on the raw error while the zero-initialized
	// weights converge (an instant handover would feed the controller a
	// near-zero signal and leave the loop effectively open)
	double frac = 0.0;
	if (can_train) {
		if (UNLIKELY(!data->t_train)) data->t_train = now;
		double up = now - data->t_train;
		frac = data->ramp > 0.0 && up < data->ramp
			? up / data->ramp : 1.0;
		if (UNLIKELY(frac > 0.0 && !data->engaged)) {
			data->engaged = true;
			log_info("kalman_filter: training started; blending "
				"in over %G s", data->ramp);
		}
	}

	for (size_t j = 0; j < s->size; j++) {
		double z = s->data[j * s->stride];
		double out = z;
		double a = frac * (j == 0 ? data->gain_y : data->gain);
		data->hist[j * data->hist_len + data->head] = z;
		double *w = data->w + j*p;
		// transient bypass: pass the raw input through and freeze
		// training while a large disturbance is in flight, and for
		// hist_len more samples so no training pair or prediction
		// window straddles the transient (see kalman_filter.h)
		if (data->bypass > 0.0 && (z > data->bypass
				|| z < -data->bypass)) {
			if (!data->bypass_cnt[j])
				log_debug("kalman_filter: bypass ON "
					"(element %zu, z = %.4f)", j, z);
			data->bypass_cnt[j] = data->hist_len + 1;
		}
		if (data->bypass_cnt[j]) {
			data->bypass_cnt[j] -= 1;
			if (!data->bypass_cnt[j])
				log_debug("kalman_filter: bypass off "
					"(element %zu)", j);
			r->data[j * r->stride] = z;
			continue;
		}
		if (can_train) {
			// train: predict the sample just written from the
			// window ending `horizon` samples before it
			kalman_filter_window(data, j, data->horizon,
				data->xbuf);
			double e = z, xx = 1e-12;
			for (size_t i = 0; i < p; i++) {
				e -= w[i] * data->xbuf[i];
				xx += data->xbuf[i] * data->xbuf[i];
			}
			if (LIKELY(isfinite(e))) {
				double step = data->mu * e / xx;
				for (size_t i = 0; i < p; i++)
					w[i] += step * data->xbuf[i];
			} else {
				log_warn("kalman_filter: non-finite update "
					"on element %zu; resetting weights",
					j);
				memset(w, 0, p * sizeof(double));
			}
		}
		if (a > 0.0) {
			// predict `horizon` ahead from the current window
			kalman_filter_window(data, j, 0, data->xbuf);
			double yhat = 0.0;
			for (size_t i = 0; i < p; i++)
				yhat += w[i] * data->xbuf[i];
			if (LIKELY(isfinite(yhat))) {
				out = a * yhat + (1.0 - a) * z;
				if (out > data->clamp) out = data->clamp;
				if (out < -data->clamp) out = -data->clamp;
			} else {
				log_warn("kalman_filter: non-finite "
					"prediction on element %zu; "
					"resetting weights", j);
				memset(w, 0, p * sizeof(double));
			}
		}
		r->data[j * r->stride] = out;
	}

	state->vector = r;
	return 0;
}


int kalman_filter_fini(struct aylp_device *self)
{
	struct aylp_kalman_filter_data *data = self->device_data;
	xfree(data->hist);
	xfree(data->w);
	xfree(data->xbuf);
	xfree(data->bypass_cnt);
	if (data->res_v) xfree_type(gsl_vector, data->res_v);
	xfree(data);
	return 0;
}
