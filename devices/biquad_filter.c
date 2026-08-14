#include <math.h>
#include <gsl/gsl_matrix.h>

#include "anyloop.h"
#include "logging.h"
#include "xalloc.h"
#include "biquad_filter.h"


// Compute normalized Direct Form I coefficients via the bilinear transform.
// Frequency pre-warping is baked in: because we build the coefficients from
// cos(w0)/sin(w0) of the *digital* frequency w0 = 2*pi*f0/fs (rather than from
// the raw analog frequency), the realized feature lands exactly on f0 even when
// f0 is a large fraction of Nyquist. See doc/devices/biquad_filter.md.
// compute one normalized coefficient set (slot k: 0=y/tip, 1=x/tilt) from f0/q
static int biquad_filter_compute_one(
	struct aylp_biquad_filter_data *data, size_t k, double f0, double q
){
	if (data->fs <= 0.0) {
		log_error("biquad_filter: fs must be > 0 (got %G)", data->fs);
		return -1;
	}
	if (f0 <= 0.0 || f0 >= data->fs / 2.0) {
		log_error("biquad_filter: f0 must be in (0, fs/2)=(0, %G); "
			"got %G", data->fs / 2.0, f0
		);
		return -1;
	}
	if (q <= 0.0) {
		log_error("biquad_filter: q must be > 0 (got %G)", q);
		return -1;
	}

	double w0 = 2.0 * M_PI * f0 / data->fs;
	double cw = cos(w0);
	double sw = sin(w0);
	double alpha = sw / (2.0 * q);

	// shared denominator (poles set width/damping)
	double a0 = 1.0 + alpha;
	double a1 = -2.0 * cw;
	double a2 = 1.0 - alpha;

	// numerator (zeros) depends on the mode
	double b0, b1, b2;
	switch (data->mode) {
	case AYLP_BIQUAD_FILTER_NOTCH:
		b0 = 1.0;
		b1 = -2.0 * cw;
		b2 = 1.0;
		break;
	case AYLP_BIQUAD_FILTER_LOWPASS:
	default:
		b0 = (1.0 - cw) / 2.0;
		b1 = 1.0 - cw;
		b2 = (1.0 - cw) / 2.0;
		break;
	}

	// normalize so a0 == 1
	data->b0[k] = b0 / a0;
	data->b1[k] = b1 / a0;
	data->b2[k] = b2 / a0;
	data->a1[k] = a1 / a0;
	data->a2[k] = a2 / a0;
	return 0;
}


// compute both coefficient sets: [1] from the base f0/q (x/tilt and matrix),
// [0] from the per-axis f0_y/q_y overrides (y/tip), falling back to the base.
static int biquad_filter_compute_coeffs(struct aylp_biquad_filter_data *data)
{
	double f0_y = data->f0_y > 0.0 ? data->f0_y : data->f0;
	double q_y = data->q_y > 0.0 ? data->q_y : data->q;
	if (biquad_filter_compute_one(data, 1, data->f0, data->q)) return -1;
	if (biquad_filter_compute_one(data, 0, f0_y, q_y)) return -1;
	return 0;
}


// (re)allocate per-element state if the pipeline element count changed
static void biquad_filter_alloc_state(
	struct aylp_biquad_filter_data *data, size_t n
){
	if (data->n == n) return;
	xfree(data->x1);
	xfree(data->x2);
	xfree(data->y1);
	xfree(data->y2);
	data->x1 = xcalloc(n, sizeof(double));
	data->x2 = xcalloc(n, sizeof(double));
	data->y1 = xcalloc(n, sizeof(double));
	data->y2 = xcalloc(n, sizeof(double));
	data->n = n;
}


// one Direct Form I step for element i, given input x; returns filter output.
// k selects the coefficient set: 0 for element 0 (y/tip), 1 otherwise (x/tilt).
static inline double biquad_filter_step(
	struct aylp_biquad_filter_data *data, size_t i, size_t k, double x
){
	double y =
		  data->b0[k] * x
		+ data->b1[k] * data->x1[i]
		+ data->b2[k] * data->x2[i]
		- data->a1[k] * data->y1[i]
		- data->a2[k] * data->y2[i]
	;
	data->x2[i] = data->x1[i];
	data->x1[i] = x;
	data->y2[i] = data->y1[i];
	data->y1[i] = y;
	return y;
}


int biquad_filter_init(struct aylp_device *self)
{
	self->proc = &biquad_filter_proc;
	self->fini = &biquad_filter_fini;
	self->device_data = xcalloc(1, sizeof(struct aylp_biquad_filter_data));
	struct aylp_biquad_filter_data *data = self->device_data;

	// defaults
	data->mode = AYLP_BIQUAD_FILTER_LOWPASS;
	data->q = M_SQRT1_2;	// ~0.707, maximally-flat (Butterworth) lowpass

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
		} else if (!strcmp(key, "mode")) {
			const char *s = json_object_get_string(val);
			if (!strcmp(s, "lowpass"))
				data->mode = AYLP_BIQUAD_FILTER_LOWPASS;
			else if (!strcmp(s, "notch"))
				data->mode = AYLP_BIQUAD_FILTER_NOTCH;
			else log_error("Unrecognized mode: %s", s);
			log_trace("mode = %s", s);
		} else if (!strcmp(key, "f0")) {
			data->f0 = json_object_get_double(val);
			log_trace("f0 = %G", data->f0);
		} else if (!strcmp(key, "q")) {
			data->q = json_object_get_double(val);
			log_trace("q = %G", data->q);
		} else if (!strcmp(key, "f0y")) {
			data->f0_y = json_object_get_double(val);
			log_trace("f0y = %G", data->f0_y);
		} else if (!strcmp(key, "qy")) {
			data->q_y = json_object_get_double(val);
			log_trace("qy = %G", data->q_y);
		} else if (!strcmp(key, "fs")) {
			data->fs = json_object_get_double(val);
			log_trace("fs = %G", data->fs);
		} else {
			log_warn("Unknown parameter \"%s\"", key);
		}
	}

	if (!data->type) {
		log_error("You must provide valid type param.");
		return -1;
	}
	if (biquad_filter_compute_coeffs(data)) return -1;

	double f0_y = data->f0_y > 0.0 ? data->f0_y : data->f0;
	double q_y = data->q_y > 0.0 ? data->q_y : data->q;
	log_info("biquad_filter: %s fs=%G Hz; "
		"x/tilt f0=%G Hz Q=%G b=[%G %G %G] a=[1 %G %G]; "
		"y/tip f0=%G Hz Q=%G b=[%G %G %G] a=[1 %G %G]",
		data->mode == AYLP_BIQUAD_FILTER_NOTCH ? "notch" : "lowpass",
		data->fs,
		data->f0, data->q,
		data->b0[1], data->b1[1], data->b2[1], data->a1[1], data->a2[1],
		f0_y, q_y,
		data->b0[0], data->b1[0], data->b2[0], data->a1[0], data->a2[0]
	);

	// set types and units
	self->type_in = data->type;
	self->units_in = AYLP_U_ANY;
	self->type_out = AYLP_T_UNCHANGED;
	self->units_out = data->units;
	return 0;
}


int biquad_filter_proc(struct aylp_device *self, struct aylp_state *state)
{
	struct aylp_biquad_filter_data *data = self->device_data;

	switch (data->type) {
	case AYLP_T_VECTOR: {
		gsl_vector *s = state->vector;
		if (UNLIKELY(!data->res_v || data->res_v->size != s->size)) {
			if (data->res_v) xfree_type(gsl_vector, data->res_v);
			data->res_v = xmalloc_type(gsl_vector, s->size);
		}
		biquad_filter_alloc_state(data, s->size);
		gsl_vector *r = data->res_v;
		for (size_t j = 0; j < s->size; j++) {
			double x = s->data[j * s->stride];
			// element 0 is y/tip (coeff set 0); the rest are x/tilt
			size_t k = j ? 1 : 0;
			r->data[j * r->stride] = biquad_filter_step(data, j, k, x);
		}
		state->vector = r;
		break;
	}
	case AYLP_T_MATRIX: {
		gsl_matrix *s = state->matrix;
		if (UNLIKELY(!data->res_m
		|| data->res_m->size1 != s->size1
		|| data->res_m->size2 != s->size2)) {
			if (data->res_m) xfree_type(gsl_matrix, data->res_m);
			data->res_m = xmalloc_type(gsl_matrix,
				s->size1, s->size2
			);
		}
		biquad_filter_alloc_state(data, s->size1 * s->size2);
		gsl_matrix *r = data->res_m;
		size_t i = 0;
		for (size_t y = 0; y < s->size1; y++) {
			for (size_t x = 0; x < s->size2; x++) {
				double in = s->data[y * s->tda + x];
				// matrix input has no axis split; use base set
				r->data[y * r->tda + x] =
					biquad_filter_step(data, i++, 1, in);
			}
		}
		state->matrix = r;
		break;
	}
	}

	return 0;
}


int biquad_filter_fini(struct aylp_device *self)
{
	struct aylp_biquad_filter_data *data = self->device_data;
	xfree(data->x1);
	xfree(data->x2);
	xfree(data->y1);
	xfree(data->y2);
	switch (data->type) {
	case AYLP_T_VECTOR:
		if (data->res_v) xfree_type(gsl_vector, data->res_v);
		break;
	case AYLP_T_MATRIX:
		if (data->res_m) xfree_type(gsl_matrix, data->res_m);
		break;
	}
	xfree(data);
	return 0;
}

