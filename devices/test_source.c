#include <math.h>
#include <limits.h>

#include "anyloop.h"
#include "logging.h"
#include "test_source.h"
#include "xalloc.h"


enum {KIND_CONSTANT=1, KIND_SINE=2, KIND_NOISE=3};

int test_source_init(struct aylp_device *self)
{
	self->proc = &test_source_proc;
	self->fini = &test_source_fini;
	self->device_data = xcalloc(1, sizeof(struct aylp_test_source_data));
	struct aylp_test_source_data *data = self->device_data;

	// set defaults
	data->frequency = 0.1;
	data->amplitude = 1.0;
	data->offset = 0.0;
	data->seed = 1;
	// multitone (kind sine): "frequencies" [Hz] + "fs" [Hz] -> freqs[] rad/call
	struct json_object *freqs_arr = NULL;
	double fs_hz = 0.0;

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
			else if (!strcmp(s, "matrix_uchar"))
				data->type = AYLP_T_MATRIX_UCHAR;
			else log_error("Unrecognized type: %s", s);
			log_trace("type = %s (0x%hhX)", s, data->type);
		} else if (!strcmp(key, "kind")) {
			const char *s = json_object_get_string(val);
			if (!strcmp(s, "constant")) data->kind = KIND_CONSTANT;
			else if (!strcmp(s, "sine")) data->kind = KIND_SINE;
			else if (!strcmp(s, "noise")) data->kind = KIND_NOISE;
			else log_error("Unrecognized kind: %s", s);
			log_trace("kind = %s (0x%hhX)", s, data->kind);
		} else if (!strcmp(key, "size1")) {
			data->size1 = json_object_get_uint64(val);
			log_trace("size1 = %zu", data->size1);
		} else if (!strcmp(key, "size2")) {
			data->size2 = json_object_get_uint64(val);
			log_trace("size2 = %zu", data->size2);
		} else if (!strcmp(key, "frequency")) {
			data->frequency = json_object_get_double(val);
			log_trace("frequency = %G", data->frequency);
		} else if (!strcmp(key, "amplitude")) {
			data->amplitude = json_object_get_double(val);
			log_trace("amplitude = %G", data->amplitude);
		} else if (!strcmp(key, "offset")) {
			data->offset = json_object_get_double(val);
			log_trace("offset = %G", data->offset);
		} else if (!strcmp(key, "seed")) {
			data->seed = (unsigned)json_object_get_uint64(val);
			log_trace("seed = %u", data->seed);
		} else if (!strcmp(key, "frequencies")) {
			freqs_arr = val;
		} else if (!strcmp(key, "fs")) {
			fs_hz = json_object_get_double(val);
			log_trace("fs = %G", fs_hz);
		} else {
			log_warn("Unknown parameter \"%s\"", key);
		}
	}
	// make sure we didn't miss any params
	if (!data->type) {
		log_error("You must provide a valid type param.");
		return -1;
	}
	if (!data->kind) {
		log_error("You must provide a valid kind param.");
		return -1;
	}
	if (!data->size1) {
		log_error("You must provide a valid size1 param.");
		return -1;
	}
	if (data->type & (AYLP_T_MATRIX|AYLP_T_MATRIX_UCHAR) && !data->size2) {
		log_error("You must provide a valid size2 param.");
		return -1;
	}
	// seed the white-noise PRNG (reproducible stimulus for cross-correlation)
	if (data->kind == KIND_NOISE) srand(data->seed);
	// build the multitone frequency list (Hz -> radians per proc() call)
	if (freqs_arr) {
		if (fs_hz <= 0.0) {
			log_error("\"frequencies\" needs \"fs\" (loop rate in Hz).");
			return -1;
		}
		data->nfreqs = json_object_array_length(freqs_arr);
		data->freqs = xcalloc(data->nfreqs, sizeof(double));
		for (size_t i = 0; i < data->nfreqs; i++) {
			double f_hz = json_object_get_double(
				json_object_array_get_idx(freqs_arr, i));
			data->freqs[i] = 2.0 * M_PI * f_hz / fs_hz;
			log_trace("multitone freq %zu = %G Hz", i, f_hz);
		}
	}
	// warn about clipping
	if (fabs(data->offset) + fabs(data->amplitude) > 1) {
		log_warn("Note that output will be clipped to ±1.0, "
			"or 0:255 for uchar types"
		);
	}

	switch (data->type) {
	case AYLP_T_VECTOR:
		data->vector = gsl_vector_alloc(data->size1);
		break;
	case AYLP_T_MATRIX:
		data->matrix = gsl_matrix_alloc(data->size1, data->size2);
		break;
	case AYLP_T_MATRIX_UCHAR:
		data->matrix_uchar = gsl_matrix_uchar_alloc(
			data->size1, data->size2
		);
		break;
	}

	// set types and units
	self->type_in = AYLP_T_ANY;
	self->units_in = AYLP_U_ANY;
	self->type_out = data->type;
	if (data->type == AYLP_T_MATRIX_UCHAR) self->units_out = AYLP_U_COUNTS;
	else self->units_out = AYLP_U_MINMAX;

	return 0;
}


int test_source_proc(struct aylp_device *self, struct aylp_state *state)
{
	struct aylp_test_source_data *data = self->device_data;

	double val = 0;
	switch (data->kind) {
	case KIND_CONSTANT:
		val = data->offset;
		break;
	case KIND_SINE:
		if (data->nfreqs) {
			// multitone: sum of equal-amplitude sines, so one in-loop run
			// integrates every tone over the whole record (max SNR/freq)
			double s = 0.0;
			for (size_t i = 0; i < data->nfreqs; i++)
				s += sin(data->freqs[i] * data->acc);
			val = data->offset
				+ data->amplitude / data->nfreqs * s;
		} else {
			val = data->offset + data->amplitude * sin(
				data->frequency * data->acc
			);
		}
		break;
	case KIND_NOISE:
		// white uniform noise in [offset-amp, offset+amp]; a white command
		// makes the command<->error cross-correlation the impulse response,
		// whose onset lag is the transport delay and peak lag the group delay
		val = data->offset + data->amplitude
			* (2.0 * ((double)rand() / RAND_MAX) - 1.0);
		break;
	}
	if (val > 1) val = 1;
	else if (val < -1) val = -1;
	data->acc += 1;

	switch (data->type) {
	case AYLP_T_VECTOR:
		gsl_vector_set_all(data->vector, val);
		state->vector = data->vector;
		state->header.log_dim.y = data->vector->size;
		state->header.log_dim.x = 1;
		break;
	case AYLP_T_MATRIX:
		gsl_matrix_set_all(data->matrix, val);
		state->matrix = data->matrix;
		state->header.type = self->type_out;
		state->header.log_dim.y = data->matrix->size1;
		state->header.log_dim.x = data->matrix->size2;
		break;
	case AYLP_T_MATRIX_UCHAR:
		gsl_matrix_uchar_set_all(data->matrix_uchar,
			(val+1)/2 * UCHAR_MAX
		);
		state->matrix_uchar = data->matrix_uchar;
		state->header.type = self->type_out;
		state->header.log_dim.y = data->matrix->size1;
		state->header.log_dim.x = data->matrix->size2;
		break;
	}
	state->header.type = self->type_out;
	state->header.units = self->units_out;

	return 0;
}


int test_source_fini(struct aylp_device *self)
{
	struct aylp_test_source_data *data = self->device_data;
	switch (data->type) {
	case AYLP_T_VECTOR:
		xfree_type(gsl_vector, data->vector);
		break;
	case AYLP_T_MATRIX:
		xfree_type(gsl_matrix, data->matrix);
		break;
	}
	xfree(data->freqs);
	xfree(data);
	return 0;
}

