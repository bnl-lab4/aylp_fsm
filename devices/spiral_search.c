// Self-contained beam finder: drives an open-loop [y, x] setpoint outward in
// a square spiral, sums the incoming frame's (background-subtracted) pixel
// values as a brightness metric, and re-centres on the brightest point found
// each pass, halving the step and radius for the next pass. Once `passes`
// complete, it parks the output at the best point seen across all passes and
// sets AYLP_DONE, so `anyloop this_config.json` runs the whole search and
// exits with the steering mirror parked at the brightest spot.
//
// Input: the raw camera frame (AYLP_T_MATRIX_UCHAR), same as
// anyloop:center_of_mass expects -- put this device right after the camera
// source. Output: AYLP_T_VECTOR of size 2 (y, x), units `units` param
// (default minmax) -- feed straight into clamp + the DAC device, no
// udp_source/udp_sink or external driver needed.
#include <math.h>
#include <stdbool.h>
#include "anyloop.h"
#include "logging.h"
#include "spiral_search.h"
#include "xalloc.h"


// sum of (pixel - threshold) over the whole frame, clamped to zero per pixel
static double frame_brightness(const gsl_matrix_uchar *mat, unsigned char threshold)
{
	double s = 0.0;
	for (size_t i = 0; i < mat->size1; i++) {
		for (size_t j = 0; j < mat->size2; j++) {
			unsigned char raw = mat->data[i*mat->tda + j];
			unsigned char el = raw > threshold ? raw - threshold : 0;
			s += el;
		}
	}
	return s;
}


// reset the incremental square-spiral generator to its starting point (0, 0).
// (0, 0) itself is emitted by goto_grid_point(data, 0, 0), not by
// spiral_next() -- see the header comment on sg_x/sg_y.
static void spiral_reset(struct aylp_spiral_search_data *data)
{
	data->sg_x = 0;
	data->sg_y = 0;
	data->sg_leg = 1;
	data->sg_direction = 0;	// 0=+x, 1=+y, 2=-x, 3=-y
	data->sg_steps_in_leg = 0;
	data->sg_turns = 0;
}


// advance the spiral generator by one grid point (the next OUTWARD step from
// wherever it last left off), writing it to *out_y/*out_x. Returns false
// (without writing) once the next point would exceed `radius` rings,
// signalling that the current pass is complete.
static bool spiral_next(struct aylp_spiral_search_data *data, int radius,
	int *out_y, int *out_x)
{
	static const int dxs[4] = {1, 0, -1, 0};
	static const int dys[4] = {0, 1, 0, -1};

	int dx = dxs[data->sg_direction];
	int dy = dys[data->sg_direction];
	int nx = data->sg_x + dx;
	int ny = data->sg_y + dy;
	if ((abs(nx) > radius) || (abs(ny) > radius)) return false;

	data->sg_x = nx;
	data->sg_y = ny;
	*out_y = ny;
	*out_x = nx;

	data->sg_steps_in_leg += 1;
	if (data->sg_steps_in_leg == data->sg_leg) {
		data->sg_steps_in_leg = 0;
		data->sg_direction = (data->sg_direction + 1) % 4;
		data->sg_turns += 1;
		if (data->sg_turns % 2 == 0) data->sg_leg += 1;
	}
	return true;
}


static double clamp_minmax(double v)
{
	if (v < -1.0) return -1.0;
	if (v > 1.0) return 1.0;
	return v;
}


// move to a fresh grid point: reset settle countdown and write cur_y/cur_x
static void goto_grid_point(struct aylp_spiral_search_data *data, int gy, int gx)
{
	data->cur_y = clamp_minmax(data->center_y + gy * data->step);
	data->cur_x = clamp_minmax(data->center_x + gx * data->step);
	data->settle_remaining = data->settle_frames;
}


int spiral_search_init(struct aylp_device *self)
{
	self->proc = &spiral_search_proc;
	self->fini = &spiral_search_fini;
	self->device_data = xcalloc(1, sizeof(struct aylp_spiral_search_data));
	struct aylp_spiral_search_data *data = self->device_data;

	// defaults
	data->init_step = 0.05;
	data->init_radius = 10;
	data->passes = 3;
	data->settle_frames = 5;
	data->threshold = 0;
	data->units = AYLP_U_MINMAX;

	if (!self->params) {
		log_error("No params object found.");
		return -1;
	}
	json_object_object_foreach(self->params, key, val) {
		if (key[0] == '_') {
			// keys starting with _ are comments
		} else if (!strcmp(key, "step")) {
			data->init_step = json_object_get_double(val);
			log_trace("step = %G", data->init_step);
		} else if (!strcmp(key, "radius")) {
			data->init_radius = json_object_get_int(val);
			log_trace("radius = %d", data->init_radius);
		} else if (!strcmp(key, "passes")) {
			data->passes = json_object_get_int(val);
			log_trace("passes = %d", data->passes);
		} else if (!strcmp(key, "settle_frames")) {
			data->settle_frames = json_object_get_int(val);
			log_trace("settle_frames = %d", data->settle_frames);
		} else if (!strcmp(key, "threshold")) {
			data->threshold = (unsigned char)json_object_get_int(val);
			log_trace("threshold = %u", data->threshold);
		} else if (!strcmp(key, "units")) {
			const char *s = json_object_get_string(val);
			data->units = aylp_units_from_string(s);
			log_trace("units = %s (0x%hhX)", s, data->units);
		} else {
			log_warn("Unknown parameter \"%s\"", key);
		}
	}
	if (data->init_step <= 0.0) {
		log_error("step must be positive");
		return -1;
	}
	if (data->init_radius < 1) {
		log_error("radius must be >= 1");
		return -1;
	}
	if (data->passes < 1) {
		log_error("passes must be >= 1");
		return -1;
	}
	if (data->settle_frames < 1) {
		log_error("settle_frames must be >= 1");
		return -1;
	}
	if (!data->units) {
		log_error("You must provide a valid units param.");
		return -1;
	}

	data->cmd = xmalloc_type(gsl_vector, 2);
	gsl_vector_set_zero(data->cmd);

	data->step = data->init_step;
	data->radius = data->init_radius;
	data->pass_num = 0;
	data->center_y = data->center_x = 0.0;
	data->pass_best_y = data->pass_best_x = 0.0;
	data->pass_best_b = -1.0;
	data->best_y = data->best_x = 0.0;
	data->best_b = -1.0;
	data->done = false;
	spiral_reset(data);
	goto_grid_point(data, 0, 0);

	log_info("spiral_search: starting (step=%G radius=%d passes=%d "
		"settle_frames=%d)", data->init_step, data->init_radius,
		data->passes, data->settle_frames);

	self->type_in = AYLP_T_MATRIX_UCHAR;
	self->units_in = AYLP_U_ANY;
	self->type_out = AYLP_T_VECTOR;
	self->units_out = data->units;
	return 0;
}


int spiral_search_proc(struct aylp_device *self, struct aylp_state *state)
{
	struct aylp_spiral_search_data *data = self->device_data;

	if (!data->done) {
		double brightness = frame_brightness(state->matrix_uchar,
			data->threshold);

		if (data->settle_remaining > 0) {
			data->settle_remaining -= 1;
		} else {
			// this frame reflects cur_y/cur_x -- trust it
			if (brightness > data->pass_best_b) {
				data->pass_best_b = brightness;
				data->pass_best_y = data->cur_y;
				data->pass_best_x = data->cur_x;
			}

			int gy, gx;
			if (spiral_next(data, data->radius, &gy, &gx)) {
				goto_grid_point(data, gy, gx);
			} else {
				// pass complete: recentre on this pass's best
				log_info("spiral_search: pass %d/%d brightest "
					"y=%+.4f x=%+.4f (brightness=%.0f)",
					data->pass_num + 1, data->passes,
					data->pass_best_y, data->pass_best_x,
					data->pass_best_b);
				data->center_y = data->pass_best_y;
				data->center_x = data->pass_best_x;
				if (data->pass_best_b > data->best_b) {
					data->best_b = data->pass_best_b;
					data->best_y = data->pass_best_y;
					data->best_x = data->pass_best_x;
				}
				data->pass_num += 1;
				if (data->pass_num >= data->passes) {
					data->done = true;
					data->cur_y = data->best_y;
					data->cur_x = data->best_x;
					log_info("spiral_search: settled at "
						"y=%+.4f x=%+.4f "
						"(brightness=%.0f)",
						data->best_y, data->best_x,
						data->best_b);
				} else {
					data->step /= 2.0;
					data->radius = data->radius / 2 > 0
						? data->radius / 2 : 1;
					data->pass_best_b = -1.0;
					spiral_reset(data);
					goto_grid_point(data, 0, 0);
				}
			}
		}
	}

	gsl_vector_set(data->cmd, 0, data->cur_y);
	gsl_vector_set(data->cmd, 1, data->cur_x);
	state->vector = data->cmd;
	state->header.type = self->type_out;
	state->header.units = self->units_out;
	state->header.log_dim.y = 2;
	state->header.log_dim.x = 1;
	if (data->done) aylp_set_done(&state->header);
	return 0;
}


int spiral_search_fini(struct aylp_device *self)
{
	struct aylp_spiral_search_data *data = self->device_data;
	xfree_type(gsl_vector, data->cmd);
	xfree(data);
	return 0;
}
