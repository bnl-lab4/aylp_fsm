#ifndef AYLP_DEVICES_SPIRAL_SEARCH_H_
#define AYLP_DEVICES_SPIRAL_SEARCH_H_

#include <stdbool.h>
#include <gsl/gsl_vector.h>
#include "anyloop.h"

struct aylp_spiral_search_data {
	// output command vector [y, x], normalised units
	gsl_vector *cmd;

	// params
	double init_step;	// initial grid step, normalised units
	int init_radius;	// initial grid radius, rings
	int passes;		// number of coarse-to-fine passes
	int settle_frames;	// frames to wait after a move before trusting
				// the brightness reading (transport delay:
				// camera exposure + mirror response)
	unsigned char threshold;	// background subtracted before summing
	aylp_units units;	// output units (default minmax)

	// current pass state
	double step;		// current pass's grid step (halved each pass)
	int radius;		// current pass's grid radius (halved each pass)
	int pass_num;		// 0-indexed current pass
	int settle_remaining;	// frames left before the current point is
				// trustworthy
	double center_y, center_x;	// current pass's spiral centre
	double cur_y, cur_x;		// currently commanded point

	// per-pass and overall best-seen point
	double pass_best_y, pass_best_x, pass_best_b;
	double best_y, best_x, best_b;

	// true once all passes are done and cmd is parked at best_y/best_x
	bool done;

	// incremental square-spiral generator state (see spiral_next()). sg_x/
	// sg_y track the last point *emitted* -- (0, 0) is emitted directly by
	// goto_grid_point() at init and at the start of each pass, never by
	// spiral_next(), so spiral_next() only ever computes the next OUTWARD
	// step from here.
	int sg_x, sg_y;
	int sg_leg, sg_direction, sg_steps_in_leg, sg_turns;
};

// initialize spiral_search device
int spiral_search_init(struct aylp_device *self);

// process spiral_search device once per loop
int spiral_search_proc(struct aylp_device *self, struct aylp_state *state);

// close spiral_search device when loop exits
int spiral_search_fini(struct aylp_device *self);

#endif
