#ifndef AYLP_DEVICES_BIQUAD_FILTER_H_
#define AYLP_DEVICES_BIQUAD_FILTER_H_

#include "anyloop.h"

// which second-order section to realize
enum aylp_biquad_filter_mode {
	AYLP_BIQUAD_FILTER_LOWPASS,
	AYLP_BIQUAD_FILTER_NOTCH,
};

struct aylp_biquad_filter_data {
	// param type in ["vector", "matrix"]
	aylp_type type;
	// units to output
	aylp_units units;
	// filter mode (lowpass or notch)
	enum aylp_biquad_filter_mode mode;
	// design parameters (base set; applies to x/tilt, i.e. elements >= 1)
	double f0;	// center (notch) / corner (lowpass) frequency, Hz
	double q;	// quality factor (sets width/damping)
	double fs;	// assumed sample rate, Hz
	// optional per-axis overrides for element 0 (y/tip); <0 means "use base"
	double f0_y;
	double q_y;
	// normalized Direct Form I coefficients (a0 has been divided out).
	// index [0] applies to element 0 (y/tip); [1] applies to elements >= 1
	// (x/tilt). For matrix input, the base set [1] is used for every element.
	double b0[2], b1[2], b2[2], a1[2], a2[2];
	// per-element filter state: two past inputs and two past outputs
	size_t n;
	double *x1, *x2, *y1, *y2;
	// place to put the result
	union {
		gsl_vector *res_v;
		gsl_matrix *res_m;
	};
};

// initialize biquad_filter device
int biquad_filter_init(struct aylp_device *self);

// process biquad_filter device once per loop
int biquad_filter_proc(struct aylp_device *self, struct aylp_state *state);

// close biquad_filter device when loop exits
int biquad_filter_fini(struct aylp_device *self);

#endif

