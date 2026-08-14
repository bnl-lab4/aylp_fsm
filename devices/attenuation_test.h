#ifndef AYLP_DEVICES_ATTENUATION_TEST_H_
#define AYLP_DEVICES_ATTENUATION_TEST_H_

#include <time.h>

#include "anyloop.h"

// measurement phases, in order
enum {
	ATTEN_PHASE_START = 0,	// start_delay: hold loop open, don't record
	ATTEN_PHASE_OPEN,	// hold loop open, record error
	ATTEN_PHASE_SETTLE,	// loop closed, let it acquire, don't record
	ATTEN_PHASE_CLOSED,	// loop closed, record error
	ATTEN_PHASE_DONE,	// analysis written; AYLP_DONE raised
};

struct aylp_attenuation_test_data {
	// parameters
	double start_delay;	// s to wait before open recording (default 0)
	double open_time;	// s of open-loop recording (required)
	double settle_time;	// s after closing before recording (default 5)
	double closed_time;	// s of closed-loop recording (required)
	size_t nfft;		// Welch segment length, power of 2 (default 4096)
	char *output_file;	// PDF filename (default "attenuation.pdf")
	char *labels;		// comma-separated element labels (optional)
	double pixel_scale;	// px per output unit for the report (default 1)
	char *config;		// config summary printed on the PDF (optional)
	// Pass the measured open-loop error downstream instead of replacing it
	// with zero. Safe only when the controller has a matching startup hold;
	// used to identify predictive controllers while their command stays zero.
	bool pass_open;

	// state
	unsigned phase;
	int analyzed;		// analysis has run (success or not)
	double t0;		// timestamp of first proc call
	time_t start_wall;	// wall-clock time of first proc call
	size_t n_elem;		// vector length, latched on first proc call

	// recorded error, interleaved [sample*n_elem + elem]
	double *open_buf;
	size_t open_n, open_cap;		// in samples
	double open_t_first, open_t_last;
	double *closed_buf;
	size_t closed_n, closed_cap;		// in samples
	double closed_t_first, closed_t_last;
};

int attenuation_test_init(struct aylp_device *self);
int attenuation_test_proc(struct aylp_device *self, struct aylp_state *state);
int attenuation_test_fini(struct aylp_device *self);

#endif
