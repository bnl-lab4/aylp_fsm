#ifndef AYLP_DEVICES_LATENCY_TEST_H_
#define AYLP_DEVICES_LATENCY_TEST_H_

#include <stdbool.h>
#include "anyloop.h"

// test stages; see latency_test.c
#define LT_WARMUP	0
#define LT_CAL		1
#define LT_MEAS		2
#define LT_DONE		3

struct aylp_latency_test_data {
	// params
	int index_cmd;		// output element to drive (0 = x_fine)
	int index_err;		// input element to watch (1 = x for CoM [y, x])
	size_t out_size;	// output vector length (4 for the dual-stage DAC map)
	double low, high;	// command levels, minmax units
	double period;		// seconds spent at each level (one half-cycle)
	double warmup;		// seconds to hold `low` before starting (CoM acquire)
	double settle_frac;	// last fraction of a half-cycle treated as settled
	size_t cal_cycles;	// full low/high cycles used for calibration
	size_t n_steps;		// command edges to measure before reporting
	char *results_file;	// append per-edge + summary results here (optional)
	char *label;		// free-text run label for the results header
	FILE *rf;		// open results file, or 0

	// output vector; element index_cmd carries the square wave, rest are 0
	gsl_vector *out;

	// state machine
	int stage;
	double t0;		// time of the first proc (s, CLOCK_MONOTONIC)
	double half_t0;		// time the current half-cycle's command was emitted
	bool level_high;	// current command level
	size_t halfcycles;	// half-cycles completed in the current stage

	// settled-window accumulators for the current half-cycle
	double wsum, wsum2;
	size_t wcnt;

	// calibration results. The step is estimated from consecutive
	// high/low window-mean differences (cancels slow beam drift) and the
	// noise from within-window scatter only (so drift doesn't count as
	// noise).
	double win_mean_prev;	// settled mean of the previous half-cycle
	bool win_prev_valid;
	double step_sum;	// sum of consecutive (high - low) window means
	size_t step_pairs;
	double var_sum;		// within-window variance, weighted by count
	size_t var_cnt;
	double msum[2];		// per-level window-mean sums, [0]=low [1]=high
	size_t mcnt[2];
	double step;		// signed settled response to a low->high step
	double sigma;		// fast (within-window) noise, standard deviation
	double gain;		// plant DC gain = step / (high - low), i.e. settled
				// error response (px) per unit minmax command.
				// Same units as bode_plot's fit_K, so the latency
				// run and a bode sweep are directly comparable.

	// baseline for the edge being measured: the settled mean of the
	// half-cycle we just left, re-measured every edge so beam drift
	// can't stale it
	double from_level;

	// per-edge measurement
	size_t edges_done;	// command edges finalized
	size_t edges_missed;	// edges where a crossing was never seen
	double *dep, *half;	// latencies (s): first movement / 50% crossing
	size_t n_dep, n_half;
	double dep_pend;	// time of a first qualifying sample (<0 = none);
				// movement needs 2 consecutive samples past the
				// threshold so single noise spikes don't count
	bool dep_found, half_found;

	// loop-rate bookkeeping (sample interval during measurement)
	double t_prev, dt_sum;
	size_t dt_cnt;
};

int latency_test_init(struct aylp_device *self);
int latency_test_proc(struct aylp_device *self, struct aylp_state *state);
int latency_test_fini(struct aylp_device *self);

#endif
