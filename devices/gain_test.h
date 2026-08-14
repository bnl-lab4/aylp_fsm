#ifndef AYLP_DEVICES_GAIN_TEST_H_
#define AYLP_DEVICES_GAIN_TEST_H_

#include <stdbool.h>
#include <stdio.h>

#include "anyloop.h"

// test stages; see gain_test.c
#define GT_WARMUP	0	// park at bias; sensor side acquires
#define GT_RAMP_IN	1	// glide bias -> first level
#define GT_SWEEP	2	// the staircase itself
#define GT_RAMP_OUT	3	// glide last level -> bias
#define GT_DONE		4

// which branch of an up/down sweep a level belongs to
#define GT_BRANCH_UP	0
#define GT_BRANCH_DOWN	1

struct aylp_gain_test_data {
	// params
	int index_cmd;		// output element to drive
	int index_err;		// input element to watch
	size_t out_size;	// output vector length
	double low, high;	// staircase endpoints, command units
	double step;		// staircase increment, command units
	double bias;		// command held before and after the sweep
	double dwell;		// seconds spent at each level
	double settle_frac;	// trailing fraction of a dwell that is averaged
	double warmup;		// seconds parked at bias before the sweep
	double ramp;		// seconds to glide into/out of the staircase
	bool updown;		// sweep back down again after reaching `high`
	bool interleave;		// alternate low/high inward to decorrelate drift
	size_t cycles;		// balanced interleaved passes over every level
	bool use_rejected;	// include tracker-held samples (diagnostic only)
	double volts_per_unit;	// DAC scale of the swept channel, for reporting
	double pixel_scale;	// px per response unit, for reporting
	// "pixel_scale": "auto" -- resolve it at init from the frame size the
	// source publishes (libaylp/timing.h) as (dim - 1) / 2 on the measured
	// axis, matching how the centroider normalizes. Necessary once the
	// source sizes its own ROI: a hardcoded scale is then a number nobody
	// can keep correct, and it silently rescales every px column in the
	// report.
	bool pixel_scale_auto;
	double resp_max;	// |response| past this means the beam is leaving
	double fit_range;	// extra fit over |cmd - bias| <= this (0 = off)
	double min_r2;		// refuse report below this fit quality (0 = off)
	char *output_file;	// PDF path; the .dat is written alongside it
	char *config;		// free-text note copied into the .dat header
	// Auto-write: install the measured gain straight into a run config's fsp
	// stage, so a gain run leaves the controller ready instead of leaving a
	// number to be copied by hand. Only ever writes the SMALL-SIGNAL slope,
	// and only when the fit passes the same gates the calibration suite
	// applies (see gt_autowrite). Null = don't write.
	char *write_config;	// path to the run config to update
	int write_axis;		// 0 = y, 1 = x; defaults to index_err

	// output vector; element index_cmd carries the staircase, rest are 0
	gsl_vector *out;

	// state machine
	int stage;
	double t0;		// time of the first proc (s, CLOCK_MONOTONIC)
	double stage_t0;	// time the current stage/level began
	double cmd;		// command currently published
	double ramp_from;	// command the active ramp started at

	// staircase position
	size_t n_levels;	// levels in one branch
	size_t n_points;	// levels in the whole sweep (2x if updown)
	size_t point;		// index of the level being measured

	// per-level settled-window accumulators (Welford, so a level that
	// barely moves doesn't lose its variance to cancellation)
	double wmean, wm2;
	size_t wcnt;
	double wfirst;		// first settled sample of the level
	bool wmoved;		// any settled sample differed from it

	// per-level results, n_points long
	double *cmds;		// commanded value, command units
	double *resp;		// settled mean response
	double *resp_sd;	// within-window standard deviation
	double *resp_sem;	// standard error of the settled mean
	double *times;		// seconds since t0 at the middle of the window
	size_t *counts;		// samples in the settled window
	unsigned char *valid;	// level is usable in the fit
	unsigned char *branch;	// GT_BRANCH_UP / GT_BRANCH_DOWN
	size_t n_done;		// levels finalized so far
	size_t n_bad_run;	// consecutive unusable levels (2 aborts the sweep)
	bool seen_valid;	// at least one level has been usable

	// loop-rate bookkeeping
	double t_prev, dt_sum;
	size_t dt_cnt;

	// fit results
	double slope, slope_err;	// response units per command unit
	double intercept;
	double r2;
	double resid_rms;
	double nonlin;			// max |residual|, response units
	double span;			// response span over the valid levels
	double slope_drift;		// slope of a fit that also carries a
	double drift_rate;		// linear-in-time term (units, units/s)
	double slope_small;		// fit restricted to fit_range
	double slope_small_err;
	size_t n_small;
	double hysteresis;		// mean down-branch minus up-branch
	bool have_fit;
};

int gain_test_init(struct aylp_device *self);
int gain_test_proc(struct aylp_device *self, struct aylp_state *state);
int gain_test_fini(struct aylp_device *self);

#endif
