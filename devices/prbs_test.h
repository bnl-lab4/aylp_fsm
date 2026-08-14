#ifndef AYLP_DEVICES_PRBS_TEST_H_
#define AYLP_DEVICES_PRBS_TEST_H_

#include <stdbool.h>
#include <stdio.h>

#include "anyloop.h"

// test stages; see prbs_test.c
#define PT_WARMUP	0	// park at bias; sensor side acquires
#define PT_BURST	1	// driving the PRBS
#define PT_QUIET	2	// parked at bias between bursts
#define PT_DONE		3

struct aylp_prbs_test_data {
	// params
	int index_cmd;		// output element to drive
	int index_err;		// input element to watch
	size_t out_size;	// output vector length
	double amplitude;	// PRBS half-swing, command units
	double bias;		// command between bursts and after the test
	unsigned order;		// LFSR order; the run is 2^order - 1 chips
	size_t chip_frames;	// loop iterations per chip
	size_t n_bursts;	// bursts to average
	size_t quiet_frames;	// frames parked at bias after each burst
	double warmup;		// seconds parked at bias before the first burst
	size_t max_lag;		// positive lags evaluated
	size_t neg_lags;	// acausal lags, which measure the noise floor
	double onset_frac;	// fraction of the peak that counts as onset
	double phase_f_lo;	// band for the cross-spectrum phase-slope fit
	double phase_f_hi;
	double phase_notch_hz;	// optional periodic-light fundamental to exclude
	double phase_notch_width_hz;
	bool use_correlation_delay;
	double volts_per_unit;	// DAC scale of the swept channel, for reporting
	double pixel_scale;	// px per response unit, for reporting
	bool use_rejected;	// correlate held samples too, as if they were data
	double min_pair_frac;	// fraction of the window a lag needs to be fit
	// Optional final-result quality gates. Zero disables a gate. These are
	// deliberately separate from the 4-sigma peak detector: finding a lobe
	// is not the same thing as proving a controller-quality delay.
	double min_peak_snr;
	double max_phase_resid_deg;
	size_t min_phase_bins;
	double min_burst_frac;
	double max_peak_mad;
	double min_live_frac;
	size_t max_holes;
	double max_phase_peak_delta;
	char *output_file;	// PDF path; the .dat is written alongside it
	char *config;		// free-text note copied into the .dat header

	// output vector; element index_cmd carries the PRBS, rest are 0
	gsl_vector *out;

	// stimulus
	size_t n_chips;		// 2^order - 1
	signed char *chips;	// +-1, one per chip
	size_t burst_frames;	// n_chips * chip_frames
	size_t win;		// burst_frames + quiet_frames
	double *cmd_win;	// the command of one whole window, mean removed

	// state machine
	int stage;
	double t0;		// time of the first proc (s, CLOCK_MONOTONIC)
	size_t burst;		// bursts started so far
	size_t frame;		// frames into the current window
	double cmd;		// command currently published

	// per-burst response, [burst*win + frame]
	double *resp;
	// 1 where the sensor published a live sample for that frame, 0 where it
	// flagged the frame and held its previous output. Tracked frame by
	// frame; nothing here assumes a duty cycle or a chop period.
	unsigned char *live;

	// loop-rate bookkeeping
	double t_prev, dt_sum, dt_sum2;
	size_t dt_cnt;

	// results
	double *rho;		// correlogram over [-neg_lags, max_lag]
	double *rbar;		// ensemble-mean response over the window
	double *rsem;		// standard error of that mean
	size_t *rcnt;		// live bursts that went into rbar[k]
	unsigned char *rlive;	// 1 where rbar[k] has at least one live burst

	// how much of the run the sensor actually delivered
	size_t n_frames;	// frames recorded into a burst window
	size_t n_live;		// of those, frames the sensor called good
	size_t n_holes;		// window positions with no live sample at all
	size_t worst_cnt;	// fewest live bursts at any window position
	size_t n_thin_lags;	// lags dropped for want of live pairs
	double fs;		// measured loop rate, Hz
	double rho_peak;	// signed peak of the correlogram
	double lag_peak;	// parabola-refined peak lag, frames
	double lag_onset;	// interpolated onset lag, frames
	double lag_centroid;	// first moment of the correlation lobe
	double rho_noise;	// rms correlation over the acausal lags
	double tau_phase_ms;	// delay from the cross-spectrum phase slope
	double phase_resid_deg;	// rms residual of that fit
	size_t phase_n;		// points it used
	double lag_peak_med;	// median over the per-burst estimates
	double lag_peak_mad;	// and their median absolute deviation
	double lag_onset_med;
	double lag_onset_mad;
	size_t n_burst_ok;
	bool have_result;
	bool quality_ok;
	char quality_reason[512];
};

int prbs_test_init(struct aylp_device *self);
int prbs_test_proc(struct aylp_device *self, struct aylp_state *state);
int prbs_test_fini(struct aylp_device *self);

#endif
