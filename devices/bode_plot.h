#ifndef AYLP_DEVICES_BODE_PLOT_H_
#define AYLP_DEVICES_BODE_PLOT_H_

#include <stdbool.h>
#include "anyloop.h"

// output modes for bode_plot
enum {
	BODE_OUTPUT_RAW  = 0,  // use element value directly as y
	BODE_OUTPUT_DIFF = 1,  // use element - element_previous as y (e.g. com-com_prev)
};

// inject modes for bode_plot
enum {
	BODE_INJECT_ADD     = 0,  // add u to existing element (closed-loop perturbation)
	BODE_INJECT_REPLACE = 1,  // replace element with u (open-loop FSM command)
};

#define BODE_MAX_SEGMENTS 32
#define BODE_MAX_AVOID 32

struct aylp_bode_plot_data {
	// required
	double sample_rate;     // Hz (nominal; the measured rate is recorded too)
	// frequency sweep range
	double freq_start;      // Hz (default: 1.0)
	double freq_end;        // Hz (default: Nyquist)
	size_t n_freqs;         // number of log-spaced frequencies (default: 30)
	// measurement quality
	size_t n_settle;        // samples to skip at each frequency (default: 50)
	double n_settle_cycles; // ALSO wait this many cycles before integrating
				// (default 2.0) -- 50 samples is 1/15 of a period
				// at 5 Hz, so a sample-only settle leaks the
				// previous frequency's transient into the point
	double n_cycles;        // sinusoidal cycles to average per frequency (default 10)
	size_t n_segments;      // split each point into this many integer-cycle
				// segments (default 7); the point estimate is the
				// time-weighted coherent mean of the segment
				// phasors and the segment scatter gives per-point
				// error bars + a quality (coherence) metric
	// injection
	double amplitude;       // injection amplitude in pipeline units (default: 0.1)
	double amplitude2;      // optional second full sweep at this amplitude for a
				// linearity/clipping check (default 0 = off)
	unsigned inject_mode;   // BODE_INJECT_ADD or BODE_INJECT_REPLACE
	double start_delay;     // seconds to inject 0 before the sweep starts (default: 0)
	// grid hygiene: nudge sweep frequencies off known disturbance lines --
	// a line within the lock-in's leakage skirt biases the point
	double avoid_freqs[BODE_MAX_AVOID];
	size_t n_avoid;
	double avoid_df;        // clearance from each avoided line, Hz (default 0.75)
	// fit band for the built-in flat-gain + pure-delay fit
	double k_max_freq;      // K = error-weighted mean |H| below this (default 30 Hz)
	// output channel
	unsigned output_mode;   // BODE_OUTPUT_RAW or BODE_OUTPUT_DIFF
	size_t element;         // linear index into vector/matrix (default: 0)
	// frequency spacing
	bool freq_linear;       // if true, linear spacing; default false (log)
	// output
	char *output_file;      // PDF filename (default: "bode_plot.pdf")
	char *config_note;      // free text recorded in the .dat header (operating
				// point: biases, alignment, camera config, ...)

	// sweep state
	double delay_start;     // timestamp of first frame (start_delay reference); 0 until seen
	bool delay_done;        // start_delay has elapsed; sweep is running
	size_t freq_idx;        // current frequency index
	size_t pass;            // 0 = main sweep, 1 = amplitude2 linearity sweep
	size_t step_count;      // samples elapsed at current frequency
	double start_time;      // real-world start timestamp of current frequency
	double prev_time;       // real-world timestamp of previous iteration
	double phase;           // current excitation phase accumulator
	double prev_value;      // previous raw value for diff mode
	double prev_cos_term;   // trapezoidal integration baseline (cos)
	double prev_sin_term;   // trapezoidal integration baseline (sin)
	double prev_phi;        // injection phase from previous step (y[n] responds to phi[n-1])
	// integer-cycle segment accumulators: each segment integrates exactly
	// cycles_per_seg cycles (boundaries on excitation-phase crossings, so
	// no partial-cycle DC/line leakage), then contributes one phasor
	bool integrating;       // settle finished; segment accumulation running
	double cycles_per_seg;  // integer cycles per segment
	double seg_start_phi;   // corr_phi at the current segment's start
	double seg_Yc, seg_Ys;  // current segment correlator integrals
	double seg_T;           // current segment duration
	size_t seg_idx;         // completed segments at this frequency
	double seg_re[BODE_MAX_SEGMENTS];   // per-segment H (re/im) and duration
	double seg_im[BODE_MAX_SEGMENTS];
	double seg_dur[BODE_MAX_SEGMENTS];
	// measured frame rate across the whole sweep
	double dt_sum;
	size_t dt_n;

	// hardware FG mode (optional — bypasses software sine injection)
	// Set fg_port to use the DAQC2 hardware function generator instead of computing
	// sine in software.  The FG runs at 200 kHz so the output waveform is clean
	// regardless of loop-rate jitter.  The piplate_bridge for the swept axis must
	// NOT be in the pipeline when fg_port is set (they would fight over the port).
	char  *fg_port;    // serial device, e.g. "/dev/ttyACM0"; NULL = software mode
	int    fg_fd;      // open serial fd, -1 when unused
	int    fg_board;   // DAQC2 board 0-7 (default 0)
	int    fg_channel; // FG output channel 1 or 2 (default 1)
	int    fg_level;   // amplitude: 1=eighth 2=quarter 3=half 4=full (default 1)
	double fg_offset;  // DC center voltage in V, 0.0-4.095 (default 0 = hardware default 2.048V)

	// results (one entry per frequency)
	double *freqs;
	double *mags;
	double *phases;
	double *mag_sem;   // standard error of the mean over segments (magnitude)
	double *ph_sem;    // standard error of the mean over segments (phase, rad)
	double *quality;   // |coherent mean| / mean |segment phasor|; 1 = clean
	double *mags2;     // amplitude2 pass results (if enabled)
	double *phases2;
	// flat-gain + pure-delay fit results (computed at end of sweep)
	double fit_K;       // error-weighted mean |H| for f <= k_max_freq
	double fit_K_err;
	double fit_tau_ms;  // from the weighted phase-slope fit (EXCLUDES the
			    // hidden frame; see the .dat header note)
	double fit_frames;  // loop delay in frames at the MEASURED rate,
			    // INCLUDING the +1 hidden frame
	double fit_mag_flat_rms;  // rms fractional |H| deviation from fit_K, full band
	double fit_ph_resid_deg;  // rms phase-fit residual, degrees
};

int bode_plot_init(struct aylp_device *self);
int bode_plot_proc(struct aylp_device *self, struct aylp_state *state);
int bode_plot_fini(struct aylp_device *self);

#endif
