#ifndef AYLP_DEVICES_FSP_H_
#define AYLP_DEVICES_FSP_H_

#include <stdbool.h>
#include <time.h>

#include "anyloop.h"

enum aylp_fsp_transient_controller {
	AYLP_FSP_TRANSIENT_PROPORTIONAL = 0,
	AYLP_FSP_TRANSIENT_INTEGRAL,
	AYLP_FSP_TRANSIENT_PI,
	AYLP_FSP_TRANSIENT_MODAL,
	AYLP_FSP_TRANSIENT_HYBRID,
};

enum aylp_fsp_transient_mode {
	AYLP_FSP_TRANSIENT_MANUAL = 0,
	AYLP_FSP_TRANSIENT_OFF,
	AYLP_FSP_TRANSIENT_AUTO,
};

// Adaptive Filtered Smith Predictor / adaptive LQG for tip-tilt beam
// stabilization. This is a COMPLETE controller: it takes the [y, x] error
// (center-of-mass output) and emits the [y, x] command, so it REPLACES both
// anyloop:kalman_filter and anyloop:pid in the pipeline (running it alongside
// either double-compensates the loop delay and will pump lines). It sits where
// pid used to: after center_of_mass and any logging sink, before clamp.
//
// Structure (per axis, independent; index 0 = y, 1 = x, matching pid's
// [y,x] convention). The plant seen by the loop is a gain K, an optional
// stable/minimum-phase frequency-shaping biquad H(z), and a transport delay
// of `delay + delay_frac` samples (camera + compute + DAC ZOH):
//
//     e(k) = phi(k) + K H(z) D(z) u(k)
//
// where phi is the disturbance (bench vibration) the beam would see open-loop,
// and u is our command. The three pieces, following Kulcsar 2006 / Petit 2008
// / Meimon 2010 / Correia 2010 (optimal/adaptive control for AO):
//
//  1. SMITH-PREDICTOR CORE. We know K, H(z), and delay, so we reconstruct the
//     disturbance from the measured error by removing our own modeled plant
//     contribution: phi_meas = e - K H(z) D(z) u. This internal model takes
//     the identified plant out of the disturbance estimate. During the
//     open-loop startup hold (command forced to 0) phi_meas is exactly the raw
//     open-loop disturbance, which is the cleanest identification data.
//
//  2. DISTURBANCE MODEL + KALMAN FILTER. phi is modeled as a sum of `n_modes`
//     narrowband AR(2) vibration modes (Meimon 2010's tip-tilt model), each a
//     damped oscillator driven by white noise:
//         x_i(k) = a1_i x_i(k-1) + a2_i x_i(k-2) + nu_i,   var(nu_i)=q_i
//         a1_i =  2 exp(-2 pi z_i f_i Ts) cos(2 pi f_i Ts sqrt(1 - z_i^2))
//         a2_i = -exp(-4 pi z_i f_i Ts)
//     The stacked state X = [x_i(k), x_i(k-1)]_i is estimated from phi_meas by
//     a Kalman filter with a STATIONARY gain L, precomputed at init by
//     iterating the Riccati recursion to convergence. The q_i/r ratio is the
//     LQG tuning knob (process vs measurement noise): larger q_i = more
//     aggressive tracking of that line.
//
//  3. delay-STEP PREDICTION + MINIMUM-VARIANCE CONTROL. The command issued now
//     lands `delay` samples from now, so we project the modal state forward by
//     running the AR recursion `delay` steps (the mean prediction, zero future
//     noise) and cancel the predicted disturbance:
//         phi_hat(k+delay|k) = C A^delay Xhat(k|k)
//         u(k) = H^-1(z) [-phi_hat(k+delay|k) / K]
//     This is what "gets around" the delay on the PREDICTABLE (modal) part of
//     phi; the unpredictable broadband floor and anything above the loop's
//     achievable bandwidth are NOT helped (a fundamental limit, not an
//     implementation one).
//
// ADAPTATION (Petit/Meimon): the vibration lines wander (~0.5 Hz/run on this
// bench) and vary several-fold in amplitude within hours, so a fixed model
// detunes. Every `adapt_period` seconds we re-identify from buffered data:
// update each mode's process-noise q_i from its recent innovation energy,
// update the measurement noise r from the broadband innovation floor,
// optionally nudge each mode center frequency f_i toward the locally measured
// line (bounded by `adapt_df_max`), rebuild the AR coefficients, and recompute
// the stationary gain L. Set adapt_period <= 0 for a fixed (non-adaptive) FSP.
//
// COMMAND ROBUSTNESS FILTER (the "F" in FSP; cmd_fc > 0): the raw
// minimum-variance command injects the estimate's broadband/noise content onto
// the beam at frequencies where prediction is worthless, which is what showed
// up as a ~3x waterbed at 120-400 Hz in the fsp_sim.py study. A 2nd-order
// low-pass at cmd_fc rolls that authority off. A plain filter would also lag
// the in-band line cancellation, so each mode is PRE-COMPENSATED at its
// center frequency: mode i's AR state is rolled forward to the real-valued
// horizon delay + delay_frac + n_i, n_i = -arg H(f_i) / omega_i, so the
// filter's phase lag at that line nets out to ~0 (the fractional remainder
// is blended between the two adjacent integer steps, same as the plain
// delay_frac case, rather than rounding n_i away), and its contribution is
// boosted by min(1/|H(f_i)|, 3) to undo the gain droop. (Do NOT "fix" this
// back to a quadrature rotation of
// the state: extracting the quadrature divides by sin(omega_i), which
// amplifies the white part of the Kalman state ~1/sin(omega_i) -- ~65x for a
// 5.6 Hz mode at 2310 Hz -- straight into the command; the fsp_sim.py study
// measured the closed loop getting WORSE that way. Rolling through the stable
// AR dynamics is a contraction: same phase advance, no noise amplification.)
// Recomputed on every adaptation tick when f_i moves. cmd_fc = 0 disables
// (raw minimum-variance command).
//
// BURST GUARD (guard_ratio > 0, on by default): any mismatch between the
// configured plant model (K, delay) and the true plant leaks the command back
// into the Smith reconstruction, closing a parasitic loop whose phase hits
// -180 deg at f = fs / (2 * (delay + delay_frac)) -- ~310-340 Hz here. When
// the margin at that frequency goes negative (K drifts with the coarse bias /
// alignment; see the 2026-07-16 "310 Hz spiral" and the 2026-07-17 run-15
// rerun), the loop emits intermittent BURSTS of ringing there that the 500 s
// PSD averages into a broad 250-450 Hz hump and a ruined RMS. The guard runs
// a per-axis band-pass at that frequency on the raw error, tracks a fast
// envelope against a slow baseline (with an absolute floor so a quiet bench
// can't lower the bar to noise), and on detection: cuts that axis's authority
// to zero, freezes NLMS training and the adaptation statistics (so the
// predictor does not learn the ringing), holds guard_hold seconds, then ramps
// authority back over guard_ramp seconds. Every activation is counted and a
// ticker line is printed every guard_tick seconds of closed-loop time. This
// keeps a marginal run alive; the FIX for recurring activations is
// re-measuring K/delay at the current operating point.
//
// STALL-GAP HANDLING (gap_trip > 0, on by default): if the source drops
// frames (a scheduler hiccup) or stalls outright (the ASI camera stream
// stall + capture restart, a ~0.3-0.6 s hole), the DAC holds the last
// command while the disturbance keeps moving, and every history in the
// predictor ends up spanning the hole. Neither re-entering at full
// authority onto stale state (burst seed) nor cutting authority (releases
// the DC correction: the FSM snaps to bias and the accumulated drift
// re-appears as a multi-px excursion -- observed 2026-07-17) is right.
// Instead the gap is simply PATCHED, sized to the frames actually missed:
// the Smith command ring is padded with the held command (exact -- that is
// what the plant received), the modal state is propagated through the gap
// with the AR recursion (the model knows how phases advance), and the NLMS
// input history is padded with a slow EWMA of the reconstructed disturbance
// -- so the prediction carries the DC correction straight through the hole
// and recovers full AC prediction as real frames refill the tap window
// (~broad_hist_len/fs, ~140 ms). No hold, no ramp, no authority change: a
// 2-frame drop costs 2 frames; a capture restart costs the hole plus the
// refill. NLMS training pauses only while fabricated samples meaningfully
// pollute the tap window (skipped entirely for drops of a few frames).
// Gap events are counted separately from burst events: bursts mean
// "re-measure K", gaps mean "the camera stalled".

#define AYLP_FSP_MAX_MODES 8
#define AYLP_FSP_MAX_DIM (2 * AYLP_FSP_MAX_MODES)

// per-axis controller state
struct aylp_fsp_axis {
	size_t n_modes;
	size_t dim;			// = 2 * n_modes
	double K;			// signed plant gain (error units per command
					// unit); sign must make the loop negative
					// feedback -- verify with a push test
	// Optional Bode-fitted plant shape
	//
	//               b0 + b1 z^-1 + b2 z^-2
	//     H(z) = --------------------------------
	//               1  + a1 z^-1 + a2 z^-2
	//
	// plant_b/plant_a default to identity. Both H and H^-1 must be stable:
	// the Smith reconstruction filters the delayed command through H, while
	// the cancellation command is filtered through H^-1 so K H u equals the
	// requested minimum-variance correction. This matched pair changes the
	// plant model without changing the disturbance predictor.
	bool plant_shaped;
	double plant_b[3];
	double plant_a[3];
	double plant_ib[3];		// normalized numerator of H^-1
	double plant_ia[3];		// normalized denominator of H^-1
	double plant_z1, plant_z2;	// delayed-command H(z) state
	double plant_iz1, plant_iz2;	// requested-command H^-1(z) state
	// Per-axis transport delay in samples. The two axes' plants differ
	// (bode 2026-07-16: x 5.62 vs y 6.37 frames at 3788 Hz), so each axis
	// may set its own "delay"/"delay_frac" inside its object; unset values
	// inherit the global ones. Resolved at init (0 / <0 = inherit).
	size_t delay;
	double delay_frac;
	// This axis takes its delay from the startup auto estimate (see
	// delay_auto in aylp_fsp_data).
	bool delay_auto;
	// The delay every length below was SIZED for. Equal to `delay`
	// normally; with an auto delay it is delay_auto_max, because the
	// estimate does not exist until the loop has been running and the rings
	// cannot be reallocated underneath it.
	size_t delay_alloc;
	// Per-axis output bounds, u held to [clamp_lo, clamp_hi]. Needed
	// separately from the global pair because the axes can have OPPOSITE
	// command->voltage signs (a steering pair typically does: +1 V/unit on
	// one, -1 V/unit on the other). An asymmetric window is mirrored by a
	// negative scale, so one global window cannot express the same VOLTAGE
	// limit on both axes -- [-12.5, +7.5] means -10..+10 V through a +1
	// scale and +15..-5 V through a -1 scale. NAN = inherit the global.
	double clamp_lo, clamp_hi;
	// per-axis command ring and full-band observer bookkeeping -- lengths
	// depend on the axis delay so they cannot be shared
	size_t uhead;			// ring write index into ucmd
	size_t broad_hist_len;		// = broad_order + delay + broad_gd + 2
	size_t broad_head;
	size_t broad_seen;
	// per-mode nominal parameters (also the adaptation targets)
	double f[AYLP_FSP_MAX_MODES];		// center frequency (Hz)
	double zeta[AYLP_FSP_MAX_MODES];	// damping ratio (0..1)
	double q[AYLP_FSP_MAX_MODES];		// process noise variance
	double r;				// measurement noise variance
	// derived AR(2) coefficients per mode
	double a1[AYLP_FSP_MAX_MODES];
	double a2[AYLP_FSP_MAX_MODES];
	// stationary variance gain of each mode: var(state)/var(drive) =
	// (1-a2)/((1+a2)(1-a1-a2)(1+a1-a2)); ~1/(4 pi zeta f Ts) for light
	// damping, i.e. HUGE for sharp lines. Converts between the state
	// energy the adaptation measures and the drive variance q the Riccati
	// solve needs -- confusing the two inflates q/r by orders of magnitude
	// and puts the loop in the waterbed-amplifying regime (see _note on
	// scaling in fsp.c / doc/devices/fsp.md).
	double Gv[AYLP_FSP_MAX_MODES];
	// filtered modal state estimate Xhat(k|k), laid out per mode as
	// [x_i(k), x_i(k-1)]
	double xhat[AYLP_FSP_MAX_DIM];
	// stationary Kalman gain (column, length dim)
	double L[AYLP_FSP_MAX_DIM];
	// Posterior (filtered) variance of each mode's position state,
	// post_var[i] = P[2i][2i] - L[2i]*PCt[2i] from the converged Riccati
	// solve (not to be confused with the "Pf" filtered-covariance scratch
	// matrix local to fsp_solve_gain_iterate). xhat[2i] alone is a biased
	// estimator of the state's second moment (E[x^2] = Var(xhat) +
	// post_var, not Var(xhat) alone); the adaptation's energy tracker
	// adds this back so q doesn't drift low whenever the posterior
	// variance isn't already negligible.
	double post_var[AYLP_FSP_MAX_MODES];
	// Amortized re-identification solve state. The Riccati gain solve is
	// spread a few iterations per frame (fsp_proc) instead of run in one
	// burst -- the burst stalled frame delivery ~7 ms every adapt_period and
	// showed up as a periodic "source hiccup". adapt_A / adapt_P persist
	// across frames while adapt_solving is set; the previous gain keeps
	// running until the new one converges.
	bool adapt_solving;
	// The normal and event-only Riccati solves share the same workspace and run
	// serially.  adapt_transient_solving identifies the second phase; modal
	// event entry stays inhibited until transient_gain_current is restored.
	bool adapt_transient_solving;
	bool transient_gain_current;
	double adapt_q_scale;
	size_t adapt_it;			// current Riccati iteration
	double adapt_prev_trace;		// trace of P at the last iteration
	double adapt_A[AYLP_FSP_MAX_DIM * AYLP_FSP_MAX_DIM];
	double adapt_P[AYLP_FSP_MAX_DIM * AYLP_FSP_MAX_DIM];
	// command delay line of length `delay + 1`, also retaining the adjacent
	// sample needed by the fractional-delay plant model
	double *ucmd;
	// adaptation accumulators (per mode): EWMA of innovation energy
	// projected onto the mode, and a quadrature demodulator for the local
	// frequency estimate
	double q_ewma[AYLP_FSP_MAX_MODES];
	double demod_re[AYLP_FSP_MAX_MODES];
	double demod_im[AYLP_FSP_MAX_MODES];
	double demod_ph[AYLP_FSP_MAX_MODES];	// running reference phase (rad)
	double r_ewma;				// EWMA innovation variance floor
	// command robustness filter: per-mode pre-compensation at the mode
	// center -- a real-valued extra prediction horizon delay + delay_frac
	// + n_i, n_i = -arg H(f_i)/omega_i cancelling the filter's phase lag,
	// plus a capped gain boost min(1/|H(f_i)|, 3) cancelling its droop.
	// comp_n[i]/comp_frac[i] are that horizon's floor/fractional part
	// (comp_n[i] = target integer step, comp_frac[i] = blend weight for
	// comp_n[i]+1), unified with the delay_frac blend so a fractional n_i
	// is not rounded away when cmd_fc is enabled -- plus the biquad state
	// (direct form II transposed)
	size_t comp_n[AYLP_FSP_MAX_MODES];
	double comp_frac[AYLP_FSP_MAX_MODES];
	double comp_g[AYLP_FSP_MAX_MODES];
	size_t max_steps;	// = max_i (comp_n[i] + (comp_frac[i] > 0))
	double lp_z1, lp_z2;
	// Optional full-band disturbance observer. This is the scalar FIR
	// realization of the delayed Wiener/Kalman predictor: it learns the
	// minimum-variance delay-step estimate directly from the Smith-
	// reconstructed disturbance, including the broadband continuum omitted
	// by the small bank of AR(2) vibration modes.
	double *broad_hist;
	double *broad_w;
	double *broad_w_next;
	// Event-time shadow learner. The production weights above remain frozen;
	// these candidates are evaluated prequentially and are promoted only after
	// sustained held-out improvement with bounded coefficient motion.
	double *shadow_w;
	double *shadow_w_next;
	// Baseline-masked copy of broad_hist: real slots hold the same phi_bl,
	// dark/padded slots hold the boxcar'd DC baseline -- never a
	// prediction, anyone's. This is the ONLY input the dark bank below
	// reads, which is what severs the feedback that made iterating the
	// main FIR on its own output unstable (2026-08-10: 1.45x vs 1.91x at a
	// 25% chop, weight blow-up at broad_mu 0.5).
	double *broad_hist_bl;
	// Direct multi-horizon dark-frame predictor bank (dark_predict), flattened;
	// filter h-1 predicts h frames ahead. On dark frame k, filter k+bd predicts
	// the command-horizon target directly from the window ending at the last
	// real sample. The ordinary observer history remains baseline-masked, so a
	// bank estimate cannot contaminate later live-frame predictions.
	double *dark_w;
	size_t dark_bank_h;		// round-robin training cursor
	double *broad_xbuf;
	// second scratch regressor: the delay+1 predictor's tap window is the
	// delay predictor's window shifted by one sample (they share broad_
	// order - 1 taps), so one ring walk fills both xbuf (delay) and
	// xbuf2 (delay+1) instead of walking the history twice
	double *broad_xbuf2;
	// Provenance of every history sample: 1 = fabricated (padded across a
	// gap or baseline-masked through a dark frame), 0 =
	// reconstructed from a real measurement. A fabricated sample is a
	// legitimate input to a PREDICTION -- it is the best estimate available
	// for that instant -- but it must never appear in an NLMS UPDATE, or
	// the filter trains toward its own output and confirms itself.
	// broad_synbuf/2 carry the same flags for the scratch regressors, so
	// the update can skip exactly those taps.
	unsigned char *broad_syn;
	unsigned char *broad_synbuf;
	unsigned char *broad_synbuf2;
	// observer prefilter ring: the last broad_lp raw phi samples (see
	// broad_lp in aylp_fsp_data)
	double *broad_lpbuf;
	size_t broad_lphead;
	// consecutive real raw samples pushed into broad_lpbuf. The NLMS target
	// phi_bl is the boxcar of the last broad_lp of them, so it is only a
	// clean training target once this reaches broad_lp.
	size_t broad_lpclean;
	// consecutive frames the sensor has flagged, and how many of those were
	// answered with a direct bank command rather than a held command
	size_t broad_dark_run;
	size_t dark_predict_max;	// per-axis cap, derived from ax->delay
	size_t dark_bank_max;		// cap plus command horizon and frac neighbor
	// burst-guard band-pass re-seed request: the first live frame after a
	// discontinuity (stall gap) pre-charges the detector to steady state
	// at that frame's value, so the jump itself cannot ring it
	bool gd_precharge;
	// Explicit slow-disturbance state.  When drift_tau > 0 the broadband
	// and modal predictors see phi - drift_hat, while drift_hat is cancelled
	// separately as a constant-over-the-horizon component.  This keeps a
	// large DC/random-walk component from consuming FIR dynamic range.
	double drift_hat;
	// Optional second slow-state component.  With drift_order=2 this is the
	// estimated disturbance slope (error units/s); drift_hat is propagated by
	// one sample before the innovation update and by the command horizon when
	// forming the cancellation request.  This embeds the ramp internal model
	// in the FSP instead of requiring an outer integral controller.
	double drift_rate;
	// Large-innovation recovery path. The normal command is cross-faded to a
	// bounded proportional/modal command while transient_active is set. The
	// broadband NLMS weights freeze for the event and return cross-fade, leaving
	// an independent stationary model for release detection, then resume on the
	// first fully normal frame.
	double transient_var;
	double transient_threshold;
	bool transient_active;
	bool transient_recovering;
	bool transient_blocked;	// timeout: wait for quiet before re-arming
	double transient_t_event;
	size_t transient_events;
	// Event-only physical-model observer gain.  It is computed from the same
	// damped modes as the normal observer, but with increased process noise so
	// a newly-started ring-down is acquired quickly.
	double transient_L[AYLP_FSP_MAX_DIM];
	enum aylp_fsp_transient_controller transient_controller;
	bool transient_modal;		// selected mode uses the fast modal observer
	double transient_i;		// event-local integral command contribution (V)
	double transient_i_prev;	// rollback point for actuator anti-windup
	double transient_i_step;	// most recent state increment (V)
	double transient_peak_i;
	double transient_start;
	double transient_error_t_last;
	double transient_peak_error;
	double transient_peak_command;
	double transient_error2;
	size_t transient_frames;
	bool shadow_active;
	bool shadow_promoted;
	bool transient_saturated;
	double shadow_primary_e2;
	double shadow_e2;
	double shadow_advantage_start;
	double shadow_error_ratio;
	size_t shadow_frames;
	size_t shadow_promotions;
	// First-order Thiran all-pass state for the fractional command delay:
	// y(k) = a*u(k) + u(k-1) - a*y(k-1).
	double frac_x1, frac_y1;
	// Open-loop operating-point estimate used by the safety detector.  The
	// beam may have a large legitimate static centroid offset, so an absolute
	// error threshold would trip as soon as the loop starts.  During the
	// startup hold this tracks that offset; closed-loop safety compares error
	// magnitude with the learned open-loop magnitude plus trip_error.
	double trip_center;
	size_t trip_count;
	bool trip_warned;
	// Burst guard (see header comment). Band-pass biquad at this axis's
	// regeneration frequency fs/(2*(delay+delay_frac)), DF2T state, fast
	// envelope and slow baseline of the band POWER, and the recovery state
	// machine: guard_gain multiplies this axis's authority (0 during the
	// hold, 0->1 over guard_ramp afterward).
	double gd_f0;			// detector center (Hz), for logging
	double gd_b0, gd_b2, gd_a1, gd_a2;	// band-pass coeffs (b1 = 0)
	double gd_z1, gd_z2;
	double gd_env;			// fast EWMA of bp^2
	double gd_base;			// slow EWMA of bp^2 (quiet frames only)
	// sustained-energy debounce: consecutive over-threshold samples so
	// long as they persist without the requirement below, a broadband
	// low-frequency transient's leakage through the band-pass skirts (a
	// single envelope spike, not a sustained tone) cannot latch a trip
	size_t gd_over_count;
	size_t gd_min_samples;		// required consecutive over-threshold
					// samples before guard_active latches
	bool guard_active;
	bool guard_ramping;		// recovery ramp has begun since the trip
	double guard_t_trip;		// time of the latest trigger (s)
	double guard_t_log;		// last per-event log line (rate limit)
	size_t guard_events;		// activations since start
	size_t guard_frames;		// frames spent at reduced authority
	// stall-gap bookkeeping: slow EWMA of the reconstructed disturbance,
	// used to pad the NLMS history across a frame gap (so the prediction
	// carries the DC correction straight through the hole), and a
	// countdown of frames until fabricated samples have left the tap
	// window (training pauses while it is nonzero).
	double phi_dc;
	size_t broad_fab;
	// diagnostics for the dark-frame path, reported at fini
	size_t n_total;			// frames this axis was visited
	size_t n_dark;			// frames the sensor flagged
	size_t n_dark_predicted;	// of those, filled by the dark bank
	size_t n_dark_cmd;		// of those, commanded through (not held)
	size_t n_train_skip_target;	// frames not trained on: dirty target
	size_t n_syn_taps;		// synthetic taps skipped in updates
};

struct aylp_fsp_data {
	// param type: only "vector" is meaningful
	aylp_type type;
	// units to output (command units, e.g. minmax)
	aylp_units units;
	// loop transport delay in samples (camera + compute + DAC ZOH); the
	// GLOBAL DEFAULT -- each axis object may override with its own
	// "delay"/"delay_frac" (see aylp_fsp_axis)
	size_t delay;
	// Fractional remainder of the measured transport delay. The plant uses a
	// first-order Thiran all-pass (unit magnitude, correct low-frequency group
	// delay); the full-band observer blends adjacent horizon predictions.
	double delay_frac;

	// AUTO DELAY. "delay": "auto" (globally or on one axis) measures the
	// transport delay from the running loop instead of taking it from a Bode
	// sweep, and sets both delay and delay_frac from the measurement.
	//
	// The delay has three parts. Two of them are measured:
	//
	//   COMPUTE. The source device publishes how long the rest of the
	//   pipeline takes per iteration -- centroiding, this controller, the
	//   DAC write, the sinks -- through libaylp/timing.h. It is the only
	//   device that can measure this: it is the one that blocks waiting for
	//   the frame, so it is the only one for which work and waiting are
	//   distinguishable. (From here the two are algebraically inseparable:
	//   the gap since our own last return is exactly the loop period minus
	//   our own duration, whatever anyone else did.) A pipeline whose source
	//   publishes nothing cannot use an auto delay, and init says so.
	//
	//   LOOP PERIOD. Measured here, as the median interval between our own
	//   proc() calls -- the delay wanted is in FRAMES, so the compute time
	//   has to be divided by the period actually achieved, not the one the
	//   config claims. A configured fs more than 10% away from it is
	//   reported: it means every other frame-referenced constant is off too.
	//
	// The third part is not measurable from inside the loop at all: the
	// sensor's own latency (integration midpoint to frame in hand -- the
	// exposure, the readout and the transfer) plus the DAC's zero-order
	// hold, plus whatever mechanical lag the actuator adds. That is
	// delay_auto_bias, in frames. The default 1.5 is one frame of sensor
	// (an exposure comparable to the frame period, pipelined readout) and
	// half a frame of ZOH. CALIBRATE IT ONCE PER BENCH: run a Bode sweep,
	// take the measured delay, subtract the compute term the run logs, and
	// what remains is this number. It is the difference between an estimate
	// and an identification, and the loop is only as good as it.
	//
	// The estimate is installed as soon as its measurement window closes,
	// inside the startup hold and well before the loop closes: the command
	// is parked at zero throughout, and the broadband filter spends the
	// REST of the hold training at the horizon it will actually run at.
	// start_delay is raised at init if it is too short to measure in.
	bool delay_auto;
	// Ceiling the auto estimate is sized and clamped to. Every delay-
	// dependent ring is allocated for this, so it costs memory (and, through
	// the dark bank, a little work) to set it far above the truth. Default
	// 16 samples.
	size_t delay_auto_max;
	// Seconds of the hold to discard before measuring: the first frames of a
	// run include the source's slow first frame, page faults and the
	// observer's own warm-up. Default 0.5.
	double delay_auto_settle;
	// Seconds of measurement after the settle, before the estimate is
	// locked in. Default 1.0.
	double delay_auto_window;
	// Frames of sensor + ZOH + actuator latency added to the measured
	// compute term. Default 1.5; see above -- this is the calibrated part.
	double delay_auto_bias;
	// "delay_auto_bias": "auto" -- derive it from delay_ident_ms instead of
	// carrying a frame count. The bias is a physical latency in SECONDS, so a
	// number of frames is only valid at the frame rate it was written for;
	// this recovers it as (identified delay - measured compute) and converts
	// at the measured rate, which keeps it right when the auto ROI changes
	// the frame rate. Needs delay_ident_ms.
	bool delay_auto_bias_auto;
	// "delay_ident_ms": "auto" -- no identification available, so model the
	// sensor path from what the source publishes (exposure/2 + one frame of
	// readout + half a frame of ZOH) and note the run is on a modelled
	// floor. The calibration suite overwrites it with a measured number the
	// first time the PRBS delay run happens, and the config self-corrects.
	bool delay_ident_auto;
	// The last IDENTIFIED end-to-end transport delay, in milliseconds -- from
	// a Bode sweep or a PRBS run, i.e. the real measurement, in the one unit
	// that does not move with the frame rate. 0 = unset.
	double delay_ident_ms;
	// "fs": "auto" -- take the loop rate from the same measurement the auto
	// delay uses instead of the configured value, and rebuild everything that
	// depends on it (mode coefficients and their Riccati gains, the burst
	// guard, the command filter, the drift/transient EWMA weights). With a
	// self-sizing camera ROI the frame rate is not knowable when the config is
	// written. The configured fs is still used until the measurement lands.
	bool fs_auto;
	// true once the estimate has been applied (or given up on)
	bool delay_auto_done;
	// Measurement state; both arrays are freed at lock-in. auto_dt holds one
	// loop-period sample per frame, auto_compute one per publication from
	// the source (each is already a mean over its own window, so only new
	// values are kept -- storing the same window once per frame would just
	// weight the medians by window length).
	double *auto_dt;
	double *auto_compute;
	size_t auto_n_dt, auto_cap_dt;
	size_t auto_n_c, auto_cap_c;
	double auto_last_c;	// last compute value recorded, to spot new ones
	// sample rate (Hz) used to turn f/zeta into AR coefficients; should
	// match the loop rate
	double fs;
	// output magnitude limit (command units); the symmetric shorthand that
	// sets clamp_lo/clamp_hi to -clamp/+clamp
	double clamp;
	// Actual output bounds, u is held to [clamp_lo, clamp_hi]. Separated
	// from `clamp` so an ASYMMETRIC limit is possible: when the actuator's
	// reachable range is not centred on the command origin (e.g. a mirror
	// spanning -10..+10 V driven from a +2.5 V bias, which is -12.5..+7.5
	// in command units), a symmetric limit either gives up part of the
	// range or lets the loop ask for a voltage the hardware cannot produce.
	// The latter matters here specifically because the clamped command is
	// ALSO what enters the delay ring below -- if something downstream cut
	// it further, this device's plant model would diverge from what the
	// actuator actually received, which is the command-echo mismatch that
	// drives regeneration. Keeping the limit here keeps the two identical.
	// NAN before config resolution means "not explicitly set".
	double clamp_lo, clamp_hi;
	// seconds to hold the command at 0 (loop open) at startup while the
	// Kalman filter converges on the clean open-loop disturbance
	double start_delay;
	// Optional DC-only startup centering. After precenter_delay, the slow
	// drift estimate is cancelled while modal/broadband authority remains off.
	double precenter_delay, precenter_ramp, precenter_clamp;
	// seconds over which the command blends 0 -> full authority after the
	// hold, so the handover is bumpless
	double ramp;
	// adaptation cadence (s); <= 0 disables adaptation (fixed FSP)
	double adapt_period;
	// max per-update center-frequency correction (Hz); 0 freezes f
	double adapt_df_max;
	// EWMA time constant for the adaptation statistics (s)
	double adapt_tau;
	// EWMA weight for the per-mode frequency-demodulator phasor. Derived
	// at init from adapt_df_max (not adapt_tau): the demodulator behaves
	// like a first-order PLL whose capture range is ~1/(2*pi*time
	// constant), so sizing it off the (typically much longer) adapt_tau
	// left it able to reliably track only a small fraction of the
	// adapt_df_max/update it was nominally allowed -- a big, sudden
	// frequency shift (e.g. right after a bias move) would read out a
	// near-arbitrary correction instead of chasing it in. A dedicated,
	// faster time constant widens the capture range to match the cap.
	double demod_beta;
	double demod_tau;	// the time constant demod_beta was derived from;
				// the angle-to-df conversion (fsp_adapt) must
				// use the SAME tau the phasor was accumulated
				// with, not adapt_tau
	// command robustness low-pass cutoff (Hz); 0 disables. See the header
	// comment: rolls off command authority where prediction is worthless,
	// with per-mode phase/gain pre-compensation so the lines still cancel.
	double cmd_fc;
	// Full-band predictor order and NLMS identification step. order=0 keeps
	// the original modal-only controller. When enabled it supplies the
	// command prediction; the modal observer continues to identify/log but
	// is not added again (which would double-count the same disturbance).
	size_t broad_order;
	double broad_mu;
	// NLMS step SCHEDULE, the fix for eigenvalue-spread-limited convergence.
	// The nominal NLMS time constant is broad_order/(broad_mu*fs) = 4.5 s at
	// 512/0.03/3788, but the measured constant on this disturbance is 125-255
	// s -- a 30-60x penalty, because LMS converges each eigenmode in
	// proportion to its own eigenvalue and this spectrum spans ~4.4 decades.
	// A single mu cannot fix that: large mu converges fast but leaves a
	// misadjustment floor of roughly mu/(2-mu), small mu is quiet but slow.
	// So schedule it (Robbins-Monro: large step early, small step late):
	//   mu(t) = broad_mu + (broad_mu_init - broad_mu)*exp(-(t-t_close)/tau)
	// starting from broad_mu_init at loop close and relaxing to broad_mu with
	// time constant broad_mu_tau. The CONVERGED solution is unchanged -- mu
	// affects only the path taken and the steady-state misadjustment, and by
	// the time the scored window opens mu is back at its nominal value. Set
	// broad_mu_init <= 0 (the default) to disable and use a constant
	// broad_mu.
	double broad_mu_init;
	double broad_mu_tau;
	double broad_mu_cur;	// current value, recomputed once per frame
	// Optional offline Wiener initialization.  The text file has one row per
	// tap: index y_w y_w_next x_w x_w_next.  It must contain broad_order rows.
	char *wiener_file;
	// Optional dump of the LEARNED weights at fini, in exactly the format
	// wiener_file reads, so a run's converged predictor can be analysed
	// offline and then replayed as the next run's initialization.  Written
	// on normal AYLP_DONE exit and on SIGINT alike (libaylp/anyloop.c
	// cleanup() runs fini on both paths), so an aborted run still yields
	// whatever the observer had learned.  NULL disables.
	char *wiener_out;
	// Optional periodic convergence trace: one text line every
	// wiener_trace_period seconds with the per-axis tap norms, the norm of
	// the CHANGE since the previous sample, drift, and the cumulative
	// guard/transient counts.  A single end-of-run dump says where the
	// observer ended up; this says whether it had stopped moving, which is
	// the question a settle_time is chosen to answer.  Tiny (~150 B per
	// sample) so the periodic write cannot disturb the RT loop the way a
	// full 512-tap snapshot could.  NULL disables.
	char *wiener_trace;
	double wiener_trace_period;	// s; <= 0 defaults to 10
	FILE *wiener_trace_fp;
	double t_wtrace;		// CLOCK_MONOTONIC of last sample (s)
	// previous-sample copies, for the ||dw|| convergence metric
	double *wtrace_prev[2];
	// Observer band-limit (broad_lp > 0, odd): boxcar prefilter length on the
	// reconstructed disturbance feeding the full-band observer. The NLMS is
	// otherwise broadband to Nyquist, and any K/delay model error leaks the
	// command back into its input as a predictable high-frequency echo that
	// it LEARNS and chases into a self-sustained ring (2026-07-22: ~380 Hz
	// x-axis limit cycle, 4.8 px rms, ignited ~1 min after close with zero
	// external trigger, even with K/delay set from same-day bodes -- the
	// mismatch is amplitude-dependent, so no static tune removes it). The
	// boxcar zeroes the observer's loop gain at fs/broad_lp and attenuates
	// the whole regeneration band; being linear-phase with exactly integer
	// group delay (broad_lp-1)/2, that delay is simply ADDED to the broad
	// prediction horizon, so in-band cancellation timing is unchanged
	// (passband droop at 30 Hz is <1%). 0 disables (raw broadband observer).
	size_t broad_lp;
	size_t broad_gd;	// = (broad_lp-1)/2, derived at init
	// --- frames the sensor could not measure (AYLP_FRAME_REJECTED) ---
	// A centroid stage that cannot fit a frame re-publishes its previous
	// coordinate, so there is no measurement for that instant. Time still
	// advanced, though, and the modal states have always been propagated
	// through it. These three govern what the BROADBAND observer does.
	//
	// dark_predict: predict through dark frames with a DIRECT multi-horizon
	// bank instead of holding. Default true; set false in the config for
	// the masking-only behavior (fill with the DC baseline, hold the
	// command). Two things happen on the k-th consecutive dark frame:
	//
	//  1. The ordinary NLMS history remains filled with the neutral DC
	//     baseline. The first implementation put predictions there and
	//     iterated the MAIN FIR on its own output; that failed twice over
	//     (measured 2026-08-10): iterating an FIR on its own output is an
	//     IIR with no guaranteed stable poles, and the training error
	//     pe = phi_bl - pred1 carried the fabricated taps' contribution,
	//     so the real weights grew to compensate and made the next hole's
	//     synthesis wilder -- 1.45x vs masking's 1.91x at a 25% chop, and
	//     a 5x AMPLIFYING blow-up at broad_mu 0.5. The bank severs both
	//     loops by construction: a direct h-step prediction reads only
	//     real samples and neutral baseline fills, never any prediction.
	//     Each bank filter trains on that same masked signal against real,
	//     clean targets only, and its output never enters observer history.
	//
	//  2. The loop keeps COMMANDING from bank horizon k+delay+broad_gd,
	//     blended with the adjacent horizon for delay_frac. This single
	//     direct prediction replaces the unstable cascade of a k-step fill
	//     followed by the main FIR. Holding instead
	//     (the old behavior, and the fallback whenever the bank cannot
	//     predict, the guard is active, or a transient is running) leaves
	//     the loop open for the whole hole -- 26% of the time at a 25%
	//     duty chop.
	bool dark_predict;
	// Frames of prediction before falling back to the DC baseline and the
	// held command. The bank also stores delay+broad_gd+1 longer horizons
	// for direct command prediction. 0 = auto, per axis:
	// 2*(delay + broad_gd), i.e. twice the horizon the main filter runs.
	size_t dark_predict_max;
	// Bank horizons trained per live frame (round-robin), bounding the
	// added per-frame cost at dark_bank_train * broad_order MACs. Default
	// 4. All stored horizons are refreshed round-robin.
	size_t dark_bank_train;
	// Fraction of the NLMS regressor that must be real measurements before
	// the filter is allowed to train on it at all. A partial-update NLMS
	// step normalizes by the norm of the taps it moves, so a regressor with
	// only a few real taps left produces an enormous step and one small
	// sample takes the whole correction. Counted in TAPS, not energy:
	// fabricated samples sit near the DC baseline and carry almost no
	// energy, so an energy test would call a window of them clean. Default
	// 0.5.
	double dark_train_min_real;
	// Consecutive dark frames after which training pauses entirely until the
	// fabricated samples have left the tap window (the broad_fab countdown).
	// This is for an OUTAGE; ordinary chops are handled tap by tap by
	// dark_train_min_real. 0 = auto, broad_hist_len/2. Do not set this to a
	// small number of frames: broad_fab counts down one frame per frame, so
	// a threshold a repeating chop keeps re-arming never reaches zero and
	// the filter never trains again.
	size_t dark_fab_frames;
	// Separate slow drift from the vibration predictor. <= 0 preserves the
	// original compound-disturbance predictor.
	double drift_tau;
	double drift_beta;	// derived from drift_tau and fs
	unsigned drift_order;	// 1 = EWMA position, 2 = alpha-beta position/rate
	// Innovation-triggered transient/reacquisition path. transient_sigma <= 0
	// disables it.  By default a bounded P servo replaces prediction.  With
	// transient_modal_q_scale > 0, a high-process-noise copy of the damped-mode
	// observer estimates and predicts ring-down motion through the plant delay;
	// transient_kp remains as a small correction for model error.
	double transient_sigma;
	double transient_floor;
	// Output-error band used only for event scoring (not event detection).
	double transient_settle_error;
	double transient_kp;
	double transient_ki;		// integral gain (1/s)
	double transient_i_leak;	// integrator leakage rate (1/s); 0 = pure
	double transient_i_limit;	// integral contribution limit (command units)
	// `manual` honors the parameter values, `off` inhibits all event recovery,
	// and `auto` enables the conservative default detector when sigma is unset.
	enum aylp_fsp_transient_mode transient_mode;
	enum aylp_fsp_transient_controller transient_controller;
	bool transient_controller_set;
	double transient_modal_q_scale;
	// Alternate P/modal recovery on successive per-axis events for a same-run
	// push A/B. Requires transient_modal_q_scale > 0.
	bool transient_modal_ab;
	double transient_tau;
	double transient_beta;	// quiet innovation-variance EWMA weight
	double transient_hold;
	double transient_ramp;
	// Optional wider command envelope used only while transient_active. NAN
	// means inherit the normal clamp. The separate trip threshold and hard
	// timeout keep this authority out of the continuously running predictor.
	double transient_clamp_lo, transient_clamp_hi;
	double transient_trip_command;
	double transient_max_duration;
	// Inhibit event detection until this many seconds after the startup hold.
	// Baseline variance and the closed-loop predictor continue learning during
	// this warm-up; only the trigger is gated.
	double transient_arm_delay;
	// Optional event-time shadow NLMS. The main FIR stays frozen. A shadow copy
	// learns conservatively and can replace it only after sustained pre-update
	// prediction improvement, bounded coefficient drift, no guard suppression,
	// and no event saturation. <=0 mu disables shadow learning/promotion.
	double transient_shadow_mu;
	double transient_shadow_tau;
	double transient_shadow_min_duration;
	double transient_shadow_hold;
	double transient_shadow_ratio;
	double transient_shadow_norm_ratio;
	double transient_shadow_beta;
	size_t shadow_promotions;
	// Optional one-row-per-event CSV summary for push/ring-down trials.
	char *transient_log;
	FILE *transient_log_fp;
	// Legacy optional mode: when true, NLMS weights freeze as soon as the
	// startup hold ends. Closed-loop-adaptive configurations set this false.
	bool broad_freeze_closed;
	// Non-latching limit diagnostics. While authority is nonzero, either error
	// magnitude exceeding the learned open-loop magnitude by trip_error, or an
	// excessive requested command, for trip_frames consecutive samples emits a
	// warning. The configured clamps remain authoritative. AYLP_BEAM_LOST holds
	// output at zero until the sensor clears it, then authority ramps back in.
	double trip_error;
	double trip_command;
	size_t trip_frames;
	bool beam_lost;
	double beam_recover_ramp;
	double beam_recover_start;
	// Burst guard tuning (see header comment). guard_ratio is an AMPLITUDE
	// ratio: trigger when the band envelope exceeds guard_ratio times the
	// quiet baseline (or the guard_floor, whichever is larger; floor is in
	// normalized error units). <= 0 disables the guard.
	double guard_ratio;
	double guard_floor;
	double guard_hold;	// s at zero authority after a trigger
	double guard_ramp;	// s to ramp authority back to 1
	double guard_tick;	// s between ticker lines; 0 silences the ticker
	double guard_beta_fast;	// EWMA weights derived at init
	double guard_beta_slow;
	// Sustained-energy requirement (cycles at this axis's gd_f0) the
	// envelope must stay over threshold before a trip latches, so a
	// single-envelope-window transient (measured cause of a past false
	// trip: low-frequency pointing content leaking through the band-pass
	// skirts) can't fire the guard the way a genuine sustained
	// regeneration ring does. <= 0 keeps the old single-sample trigger.
	double guard_min_cycles;
	double t_tick;		// last ticker print (s)
	size_t guard_events;	// total activations, both axes
	// Stall-gap handling (see header comment). Patch histories when the
	// wall-clock gap between consecutive proc calls exceeds gap_trip
	// seconds; <= 0 disables.
	double gap_trip;
	// "gap_trip": "auto" -- 1.5 measured frame periods. A trip level in
	// SECONDS only means one thing at one frame rate: too tight and every
	// normal frame is a "gap", too loose and a real dropped frame is
	// invisible. With a self-sizing ROI the frame rate is not knowable when
	// the config is written, so this follows the same measurement fs does.
	bool gap_trip_auto;
	double gap_dc_beta;	// EWMA weight for the per-axis phi_dc estimate
	double t_last;		// time of the previous proc call (s)
	size_t gap_events;	// stall gaps detected since start
	size_t n_closed;	// frames processed since the loop first closed
	// shared biquad coefficients for the command filter (normalized)
	double lp_b0, lp_b1, lp_b2, lp_a1, lp_a2;

	// per-axis controllers: [0] = y (element 0), [1] = x (element 1)
	struct aylp_fsp_axis axis[2];

	// timing / handover bookkeeping
	double t0;		// CLOCK_MONOTONIC of first proc (s)
	double t_close;		// time the loop first closed (s); 0 until then
	double t_adapt;		// time of last adaptation (s)
	bool closed;		// true once past start_delay (logged once)
	size_t n_seen;		// frames processed since start

	// number of pipeline elements the last proc saw (for the passthrough
	// guard when the vector isn't 2-wide)
	size_t n_elem;
	// result vector
	gsl_vector *res_v;
};

// initialize fsp device
int fsp_init(struct aylp_device *self);

// process fsp device once per loop
int fsp_proc(struct aylp_device *self, struct aylp_state *state);

// close fsp device when loop exits
int fsp_fini(struct aylp_device *self);

#endif
