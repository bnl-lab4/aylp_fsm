#ifndef AYLP_DEVICES_KALMAN_FILTER_H_
#define AYLP_DEVICES_KALMAN_FILTER_H_

#include <stdbool.h>
#include <time.h>

#include "anyloop.h"

// Adaptive linear predictor: replaces each element of the state vector with a
// prediction of itself `horizon` samples in the future, so a downstream
// controller acts on where the error is GOING to be instead of where it was
// one loop-delay ago. Set horizon to the total command->error latency in
// samples and the predictable part of the disturbance (vibration lines,
// resonant humps -- anything with autocorrelation longer than the delay) is
// corrected as if the delay weren't there; the unpredictable white part is
// passed through untouched (the predictor's least-squares solution learns to
// ignore it). Measured on the FSM loop 2026-07-10: the closed-loop error is
// ~50% predictable at the 2 ms loop delay, so the ceiling is ~2x.
//
// Each element gets an order-p AR model of the signal, predicting s(t+k) from
// [s(t) ... s(t-p+1)], with the weights adapted online by NLMS. (True
// RLS/Kalman weight estimation was tried first and is the natural fit on
// paper, but with exponential forgetting on a narrowband signal its
// covariance grows without bound in the unexcited directions -- classic
// covariance wind-up; it diverged within seconds on recorded loop data.
// NLMS is the robust stochastic approximation of the same estimator, costs
// O(order) per element, and reached the ideal block-least-squares prediction
// bound on the same data with mu ~= 0.02-0.05.)
//
// The handover is ramped: output = a*prediction + (1-a)*input with `a`
// rising linearly from 0 to `gain` over `ramp` seconds once training starts.
// Weights start at zero, so an instant handover would feed the controller a
// near-zero error signal (loop effectively open) until the model converges;
// the ramp keeps the loop closed on the raw error while the predictor earns
// its authority.
//
// Put it after the sensor (center_of_mass) and any logging sink, before
// anyloop:pid. It trains continuously on its own input, so it keeps tracking
// the closed-loop statistics as the loop and environment change.
struct aylp_kalman_filter_data {
	// param type: only "vector" makes sense here
	aylp_type type;
	// units to output (default unchanged)
	aylp_units units;
	// p: AR order (taps of history per element)
	size_t order;
	// k: samples ahead to predict; set to the loop's total latency
	size_t horizon;
	// NLMS step size, stable for 0 < mu < 2; smaller = lower steady-state
	// misadjustment but slower convergence/tracking
	double mu;
	// final blend: output = gain*prediction + (1-gain)*input. 1 = full
	// prediction; use <1 for margin while validating. Base gain is the x
	// axis (elements >= 1); gain_y is element 0 (y), inheriting the base
	// gain if not given, same convention as pid. 0 disables prediction on
	// that axis (training continues either way). Measured 2026-07-10: the
	// predictor helps x (lower crossover, line-dominated) but pumps y's
	// 30.8 Hz line 3x -- y's high-gain loop (~62 Hz crossover) reacts
	// badly to the extra phase dynamics, so run x-only with gainy 0.
	double gain, gain_y;
	// seconds over which the blend ramps 0 -> gain after training starts
	double ramp;
	// output magnitude limit (predictions can overshoot on transients)
	double clamp;
	// seconds to pass input through untouched (no training) at startup;
	// match the pid/DAC start_delay so the open-loop settling transient
	// doesn't poison the weights
	double start_delay;
	// CLOCK_MONOTONIC time of the first proc (s); 0 until seeded
	double t0;
	// CLOCK_MONOTONIC time of the first training update; 0 until seeded
	double t_train;
	// number of elements (allocated size of the per-element arrays)
	size_t n;
	// per-element ring buffer of the last (order + horizon) inputs,
	// laid out [n][order+horizon]
	double *hist;
	size_t hist_len;	// = order + horizon
	size_t head;		// ring write index (shared by all elements)
	size_t n_seen;		// samples since training started
	bool engaged;		// true once the ramp has started (logged once)
	// per-element NLMS weights, [n][order]
	double *w;
	// scratch window copy, [order]
	double *xbuf;
	// place to put the result
	gsl_vector *res_v;
	// transient bypass: while |input| > bypass, and for hist_len samples
	// after it last was (until the regression window has flushed the
	// transient), pass the RAW input through and freeze training. A large
	// disturbance (bumped mount, acquisition) is out-of-distribution for
	// the AR model: the prediction collapses toward zero while the beam is
	// still displaced -- which fools any downstream large-error logic (the
	// pid's fast mode) into disengaging early -- and training on it
	// scrambles the weights, parking the loop off-null until they
	// re-converge. <= 0 disables (default). Units = input error units.
	double bypass;
	// per-element countdown of remaining bypassed samples
	size_t *bypass_cnt;
};

// initialize kalman_filter device
int kalman_filter_init(struct aylp_device *self);

// process kalman_filter device once per loop
int kalman_filter_proc(struct aylp_device *self, struct aylp_state *state);

// close kalman_filter device when loop exits
int kalman_filter_fini(struct aylp_device *self);

#endif
