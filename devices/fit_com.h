#ifndef AYLP_DEVICES_FIT_COM_H_
#define AYLP_DEVICES_FIT_COM_H_

#include <stdbool.h>
#include "anyloop.h"

// Number of model parameters: y0, x0, slope_y, slope_x, sigma, amplitude,
// background. Fixed, so the normal equations live on the stack.
//
// The beam is modelled as CIRCULAR: one sigma, shared by both axes. Splitting
// it into sigma_y/sigma_x means a real elliptical fit -- separate decay
// constants through the rowwalk recurrences and a second sigma column in the
// Jacobian -- not just two enumerators. NP also sizes NGRAM, JtJ and the
// Cholesky, so bumping it alone widens the normal equations by a parameter
// nothing fits.
#define AYLP_FIT_COM_NP 7
enum {
	AYLP_FIT_P_Y0 = 0,
	AYLP_FIT_P_X0,
	AYLP_FIT_P_SY,		// px of y motion per image row of readout
	AYLP_FIT_P_SX,		// px of x motion per image row of readout
	AYLP_FIT_P_SIGMA,
	AYLP_FIT_P_AMP,
	AYLP_FIT_P_BG,
};

// Joint spatio-temporal beam fit -- see doc/devices/fit_com.md.
//
// A rolling-shutter frame does not have "a" beam position: every row is
// captured at a different instant, so a frame in which the beam moves is
// sheared, and any single position extracted from it is ill-defined. Rather
// than measure one position and then try to undo the shear afterwards (which
// is what anyloop:wfs_com's row regression attempts, from 2-3 row-group means
// -- far too few points to condition a slope), this device puts the motion
// INSIDE the model and fits it directly:
//
//   I(row,col) = bg + amp * exp(-[(row-y(row))^2 + (col-x(row))^2] / 2 sigma^2)
//     y(row) = y0 + slope_y * (row - ref_row)
//     x(row) = x0 + slope_x * (row - ref_row)
//
// Seven parameters against every pixel in the window -- overdetermined by
// orders of magnitude, where the row-group regression had two parameters from
// three points. Shear stops being a defect to be detected and rejected, and
// becomes the signal that identifies the slope. Nothing is ever discarded for
// moving, so there is no shear gate, and therefore none of the
// velocity-correlated frame dropout that a shear gate necessarily produces (a
// shear threshold is a velocity threshold, so it preferentially drops frames
// at the fast part of every oscillation and keeps the turning points).
struct aylp_fit_com_data {
	// --- model state ---
	double p[AYLP_FIT_COM_NP];	// current parameter vector
	bool active[AYLP_FIT_COM_NP];	// slopes are inactive if !fit_slope
	bool acquired;

	// --- window (0 => whole image) ---
	size_t window_h, window_w;	// as configured
	size_t win_h, win_w;		// actual, clamped to the image
	size_t win_y, win_x;		// window centre, image coords
	long init_y, init_x;		// <0 => acquire from the brightest pixel

	// --- core box: the part of the window the solver actually iterates on
	// Past a few sigma the model is flat to well below one count, so those
	// pixels say nothing about position, width or amplitude -- only about
	// the background, and that contribution is summable in closed form
	// (tail_* below). Iterating on the core instead of the window is what
	// decouples per-frame cost from sensor size: a 248x248 frame and a
	// 32x32 one do the same solver work. Re-planned every frame from the
	// warm-started sigma and slopes, so it always contains the beam.
	size_t core_y, core_x;		// origin, image coords
	size_t core_h, core_w;		// extent
	size_t core_stride;		// core_w padded up to a vector width
	size_t core_cap;		// allocated doubles per buffer
	size_t core_cap_stride;		// stride the scratch was sized for
	double fit_radius;		// core half-width in sigmas
	size_t max_core;		// hard cap on either core dimension

	// Window minus core, which enters the fit only through bg: with the
	// model flat there, sum((bg-I)^2) = n*bg^2 - 2*bg*sum(I) + sum(I^2),
	// so three numbers gathered once per frame stand in for every one of
	// those pixels at every iteration. Sampled over tail_rows rows and
	// scaled, so the gather costs the same on a 248x248 frame as on a
	// 32x32 one.
	double tail_n, tail_s, tail_s2;
	size_t tail_rows;

	// --- solver tuning ---
	size_t max_iter;
	double tol;			// relative cost improvement to stop at
	double max_us;			// latency guard; 0 = bounded by max_iter
	double lambda;			// carried between frames as a warm start
	double robust_k;		// Tukey biweight cutoff in residual
					// sigmas; 0 = no reweighting
	size_t robust_iter;		// plain iterations before reweighting starts

	// --- model bounds / validity gates ---
	double sigma_init, sigma_min, sigma_max;
	double min_amplitude;		// fitted amp below this => no beam
	double max_residual;		// rms residual (counts) above this => reject
	double max_step;		// maximum accepted centre motion per frame (px), 0 off
	bool moment_output;		// report robust intensity centroid
	bool fit_gaussian;		// run Gaussian LM solver; false requires moment_output
	double moment_cut;		// counts above fitted background used by centroid
	size_t moment_col_stride;	// sample every Nth column in moment pass
	double moment_max_y_skew;	// identify rolling-shutter fragments; 0 disables
	double moment_min_y_width;	// fragment if measured/fitted row width is below this
	bool moment_reject_inferred;	// publish inferred centre but exclude it from control
	double pwm_period_frames;	// rolling-shutter PWM period; 0 disables phase tracker
	double pwm_dark_flux, pwm_bright_flux;
	size_t pwm_full_start, pwm_full_end_margin, pwm_phase;
	bool pwm_prev_dark, pwm_phase_locked, pwm_filter_have;
	double pwm_y, pwm_x, pwm_py, pwm_px;
	double pwm_qy, pwm_qx, pwm_ry_full, pwm_rx_full;
	double moment_dy, moment_dx;	// filtered full-frame motion for sliced-frame prediction
	double moment_full_y, moment_full_x;
	size_t moment_full_age;
	bool moment_have_full;
	double *moment_row_template;
	size_t moment_template_cap, moment_template_n, moment_template_org_y;
	double moment_template_y;
	size_t reacquire_after;
	size_t lost;

	// --- timing ---
	// Seconds per image row of readout. The POSITION fit does not need it --
	// slope is fitted in px/row and the reported position is the model
	// evaluated at ref_row, both of which are independent of it. It converts
	// the slope to a physical velocity and fixes the epoch the reported
	// position belongs to, which is what a delay-modelling controller
	// downstream needs. 0 => velocity is reported as 0 and nothing else
	// changes.
	double row_time;
	double ref_row;			// fixed epoch, image coords
	double vel_y, vel_x;		// px/s, last good fit (diagnostic)

	// --- last good output, normalized [-1,1] ---
	double last_y, last_x;
	double last_rms;		// residual rms of the last accepted fit
	bool inferred_last;		// accepted output inferred a shutter-sliced centre
	size_t n_iter_last;
	bool budget_hit;		// max_us cut the last frame short

	// --- scratch, core-sized, 32-byte aligned, row stride core_stride ---
	// The pad columns beyond core_w carry pixel 0 and weight 0 so the
	// kernels can run whole vectors off the end of a row without a tail
	// loop and without contributing anything.
	double *pix;			// core pixels as doubles, gathered once
	// The ROOT of each robust weight, not the weight. A weighted
	// least-squares problem is the unweighted one on rows scaled by
	// sqrt(w), which is what lets eval_jac build the normal equations as a
	// plain Gram matrix; and the Tukey weight is a square by construction,
	// so its root costs nothing to keep instead.
	double *weights;
	double *resid;			// per-pixel residual, core_h*core_stride
	// spare residual buffer: a trial step writes here and, if accepted,
	// swaps rather than recomputing what it already walked the core for
	double *resid_alt;
	// whitened design matrix for one block of core rows, NCOL*JAC_BLOCK
	// rows' worth; see eval_jac
	double *scratch;

	gsl_vector *com;		// output [y, x]
};

int fit_com_init(struct aylp_device *self);
int fit_com_proc(struct aylp_device *self, struct aylp_state *state);
int fit_com_fini(struct aylp_device *self);

#endif
