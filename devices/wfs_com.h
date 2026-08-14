#ifndef AYLP_DEVICES_WFS_COM_H_
#define AYLP_DEVICES_WFS_COM_H_

#include <stdbool.h>
#include "anyloop.h"

// Exhaustive integer-pixel correlation search radius is capped at compile
// time so the per-subaperture NCC grid can live on the stack instead of
// being heap-allocated every frame.
#define AYLP_WFS_COM_MAX_SEARCH_RADIUS 6
#define AYLP_WFS_COM_MAX_GRID (2 * AYLP_WFS_COM_MAX_SEARCH_RADIUS + 1)

// How far, in pixels, a subaperture's own shift estimate may sit from the
// frame's consensus before it is treated as distorted rather than displaced:
// excluded from the shift that re-registers the reference templates, and left
// frozen rather than re-registered itself. This is a consistency tolerance on
// a RIGID-TRANSLATION model, not a motion limit -- a real translation of any
// size moves every subaperture alike and so produces no disagreement at all.
// Half a pixel is the scale below which the sub-pixel fit is trustworthy;
// looser values readmit exactly the ill-conditioned estimates the tolerance
// exists to keep out of the templates.
#define AYLP_WFS_COM_CONSENSUS_TOL 0.5

// Subaperture (Shack-Hartmann-style) matched-filter tip/tilt tracker.
// See doc/devices/wfs_com.md for the algorithm and its scope/limitations.
struct aylp_wfs_com_data {
	// --- subaperture grid geometry (all fixed at init from params) ---
	size_t subap_h, subap_w;	// pixel size of one subaperture
	size_t subap_rows, subap_cols;	// grid shape
	size_t n_subaps;		// = subap_rows * subap_cols
	size_t search_radius;		// max integer-pixel correlation offset
	bool reject_edge_matches;	// discard NCC peaks on search boundary
	size_t ext_h, ext_w;		// = subap_h/w + 2*search_radius
	size_t win_h, win_w;		// total window incl. search margin;
					// = subap_rows*subap_h + 2*search_radius
					// (and the w equivalent)
	unsigned char threshold;

	// --- matched-filter / confidence tuning ---
	double ref_beta;		// EWMA rate for reference template update
	// NLMS-style step SCHEDULE (mirrors anyloop:fsp's broad_mu_init /
	// broad_mu_tau): a fixed ref_beta has to trade fast initial
	// convergence against steady-state quietness, same as any EWMA. Set
	// ref_beta_init > ref_beta to start fast right after (re)acquisition
	// and relax to the steady-state ref_beta with time constant
	// ref_beta_tau; ref_beta_init <= 0 (the default) disables the
	// schedule and uses a constant ref_beta throughout, as before.
	double ref_beta_init;
	double ref_beta_tau;
	double t_acquire;		// CLOCK_MONOTONIC time of last (re)acquire
	double min_confidence;		// per-subaperture normalized-correlation
					// floor to be included in reconstruction
	double flux_floor;		// per-subaperture minimum flux (sum of
					// thresholded counts) to attempt a match
	size_t min_valid_subaps;	// frame is "no signal" below this many
					// confident subapertures

	// --- rolling-shutter correction ---
	// Explicit enable for the row-index reconstruction. Prefer the
	// rolling_shutter JSON boolean; row_time remains a backwards-compatible
	// alias whose positive magnitude has no calibration meaning.
	bool rolling_shutter;
	// A positive legacy value enables it: subapertures are first averaged
	// PER ROW GROUP, then a weighted linear fit of (row-group shift) vs
	// row-group row index is extrapolated to the window's vertical-centre
	// row instead of blending
	// every row's measurement -- taken at genuinely different instants
	// under a rolling shutter -- into one number.
	//
	// The MAGNITUDE of row_time does not matter, only its sign: both the
	// the old regression's independent variable and evaluation point scaled
	// by the same factor, which cancelled exactly. This was verified in
	// contrib/diagnostics/test_wfs_com_rolling_shutter.c -- every positive value
	// tested gave identical tracking error). So this is currently a
	// row-index-space correction, not a genuine physical-time one, and
	// row_time is retained only as a legacy enable switch. See the docs.
	double row_time;
	// per-row-group scratch accumulators, sized subap_rows at init
	double *row_sum_w;
	double *row_sum_wt;		// weight * actual matched row coordinate
	double *row_sum_wdy;
	double *row_sum_wdx;

	// --- rolling-shutter REJECTION (independent of correction above) ---
	// <= 0 (default) disables this. > 0: every frame, the per-row-group
	// weighted-mean shift (the same row_sum_w/wdy/wdx the row regression
	// uses) is computed regardless of correction enablement, and if the
	// spread between row-groups' shift estimates exceeds max_row_shear
	// pixels in y or x, the WHOLE frame is held (like a lost-signal frame)
	// instead of reporting a shear-distorted position. Unlike the
	// lost-signal path, this does NOT count toward reacquire_after: the
	// geometry hasn't changed, only this frame's readout timing, so a
	// sustained rolling-shutter condition (e.g. continuous vibration) must
	// not force a repeated window/template reset.
	double max_row_shear;
	size_t shear_rejected;		// running count, for diagnostics/logging

	// --- window/acquisition state (mirrors center_of_mass track mode) ---
	long init_y, init_x;		// <0 => acquire from brightest pixel;
					// configured point is used on every reacquire
	bool acquired;
	bool just_acquired;		// true for the one warm-up frame right
					// after (re)acquisition, when every
					// subaperture is still seeding its
					// reference template and so has no
					// shift estimate yet -- must not count
					// as a lost-signal frame
	size_t win_y, win_x;		// integer window centre, image coords
	double abs_y, abs_x;		// running sub-pixel absolute position,
					// image coords
	double ref_off_y, ref_off_x;	// acquisition centroid minus integer
					// window centre; correlation residuals are
					// measured relative to this fixed anchor
	double last_y, last_x;		// last valid normalized [-1,1] output
	size_t lost;			// consecutive no-signal frames
	size_t reacquire_after;

	// --- per-subaperture matched-filter state ---
	// self-calibrating reference templates, one subap_h*subap_w block per
	// subaperture, flattened row-major and concatenated by subaperture
	// index (r*subap_cols+c)
	double *ref;
	// Snapshot taken before processing an established frame. Reference
	// updates are speculative until the frame-level validity/shear gates
	// pass; rejected or incomplete frames restore this copy.
	double *ref_backup;
	bool *ref_set;		// whether each subaperture's template is seeded
	bool *ref_set_backup;
	// Which subapertures produced a usable match this frame, i.e. cleared
	// every per-subaperture gate, and what each one measured. Their
	// templates are re-registered in the second pass once the frame's
	// consensus shift is known; everyone else's is left frozen.
	bool *matched;
	double *sub_dy, *sub_dx;

	// scratch: the extended (subap + 2*search_radius margin) thresholded
	// block for whichever subaperture is currently being processed
	double *ext;

	gsl_vector *com;	// output [y, x]
};

int wfs_com_init(struct aylp_device *self);
int wfs_com_proc(struct aylp_device *self, struct aylp_state *state);
int wfs_com_fini(struct aylp_device *self);

#endif
