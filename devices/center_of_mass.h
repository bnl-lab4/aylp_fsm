#ifndef AYLP_DEVICES_CENTER_OF_MASS_H_
#define AYLP_DEVICES_CENTER_OF_MASS_H_

#include "anyloop.h"
#include "thread_pool.h"

// src bundle passed to com_mat_uchar via the task runner
struct aylp_com_src {
	gsl_matrix_uchar mat;
	unsigned char threshold;
};

struct aylp_center_of_mass_data {
	// param: height of regions/subapertures
	size_t region_height;
	bool auto_region_height;
	// param: width of regions/subapertures
	size_t region_width;
	bool auto_region_width;
	// param; set to 1 for no multithreading
	size_t thread_count;
	// param: subtract this value from each pixel before computing CoM
	unsigned char threshold;
	bool auto_threshold;
	// param (track mode): the brightest pixel in the window must reach this
	// for the frame to count as containing the beam. 0 disables the test,
	// which is the pre-2026-08-02 behaviour: any sum above zero was treated
	// as a beam, so read noise a few counts above `threshold` produced a
	// centroid of pure noise. Distinct from `threshold`, which only shapes
	// the weighting -- this decides whether the frame is used at all.
	unsigned char min_peak;
	bool auto_min_peak;
	double auto_mean_peak;
	size_t auto_peak_samples;
	// param (track mode): reject a frame when its dimmest significant row
	// falls below this fraction of the frame's own typical row, both
	// measured against a learned reference profile. 0 disables the test.
	//
	// This is the gate for a partially exposed beam. A rolling shutter
	// reading out across a chopped source produces frames in which only some
	// rows caught the source's on-window, so the beam comes back cut by a
	// hard horizontal edge. Neither brightness test sees that: the rows that
	// were lit are at full brightness, so `min_peak` passes, and a fragment
	// still carries plenty of flux.
	//
	// The test asks the question directly -- are some rows lit while others
	// are not? -- by keeping a reference row profile of what a whole beam
	// looks like and dividing this frame's profile by it. On a whole beam
	// every row comes out at the same ratio, whatever that ratio happens to
	// be, so the frame's overall brightness divides out exactly rather than
	// approximately: a dimmer frame is uniformly dimmer. On a cut beam the
	// ratios split in two, near the frame's level on the rows that were lit
	// and near zero on the rows that were not, and the gap between them is
	// what this parameter thresholds. Because the reference is the beam's
	// own measured profile, no assumption about beam shape or size enters,
	// and a tight spot is handled the same as a broad one.
	//
	// Assumes the shutter rolls along rows, which is how the ASI sensors
	// read out. The centroid of a cut frame is pulled toward the surviving
	// rows -- a systematic error along that same axis, not noise that
	// averages out.
	double ref_cut;
	// param: frames used to learn the reference before the gate switches on.
	// The bootstrap takes the row-wise MAXIMUM over these frames rather than
	// the mean, because the beam may already be chopping while it runs: a
	// cut only ever removes light, so the largest value each row reaches
	// over enough frames is that row's uncut value. A mean would learn a
	// blend of whole and cut beams and set the reference too low.
	size_t ref_warmup;
	// param: EMA rate at which accepted frames update the reference
	// afterwards, so it follows slow drift in power, focus and alignment.
	// Only accepted frames contribute, which is what keeps cut frames from
	// pulling the reference down toward themselves.
	double ref_rate;
	// param: rows whose reference is below this fraction of the brightest
	// reference row are ignored. Their ratio is a small number over a small
	// number, i.e. noise, and including them would swamp the test. It also
	// keeps the test off the profile's outer skirts, where beam motion moves
	// a row's flux the most in relative terms.
	double ref_floor;
	// learned reference row profile, region_height entries, in window
	// coordinates -- the window follows the beam, so the beam sits at a
	// stable place in it and the profile is stationary there
	double *ref;
	// this frame's row profile, and scratch for the per-row ratios
	double *rows;
	double *rho;
	// frames folded into the bootstrap so far, and whether it has finished
	size_t ref_seen;
	bool ref_ready;

	// param: confine the sum to a single region that follows the previous
	// center of mass, ignoring everything outside it
	bool track;
	// Optional pattern-agnostic translation output. A rolling keyframe supplies
	// local spatial structure; per-frame gain and offset are nuisance parameters,
	// so brightness and background changes cannot masquerade as motion.
	bool registration;
	double *registration_ref;
	size_t *registration_sample_y, *registration_sample_x;
	double *registration_sample_ref;
	double *registration_sample_gy, *registration_sample_gx;
	size_t registration_n_samples;
	size_t registration_sample_capacity;
	double *registration_residuals;
	double *registration_errors;
	double registration_info_need;
	double registration_condition;
	bool registration_ready;
	double registration_quality;
	double registration_noise, registration_residual_scale;
	size_t registration_rolls;
	double registration_anchor_y[31], registration_anchor_x[31];
	size_t registration_anchor_n, registration_anchor_pos;
	double registration_y, registration_x;
	double registration_ref_y, registration_ref_x;
	// param: initial window centre in image pixels; <0 means acquire from
	// the brightest pixel on the first frame
	long init_y;
	long init_x;
	// param: consecutive frames of zero signal before we re-acquire
	size_t reacquire_after;
	bool auto_reacquire;
	// param: seconds to run with the wide acquisition window before narrowing
	// to region_height/region_width (0 disables the acquisition phase)
	double acquire_seconds;
	// param: acquisition window size; 0 means the whole image
	size_t acquire_height;
	size_t acquire_width;
	// true while the wide acquisition window is in effect
	bool acquiring;
	// CLOCK_MONOTONIC time the current acquisition phase started (s)
	double acquire_t0;
	// current window centre, in image pixel coordinates
	size_t win_y;
	size_t win_x;
	// whether the window has been placed yet
	bool acquired;
	// last valid output, held while the beam is lost
	double last_y;
	double last_x;
	// consecutive frames with no signal inside the window
	size_t lost;
	// whether the previous frame was judged to hold the beam, so the
	// lost/reacquired transitions are logged once instead of every frame
	bool had_beam;
	// frames seen and frames whose centroid was held (track mode), so fini
	// can report the reject rate. A gate meant to drop chopped frames trips
	// by design, and the fraction it drops is the number that decides
	// whether gating is viable at all: every held frame feeds the loop a
	// stale error that the integrator keeps acting on, which reads
	// downstream as added delay.
	size_t n_frames;
	size_t n_held;
	// number of distinct loss episodes, including the ones whose log line
	// was suppressed by the rate limiter
	size_t n_episodes;
	// CLOCK_MONOTONIC time of the last beam-lost/beam-back log line (s).
	// Those transitions are once-per-episode, which is quiet for a beam that
	// occasionally drops out and a flood for a chopped one -- at a few
	// hundred episodes a second the logging alone would cost more than the
	// centroid. Rate-limited to COM_LOG_INTERVAL.
	double last_log_t;
	// array of threads
	pthread_t *threads;
	// array of tasks
	struct aylp_task *tasks;
	// number of tasks
	size_t n_tasks;
	// array of com_src inputs for the tasks (matrix view + threshold)
	struct aylp_com_src *com_srcs;
	// queue of tasks
	struct aylp_queue queue;

	// center of mass result (contiguous vector)
	gsl_vector *com;
};

// initialize center_of_mass device
int center_of_mass_init(struct aylp_device *self);

// process center_of_mass device once per loop
int center_of_mass_proc(struct aylp_device *self, struct aylp_state *state);
// multithreaded version of proc function
int center_of_mass_proc_threaded(
	struct aylp_device *self, struct aylp_state *state
);
// single-region version whose window follows the previous center of mass
int center_of_mass_proc_track(
	struct aylp_device *self, struct aylp_state *state
);

// close center_of_mass device when loop exits
int center_of_mass_fini(struct aylp_device *self);
int center_of_mass_fini_threaded(struct aylp_device *self);

#endif
