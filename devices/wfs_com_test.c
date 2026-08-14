// Regression test: quantifies wfs_com against a reimplementation of
// center_of_mass's plain first-moment tracking (same math as
// center_of_mass_proc_track in devices/center_of_mass.c) on several
// synthetic scenarios, and validates the row_time rolling-shutter
// correction. Fails (nonzero exit) if wfs_com doesn't outperform the plain
// first moment on the scenarios it's specifically meant to help with, or if
// the rolling-shutter correction doesn't reduce error relative to leaving it
// off.
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <gsl/gsl_matrix.h>

#include "anyloop.h"
#include "wfs_com.h"
#include "xalloc.h"

static int failures = 0;
#define CHECK(cond, ...) do { \
	if (!(cond)) { \
		printf("FAIL: " __VA_ARGS__); \
		printf("\n"); \
		failures++; \
	} \
} while (0)

// ---- synthetic scene rendering -------------------------------------------

static void add_gaussian(gsl_matrix_uchar *img, double cy, double cx,
	double sigma, double amp
) {
	for (size_t i = 0; i < img->size1; i++) {
		for (size_t j = 0; j < img->size2; j++) {
			double dy = (double)i - cy, dx = (double)j - cx;
			double v = amp * exp(-(dy*dy+dx*dx)/(2*sigma*sigma));
			double cur = img->data[i*img->tda+j];
			double sum = cur + v;
			img->data[i*img->tda+j] =
				(unsigned char)(sum > 255.0 ? 255.0 : sum);
		}
	}
}

static void clear_img(gsl_matrix_uchar *img) { gsl_matrix_uchar_set_zero(img); }

// zero out a rectangular region (simulates partial beam clipping / an
// obstruction transiting the beam)
static void occlude(gsl_matrix_uchar *img, size_t y0, size_t y1,
	size_t x0, size_t x1
) {
	for (size_t i = y0; i < y1 && i < img->size1; i++)
		for (size_t j = x0; j < x1 && j < img->size2; j++)
			img->data[i*img->tda+j] = 0;
}

// ---- plain first-moment baseline (matches center_of_mass_proc_track) -----

struct plain_com_state {
	size_t win_y, win_x;	// window centre, image coords
	unsigned char threshold;
};

static void plain_com_step(struct plain_com_state *s, gsl_matrix_uchar *img,
	size_t win_h, size_t win_w, double *out_y, double *out_x
) {
	size_t half_y = win_h/2, half_x = win_w/2;
	if (s->win_y < half_y) s->win_y = half_y;
	if (s->win_x < half_x) s->win_x = half_x;
	if (s->win_y > img->size1 - win_h + half_y) s->win_y = img->size1 - win_h + half_y;
	if (s->win_x > img->size2 - win_w + half_x) s->win_x = img->size2 - win_w + half_x;
	size_t org_y = s->win_y - half_y, org_x = s->win_x - half_x;

	double y = 0.0, x = 0.0, sum = 0.0;
	for (size_t i = 0; i < win_h; i++) {
		for (size_t j = 0; j < win_w; j++) {
			unsigned char raw = img->data[(org_y+i)*img->tda + org_x+j];
			unsigned char el = raw > s->threshold ? raw - s->threshold : 0;
			y += (org_y+i)*el;
			x += (org_x+j)*el;
			sum += el;
		}
	}
	if (sum <= 0.0) { *out_y = -1; *out_x = -1; return; }	// held/lost, caller skips
	double abs_y = y/sum, abs_x = x/sum;
	s->win_y = (size_t)llround(abs_y);
	s->win_x = (size_t)llround(abs_x);
	*out_y = abs_y;
	*out_x = abs_x;
}

// ---- wfs_com harness -------------------------------------------------------

static void wfs_setup(struct aylp_wfs_com_data *data, size_t init_y,
	size_t init_x, double row_time
) {
	memset(data, 0, sizeof(*data));
	data->subap_h = 8; data->subap_w = 8;
	data->subap_rows = 4; data->subap_cols = 4;
	data->threshold = 5;
	data->search_radius = 3;
	data->ref_beta = 0.05;
	// Fast-then-slow template convergence (mirrors fsp.c's broad_mu_init/
	// broad_mu_tau): tau is scaled to this test's own short (paced, but
	// not real-time-realistic) run length, not a production value.
	data->ref_beta_init = 0.35;
	data->ref_beta_tau = 0.006;
	data->min_confidence = 0.4;
	data->flux_floor = 15.0;
	data->min_valid_subaps = 4;
	data->row_time = row_time;
	data->rolling_shutter = row_time > 0.0;
	data->init_y = (long)init_y;
	data->init_x = (long)init_x;
	data->reacquire_after = 10;
	data->n_subaps = data->subap_rows * data->subap_cols;
	data->ext_h = data->subap_h + 2*data->search_radius;
	data->ext_w = data->subap_w + 2*data->search_radius;
	data->win_h = data->subap_rows*data->subap_h + 2*data->search_radius;
	data->win_w = data->subap_cols*data->subap_w + 2*data->search_radius;
	data->ref = xcalloc(data->n_subaps * data->subap_h * data->subap_w, sizeof(double));
	data->ref_backup = xcalloc(data->n_subaps * data->subap_h * data->subap_w, sizeof(double));
	data->ref_set = xcalloc(data->n_subaps, sizeof(bool));
	data->ref_set_backup = xcalloc(data->n_subaps, sizeof(bool));
	data->matched = xcalloc(data->n_subaps, sizeof(bool));
	data->sub_dy = xcalloc(data->n_subaps, sizeof(double));
	data->sub_dx = xcalloc(data->n_subaps, sizeof(double));
	data->ext = xcalloc(data->ext_h * data->ext_w, sizeof(double));
	data->row_sum_w = xcalloc(data->subap_rows, sizeof(double));
	data->row_sum_wt = xcalloc(data->subap_rows, sizeof(double));
	data->row_sum_wdy = xcalloc(data->subap_rows, sizeof(double));
	data->row_sum_wdx = xcalloc(data->subap_rows, sizeof(double));
	data->com = xmalloc_type(gsl_vector, 2);
}

static void wfs_teardown(struct aylp_wfs_com_data *data)
{
	xfree(data->ref); xfree(data->ref_backup); xfree(data->ref_set); xfree(data->ref_set_backup); xfree(data->matched); xfree(data->sub_dy); xfree(data->sub_dx); xfree(data->ext);
	xfree(data->row_sum_w); xfree(data->row_sum_wt); xfree(data->row_sum_wdy); xfree(data->row_sum_wdx);
	xfree_type(gsl_vector, data->com);
}

// ---- scenario A/B/C/D: distortion robustness, vs plain first moment ------

typedef void (*render_fn)(gsl_matrix_uchar *img, double true_y, double true_x);

static double run_scenario(const char *name, render_fn render,
	size_t img_sz, int n_frames, double vy, double vx,
	double *baseline_rms, double *wfs_rms
) {
	gsl_matrix_uchar *img = gsl_matrix_uchar_alloc(img_sz, img_sz);
	double true_y = img_sz/2.0, true_x = img_sz/2.0;

	struct aylp_wfs_com_data wdata;
	wfs_setup(&wdata, (size_t)true_y, (size_t)true_x, 0.0);
	struct aylp_device self = {0};
	self.device_data = &wdata;
	self.type_out = AYLP_T_VECTOR;
	self.units_out = AYLP_U_MINMAX;
	struct aylp_state state = {0};

	struct plain_com_state pdata = { .win_y = (size_t)true_y,
		.win_x = (size_t)true_x, .threshold = 5 };
	size_t plain_win = wdata.subap_rows*wdata.subap_h;	// same core size

	double se_base = 0.0, se_wfs = 0.0;
	int n_scored = 0;
	for (int f = 0; f < n_frames; f++) {
		true_y += vy; true_x += vx;
		render(img, true_y, true_x);
		state.matrix_uchar = img;
		// ref_beta_init/tau decays against real wall-clock time (see
		// wfs_com.h), so a back-to-back call loop would never let it
		// decay at all -- pace frames in real time to actually exercise
		// the schedule, same as a live camera loop would.
		usleep(1000);
		wfs_com_proc(&self, &state);
		double py, px;
		plain_com_step(&pdata, img, plain_win, plain_win, &py, &px);
		if (f < 8) continue;	// skip acquisition/warm-up transient
		double wy = (wdata.last_y+1.0)/2.0*(img_sz-1);
		double wx = (wdata.last_x+1.0)/2.0*(img_sz-1);
		if (py >= 0) {
			double eb = hypot(py-true_y, px-true_x);
			se_base += eb*eb;
		}
		double ew = hypot(wy-true_y, wx-true_x);
		se_wfs += ew*ew;
		n_scored++;
	}
	*baseline_rms = sqrt(se_base/n_scored);
	*wfs_rms = sqrt(se_wfs/n_scored);
	printf("%-28s baseline RMS %.3f px   wfs_com RMS %.3f px\n",
		name, *baseline_rms, *wfs_rms);

	wfs_teardown(&wdata);
	gsl_matrix_uchar_free(img);
	return 0;
}

static void render_clean(gsl_matrix_uchar *img, double y, double x)
{ clear_img(img); add_gaussian(img, y, x, 2.5, 200.0); }

static void render_asymmetric(gsl_matrix_uchar *img, double y, double x)
{
	clear_img(img);
	add_gaussian(img, y, x, 2.5, 200.0);
	// off-centre "hot spot" (scintillation / partial reflection), enough
	// energy to pull a naive first moment noticeably off the true centre
	add_gaussian(img, y-6.0, x+5.0, 1.5, 150.0);
}

static void render_occluded(gsl_matrix_uchar *img, double y, double x)
{
	clear_img(img);
	add_gaussian(img, y, x, 2.5, 200.0);
	// clip the bottom-right quadrant of the spot (partial beam clipping)
	occlude(img, (size_t)y, (size_t)y+6, (size_t)x, (size_t)x+6);
}

// A broad spot (wide enough that every subaperture row group carries flux)
// whose x-centre ramps linearly with row index -- a static, deterministic
// stand-in for rolling-shutter skew. `slope` px of x-offset per row; 0
// renders the unskewed reference shape. Used by the interleave test, which
// needs shear that RELIABLY exceeds max_row_shear while still being small
// enough that every subaperture matches confidently within search_radius
// (a sinusoidal scene produces shear that varies frame to frame, which
// muddles "rejected for shear" with "dropped for low confidence").
static void render_broad_skewed(gsl_matrix_uchar *img, double y0, double x0,
	double sigma, double amp, double slope
) {
	clear_img(img);
	for (size_t row = 0; row < img->size1; row++) {
		double cx = x0 + slope*((double)row - y0);
		for (size_t col = 0; col < img->size2; col++) {
			double dy = (double)row - y0, dx = (double)col - cx;
			double v = amp * exp(-(dy*dy+dx*dx)/(2*sigma*sigma));
			img->data[row*img->tda+col] =
				(unsigned char)(v > 255.0 ? 255.0 : v);
		}
	}
}

static void render_stray(gsl_matrix_uchar *img, double y, double x)
{
	clear_img(img);
	add_gaussian(img, y, x, 2.5, 200.0);
	// a STATIC stray reflection elsewhere in the window, same rationale
	// as doc/devices/center_of_mass.md's tracking-window section
	add_gaussian(img, 20.0, 44.0, 2.0, 180.0);
}

// Frames from acquisition until tracking error drops below thresh_px and
// STAYS below it for `hold` consecutive frames (a single lucky dip doesn't
// count as "converged"). Real-time-paced at ~1 kHz so ref_beta_init/tau
// (which decays against wall-clock time) has something to act against.
// Returns -1 if it never converges within n_frames.
static int measure_convergence(bool use_schedule, double thresh_px, int hold)
{
	const size_t IMG = 64;
	const int N_FRAMES = 60;
	gsl_matrix_uchar *img = gsl_matrix_uchar_alloc(IMG, IMG);
	struct aylp_wfs_com_data data;
	wfs_setup(&data, IMG/2, IMG/2, 0.0);
	if (!use_schedule) data.ref_beta_init = 0.0;	// constant ref_beta only
	struct aylp_device self = {0};
	self.device_data = &data;
	self.type_out = AYLP_T_VECTOR;
	self.units_out = AYLP_U_MINMAX;
	struct aylp_state st = {0};

	double true_y = IMG/2.0, true_x = IMG/2.0;
	int converge_frame = -1, hold_count = 0;
	for (int f = 0; f < N_FRAMES && converge_frame < 0; f++) {
		true_y += 0.06; true_x += 0.04;
		render_clean(img, true_y, true_x);
		st.matrix_uchar = img;
		usleep(1000);
		wfs_com_proc(&self, &st);
		double ry = (data.last_y+1.0)/2.0*(IMG-1);
		double rx = (data.last_x+1.0)/2.0*(IMG-1);
		double err = hypot(ry-true_y, rx-true_x);
		if (f >= 1 && err < thresh_px) {
			if (++hold_count >= hold) converge_frame = f - hold + 1;
		} else {
			hold_count = 0;
		}
	}
	wfs_teardown(&data);
	gsl_matrix_uchar_free(img);
	return converge_frame;
}

int main(void)
{
	printf("=== scenario comparison: wfs_com vs plain first-moment ===\n");
	double br, wr;

	run_scenario("clean symmetric blob", render_clean, 64, 60, 0.06, 0.04, &br, &wr);
	// A noiseless symmetric Gaussian's exact first moment IS its true
	// centre, so the plain baseline is near machine precision here (this
	// is its best case, not a representative one). wfs_com's correlation
	// + parabolic sub-pixel fit carries a known small inherent bias even
	// with nothing wrong in the scene -- compare against an absolute
	// bound, not against the near-zero baseline, so this check reflects
	// "is wfs_com's clean-scene overhead acceptable" rather than
	// "did the baseline's best-case number regress".
	CHECK(wr < 1.0, "wfs_com's clean-scene error should stay under 1 px "
		"(baseline %.3f, wfs_com %.3f)", br, wr);

	run_scenario("asymmetric hot spot", render_asymmetric, 64, 60, 0.06, 0.04, &br, &wr);
	CHECK(wr < br, "wfs_com should beat plain first-moment under "
		"asymmetric distortion (baseline %.3f, wfs_com %.3f)", br, wr);

	run_scenario("partial occlusion", render_occluded, 64, 60, 0.06, 0.04, &br, &wr);
	CHECK(wr < br, "wfs_com should beat plain first-moment under partial "
		"occlusion (baseline %.3f, wfs_com %.3f)", br, wr);

	run_scenario("static stray reflection", render_stray, 64, 60, 0.06, 0.04, &br, &wr);
	CHECK(wr < br, "wfs_com should beat plain first-moment with a stray "
		"reflection in-window (baseline %.3f, wfs_com %.3f)", br, wr);

	// ---- catastrophic apparent 10 px jump: reject and hold ----
	// Models the operational failure this tracker must contain: one rolling-
	// shutter-corrupted frame appears displaced far beyond the normal per-frame
	// motion. Boundary-peak rejection must not quantize that into repeated
	// search-radius steps and walk the control loop toward the artifact.
	printf("\n=== catastrophic 10 px apparent jump ===\n");
	{
		const size_t IMG = 32;
		gsl_matrix_uchar *img = gsl_matrix_uchar_alloc(IMG, IMG);
		json_object *params = json_tokener_parse(
			"{\"subap_height\":6,\"subap_width\":6,"
			"\"subap_rows\":3,\"subap_cols\":3,\"search_radius\":2,"
			"\"reject_edge_matches\":true,\"threshold\":5,"
			"\"ref_beta\":0.01,\"ref_beta_init\":0.25,"
			"\"ref_beta_tau\":2.0,\"min_confidence\":0.55,"
			"\"flux_floor\":20,\"min_valid_subaps\":6,"
			"\"rolling_shutter\":true,\"max_row_shear\":2.5,"
			"\"init_y\":16,\"init_x\":16,\"reacquire_after\":10}");
		struct aylp_device self_j = {0};
		self_j.params = params;
		CHECK(wfs_com_init(&self_j) == 0,
			"production 32x32 WFS configuration should initialize");
		struct aylp_wfs_com_data *d = self_j.device_data;
		struct aylp_state st = {0};
		for (int f = 0; f < 6; f++) {
			render_clean(img, IMG/2.0, IMG/2.0);
			st.matrix_uchar = img;
			st.header.status = 0;
			usleep(1000);
			wfs_com_proc(&self_j, &st);
		}
		// A late-seeding subaperture on a frame that fails the frame-level
		// count gate must roll back both its numeric template and ref_set bit.
		size_t late = 4;
		d->ref_set[late] = false;
		memset(d->ref + late*d->subap_h*d->subap_w, 0,
			d->subap_h*d->subap_w*sizeof(double));
		d->min_valid_subaps = d->n_subaps;
		render_clean(img, IMG/2.0, IMG/2.0);
		st.matrix_uchar = img;
		st.header.status = 0;
		wfs_com_proc(&self_j, &st);
		CHECK((st.header.status & AYLP_FRAME_REJECTED) && !d->ref_set[late],
			"rejected late template seed should roll back ref_set");
		d->min_valid_subaps = 6;
		st.matrix_uchar = img;
		st.header.status = 0;
		wfs_com_proc(&self_j, &st);
		// Persistent subpixel residual: it must be reconstructed afresh, not
		// accumulated once per frame into a fictitious drift.
		for (int f = 0; f < 20; f++) {
			render_clean(img, IMG/2.0, IMG/2.0 + 0.4);
			st.matrix_uchar = img;
			st.header.status = 0;
			wfs_com_proc(&self_j, &st);
		}
		double measured_x = (d->last_x+1.0)/2.0*(IMG-1);
		CHECK(!(st.header.status & AYLP_FRAME_REJECTED),
			"ordinary persistent +0.4 px motion should remain valid");
		CHECK(fabs(measured_x-(IMG/2.0+0.4)) < 0.35,
			"persistent +0.4 px residual should remain near +0.4, got %.3f px",
			measured_x-IMG/2.0);
		double y_hold = d->last_y, x_hold = d->last_x;
		render_clean(img, IMG/2.0, IMG/2.0 - 10.0);
		st.matrix_uchar = img;
		st.header.status = 0;
		usleep(1000);
		wfs_com_proc(&self_j, &st);
		CHECK(st.header.status & AYLP_FRAME_REJECTED,
			"a 10 px out-of-range apparent jump should be rejected");
		CHECK(d->last_y == y_hold && d->last_x == x_hold,
			"a rejected 10 px apparent jump should hold output exactly");
		printf("10 px apparent jump rejected; output held exactly\n");
		wfs_com_fini(&self_j);
		json_object_put(params);
		gsl_matrix_uchar_free(img);
	}

	// ---- rolling-shutter correction: corrected must beat uncorrected ----
	printf("\n=== rolling-shutter correction ===\n");
	{
		const size_t IMG = 96;
		const double FRAME_PERIOD = 1.0/3788.0;	// s, this rig's fs
		// Self-consistent: reading out all IMG rows takes exactly one
		// frame period (no blanking, for simplicity).
		const double ROW_TIME_TRUE = FRAME_PERIOD / (double)IMG;
		// Bounded, oscillatory true motion -- a runaway constant
		// velocity (the first version of this test) walks the target
		// off the image within a few dozen frames, which produced a
		// meaningless 170+ px "error" that reflected the test's own
		// bug, not wfs_com's behavior. Amplitude/frequency are chosen
		// so the intra-window (rolling-shutter) skew is a clear ~1-2
		// px over one frame's readout -- deliberately larger than this
		// rig's real vibration lines, to make the effect unambiguous
		// in a short run, not a claim about real motion amplitude.
		const double AY = 10.0, FY = 250.0;		// px, Hz
		const double AX = 4.0,  FX = 180.0, PHX = 1.0;	// px, Hz, rad
		gsl_matrix_uchar *img = gsl_matrix_uchar_alloc(IMG, IMG);

		struct aylp_wfs_com_data unc, cor;
		wfs_setup(&unc, IMG/2, IMG/2, 0.0);
		wfs_setup(&cor, IMG/2, IMG/2, ROW_TIME_TRUE);
		struct aylp_device self_u = {0}, self_c = {0};
		self_u.device_data = &unc; self_u.type_out = AYLP_T_VECTOR; self_u.units_out = AYLP_U_MINMAX;
		self_c.device_data = &cor; self_c.type_out = AYLP_T_VECTOR; self_c.units_out = AYLP_U_MINMAX;
		struct aylp_state st = {0};

		double y0 = IMG/2.0, x0 = IMG/2.0, t0 = 0.0;
		double bu_y=0, bu_x=0, bc_y=0, bc_x=0;
		double su_y=0, su_x=0, sc_y=0, sc_x=0;
		int n_scored = 0;
		for (int f = 0; f < 80; f++) {
			clear_img(img);
			for (size_t row = 0; row < IMG; row++) {
				double t = t0 + ROW_TIME_TRUE*(double)row;
				double cy = y0 + AY*sin(2*M_PI*FY*t);
				double cx = x0 + AX*sin(2*M_PI*FX*t + PHX);
				for (size_t col = 0; col < IMG; col++) {
					double dy = (double)row-cy, dx = (double)col-cx;
					double v = 200.0*exp(-(dy*dy+dx*dx)/(2*2.5*2.5));
					img->data[row*img->tda+col] = (unsigned char)(v > 255 ? 255 : v);
				}
			}
			bool score = f >= 8;
			// Capture each instance's window position as it stood
			// BEFORE this call -- that is what was actually used to
			// extract this frame's subapertures, and hence the frame
			// of reference each one's output should be scored against.
			size_t org_y_u = score ? unc.win_y - unc.win_h/2 : 0;
			size_t org_y_c = score ? cor.win_y - cor.win_h/2 : 0;
			double t0_frame = t0;

			st.matrix_uchar = img;
			usleep(1000);	// pace in real time; see the note above
			wfs_com_proc(&self_u, &st);
			st.matrix_uchar = img;
			wfs_com_proc(&self_c, &st);
			t0 += FRAME_PERIOD;

			if (score) {
				// Ground truth for each instance is evaluated at ITS
				// OWN window-centre reference time -- the exact
				// instant its output is meant to represent (for the
				// corrected instance, the regression target; for the
				// uncorrected one, the same instant, so its deviation
				// from it is exactly the rolling-shutter smear this
				// test exists to reveal).
				double tru = t0_frame + ROW_TIME_TRUE
					*((double)org_y_u + unc.subap_rows*unc.subap_h/2.0);
				double trc = t0_frame + ROW_TIME_TRUE
					*((double)org_y_c + cor.subap_rows*cor.subap_h/2.0);
				double tuy = y0+AY*sin(2*M_PI*FY*tru);
				double tux = x0+AX*sin(2*M_PI*FX*tru+PHX);
				double tcy = y0+AY*sin(2*M_PI*FY*trc);
				double tcx = x0+AX*sin(2*M_PI*FX*trc+PHX);
				double uy = (unc.last_y+1.0)/2.0*(IMG-1);
				double ux = (unc.last_x+1.0)/2.0*(IMG-1);
				double cyv = (cor.last_y+1.0)/2.0*(IMG-1);
				double cxv = (cor.last_x+1.0)/2.0*(IMG-1);
				// Accumulate per-axis error and its MEAN separately so
				// the RMS can be debiased below. Without that, this
				// number is dominated by a constant offset that has
				// nothing to do with rolling shutter (see below).
				bu_y += uy-tuy; bu_x += ux-tux;
				bc_y += cyv-tcy; bc_x += cxv-tcx;
				su_y += (uy-tuy)*(uy-tuy); su_x += (ux-tux)*(ux-tux);
				sc_y += (cyv-tcy)*(cyv-tcy); sc_x += (cxv-tcx)*(cxv-tcx);
				n_scored++;
			}
		}
		// DEBIASED RMS. The raw RMS here is dominated by a constant
		// offset: init_y/init_x define zero error, so whatever
		// displacement the beam happens to have at the acquisition
		// instant is absorbed and shows up as a fixed bias on every
		// later frame (this scene's x phase puts it ~1.3 px off). That
		// bias swamps the shear term and made the correction look like
		// it delivered ~13% when measured against the raw total. Only
		// the varying part says anything about rolling shutter.
		double n = n_scored;
		double rms_unc = sqrt(su_y/n - (bu_y/n)*(bu_y/n)
			+ su_x/n - (bu_x/n)*(bu_x/n));
		double rms_cor = sqrt(sc_y/n - (bc_y/n)*(bc_y/n)
			+ sc_x/n - (bc_x/n)*(bc_x/n));
		printf("uncorrected (row_time=0) RMS: %.4f px (debiased)\n", rms_unc);
		printf("corrected   (row_time=true) RMS: %.4f px (debiased)\n", rms_cor);
		// The correction is NOT asserted to help. Measured across
		// amplitude/frequency it does not: with only subap_rows row
		// groups, of which a real beam lights 2-3, the fitted intra-frame
		// slope carries little signal and a lot of noise. What must hold
		// is that it never does appreciable HARM -- that is what the
		// significance shrinkage in wfs_com.c guarantees, and what
		// regressed badly before it existed (0.251 -> 0.519 px here).
		CHECK(rms_cor < 1.25*rms_unc, "row_time correction must not "
			"appreciably increase tracking error (uncorrected %.4f, "
			"corrected %.4f)", rms_unc, rms_cor);

		wfs_teardown(&unc);
		wfs_teardown(&cor);
		gsl_matrix_uchar_free(img);
	}

	// ---- rolling-shutter REJECTION: detects shear, holds output ----
	// Same synthetic rolling-shutter scene as above, but tests
	// max_row_shear independently of row_time (row_time left at 0): every
	// frame with detectable row-group disagreement should be flagged
	// AYLP_FRAME_REJECTED and hold the previous frame's output exactly,
	// rather than reporting a shear-distorted position.
	printf("\n=== rolling-shutter shear rejection (max_row_shear) ===\n");
	{
		const size_t IMG = 96;
		const double FRAME_PERIOD = 1.0/3788.0;
		const double ROW_TIME_TRUE = FRAME_PERIOD / (double)IMG;
		const double AY = 10.0, FY = 250.0;
		const double AX = 4.0,  FX = 180.0, PHX = 1.0;
		gsl_matrix_uchar *img = gsl_matrix_uchar_alloc(IMG, IMG);

		struct aylp_wfs_com_data rej;
		wfs_setup(&rej, IMG/2, IMG/2, 0.0);
		// px; well under this scenario's several-px intra-window skew
		rej.max_row_shear = 1.5;
		struct aylp_device self_r = {0};
		self_r.device_data = &rej;
		self_r.type_out = AYLP_T_VECTOR;
		self_r.units_out = AYLP_U_MINMAX;
		struct aylp_state st = {0};

		double y0 = IMG/2.0, x0 = IMG/2.0, t0 = 0.0;
		int n_scored = 0, n_rejected = 0, n_held_correctly = 0;
		int n_shear_rejected = 0;
		double prev_y = 0.0, prev_x = 0.0;
		size_t ref_bytes = rej.n_subaps * rej.subap_h * rej.subap_w
			* sizeof(double);
		double *ref_before = xmalloc(ref_bytes);
		int n_refs_held = 0;
		for (int f = 0; f < 80; f++) {
			clear_img(img);
			for (size_t row = 0; row < IMG; row++) {
				double t = t0 + ROW_TIME_TRUE*(double)row;
				double cy = y0 + AY*sin(2*M_PI*FY*t);
				double cx = x0 + AX*sin(2*M_PI*FX*t + PHX);
				for (size_t col = 0; col < IMG; col++) {
					double dy = (double)row-cy, dx = (double)col-cx;
					double v = 200.0*exp(-(dy*dy+dx*dx)/(2*2.5*2.5));
					img->data[row*img->tda+col] = (unsigned char)(v > 255 ? 255 : v);
				}
			}
			st.matrix_uchar = img;
			st.header.status = 0;
			size_t shear_before = rej.shear_rejected;
			memcpy(ref_before, rej.ref, ref_bytes);
			usleep(1000);	// pace in real time; see the note above
			wfs_com_proc(&self_r, &st);
			bool rejected = st.header.status & AYLP_FRAME_REJECTED;
			if (rej.shear_rejected > shear_before) {
				n_shear_rejected++;
				if (!memcmp(ref_before, rej.ref, ref_bytes)) n_refs_held++;
			}
			if (f >= 8) {
				n_scored++;
				if (rejected) {
					n_rejected++;
					if (rej.last_y == prev_y && rej.last_x == prev_x)
						n_held_correctly++;
				}
			}
			prev_y = rej.last_y; prev_x = rej.last_x;
			t0 += FRAME_PERIOD;
		}
		printf("rejected %d/%d scored frames for shear\n", n_rejected, n_scored);
		CHECK(n_rejected > 0, "max_row_shear should reject at least some "
			"frames on this deliberately-skewed scenario (rejected "
			"0/%d)", n_scored);
		CHECK(rej.shear_rejected == (size_t)n_shear_rejected,
			"shear_rejected counter (%zu) should match frames flagged "
			"specifically by the shear gate (%d)", rej.shear_rejected,
			n_shear_rejected);
		CHECK(n_held_correctly == n_rejected, "every rejected frame should "
			"hold the previous frame's output unchanged (%d/%d did)",
			n_held_correctly, n_rejected);
		CHECK(n_refs_held == n_shear_rejected, "every shear-rejected frame should "
			"leave reference templates unchanged (%d/%d did)", n_refs_held,
			n_shear_rejected);

		xfree(ref_before);
		wfs_teardown(&rej);
		gsl_matrix_uchar_free(img);
	}

	// ---- interleaved low-signal / shear-rejected frames ----
	// A shear-rejected frame only reaches the rejection branch by PASSING
	// the min_valid_subaps gate -- i.e. the beam was demonstrably present
	// and well-matched that frame. So it must reset the consecutive-loss
	// counter, exactly as a normally-tracked frame does. If it merely
	// leaves the counter frozen, low-signal frames interleaved with
	// shear-rejected ones accumulate toward reacquire_after across frames
	// where the beam was repeatedly seen just fine, and eventually fire a
	// spurious AYLP_BEAM_LOST + window/template reset.
	//
	// This alternates dark frames (no signal -> increments the counter)
	// with high-shear frames (good signal, rejected for shear). With a
	// correct reset, the counter never exceeds 1 and no loss is ever
	// asserted; without it, the dark frames alone reach reacquire_after.
	printf("\n=== interleaved low-signal / shear rejection ===\n");
	{
		const size_t IMG = 96;
		const double SIGMA = 8.0, AMP = 200.0;
		// ~2.9 px of measured row-group spread across the 32-row window:
		// comfortably above max_row_shear (1.5) but with per-row-group
		// offsets (~+-1.9 px) still well inside search_radius (3), so
		// these frames are rejected for SHEAR, never dropped for low
		// confidence. Note the measured spread runs ~75% of the rendered
		// skew -- the correlation estimator under-responds slightly on a
		// blob this smooth -- so this is calibrated from the observed
		// number, not from the render geometry.
		const double SLOPE = 5.0/32.0;
		const int N_SEED = 12;		// clean frames to establish templates
		const int N_ALTERNATE = 40;	// >= 2*reacquire_after dark frames
		gsl_matrix_uchar *img = gsl_matrix_uchar_alloc(IMG, IMG);

		struct aylp_wfs_com_data d;
		wfs_setup(&d, IMG/2, IMG/2, 0.0);
		d.max_row_shear = 1.5;
		struct aylp_device self_i = {0};
		self_i.device_data = &d;
		self_i.type_out = AYLP_T_VECTOR;
		self_i.units_out = AYLP_U_MINMAX;
		struct aylp_state st = {0};

		double y0 = IMG/2.0, x0 = IMG/2.0;

		// seed templates on the same broad spot, unskewed
		for (int f = 0; f < N_SEED; f++) {
			render_broad_skewed(img, y0, x0, SIGMA, AMP, 0.0);
			st.matrix_uchar = img;
			st.header.status = 0;
			usleep(1000);
			wfs_com_proc(&self_i, &st);
		}

		int n_dark = 0, n_shear = 0, n_other = 0;
		size_t max_lost = 0;
		bool beam_lost_ever = false;
		for (int f = 0; f < N_ALTERNATE; f++) {
			bool dark = (f % 2) == 1;
			if (dark) {
				// entirely below flux_floor -> no valid subapertures
				clear_img(img);
			} else {
				// Good signal, sheared across rows. The sign alternates
				// so this exercises both signs of the rejection gate.
				double slope = ((f/2) % 2 == 0) ? SLOPE : -SLOPE;
				render_broad_skewed(img, y0, x0, SIGMA, AMP, slope);
			}
			st.matrix_uchar = img;
			st.header.status = 0;
			usleep(1000);
			wfs_com_proc(&self_i, &st);

			if (st.header.status & AYLP_BEAM_LOST) beam_lost_ever = true;
			if (d.lost > max_lost) max_lost = d.lost;
			if (dark) n_dark++;
			else if (st.header.status & AYLP_FRAME_REJECTED) n_shear++;
			else {
				n_other++;
				// Report the spread this frame actually produced --
				// the scenario is only valid if the rendered skew
				// really does exceed max_row_shear, and guessing at
				// that from the render geometry alone has already
				// been wrong once.
				if (n_other <= 2) {
					double ymin = 0, ymax = 0, xmin = 0, xmax = 0;
					size_t seen = 0;
					for (size_t r = 0; r < d.subap_rows; r++) {
						double w = d.row_sum_w[r];
						if (w <= 0.0) continue;
						double yr = d.row_sum_wdy[r]/w;
						double xr = d.row_sum_wdx[r]/w;
						if (!seen || yr < ymin) ymin = yr;
						if (!seen || yr > ymax) ymax = yr;
						if (!seen || xr < xmin) xmin = xr;
						if (!seen || xr > xmax) xmax = xr;
						seen++;
					}
					printf("  (unrejected frame %d: %zu row groups, "
						"spread y=%.2f x=%.2f px vs threshold %.2f)\n",
						f, seen, ymax-ymin, xmax-xmin, d.max_row_shear);
				}
			}
		}
		printf("%d dark frames, %d shear-rejected, %d normally tracked; "
			"peak consecutive-loss counter %zu\n",
			n_dark, n_shear, n_other, max_lost);

		// The scenario is only meaningful if it actually produced the
		// interleaving it intends to -- otherwise the assertions below
		// would pass vacuously.
		CHECK(n_other == 0 && n_shear == N_ALTERNATE/2, "every non-dark "
			"frame should be rejected for SHEAR specifically, not "
			"dropped for low confidence -- otherwise this tests the "
			"ordinary lost-signal path rather than the shear/loss "
			"interaction (%d shear-rejected, %d normally tracked, of %d)",
			n_shear, n_other, N_ALTERNATE/2);
		CHECK(n_dark > (int)d.reacquire_after, "interleave scenario should "
			"produce more dark frames (%d) than reacquire_after (%zu), "
			"so a missing reset would actually fire", n_dark,
			d.reacquire_after);

		CHECK(!beam_lost_ever, "a shear-rejected frame proves the beam is "
			"present, so interleaved dark frames must not accumulate "
			"into a spurious AYLP_BEAM_LOST");
		CHECK(max_lost <= 1, "the consecutive-loss counter should reset on "
			"every shear-rejected frame, so it should never exceed 1 "
			"in this alternating scenario (peaked at %zu)", max_lost);

		wfs_teardown(&d);
		gsl_matrix_uchar_free(img);
	}

	// ---- template convergence SPEED: schedule vs fixed ref_beta ----
	// The scenario RMS numbers above are steady-state metrics and
	// deliberately exclude the first several frames as warm-up -- which
	// hides exactly what a fast-then-slow schedule is supposed to
	// improve. Measure it directly instead: frames from acquisition
	// until error first drops below and stays below a threshold.
	printf("\n=== template convergence speed (schedule vs fixed ref_beta) ===\n");
	{
		int f_fixed = measure_convergence(false, 0.6, 5);
		int f_sched = measure_convergence(true, 0.6, 5);
		printf("fixed ref_beta:     converged at frame %d\n", f_fixed);
		printf("ref_beta schedule:  converged at frame %d\n", f_sched);
		CHECK(f_sched >= 0, "scheduled ref_beta should converge within "
			"the test window");
		if (f_sched >= 0 && f_fixed >= 0) {
			CHECK(f_sched <= f_fixed, "the schedule should converge "
				"no slower than a fixed ref_beta (fixed: frame "
				"%d, scheduled: frame %d)", f_fixed, f_sched);
		}
	}

	// ---- static-scene invariant: no self-generated drift ----
	// The strongest statement that can be made about a tracker: feed it the
	// SAME BYTES every frame and its output must not move. Any motion is
	// manufactured by the device's own template feedback. This is the
	// scenario that caught per-subaperture template registration diverging:
	// outer subapertures see a smooth gradient, score ~1 NCC at every offset,
	// and used to ratchet their own templates outward on their own noise
	// until their peaks pinned to the search boundary -- at which point the
	// frame validity gate failed on EVERY frame and the device asserted
	// AYLP_BEAM_LOST and re-acquired roughly every 34 frames, on a
	// motionless beam. None of the moving scenarios above detect this,
	// because there real motion masks the fabricated part.
	printf("\n=== static-scene invariant (no self-generated drift) ===\n");
	{
		const size_t IMG = 32;
		const int N = 400;
		gsl_matrix_uchar *img = gsl_matrix_uchar_alloc(IMG, IMG);
		// deliberately off the integer window centre: a sub-pixel offset is
		// what the sub-pixel refinement feeds back, so 16.4 exercises the
		// failure that an exactly-centred beam would hide
		render_clean(img, IMG/2.0, IMG/2.0 + 0.4);
		json_object *params = json_tokener_parse(
			"{\"subap_height\":6,\"subap_width\":6,"
			"\"subap_rows\":3,\"subap_cols\":3,\"search_radius\":2,"
			"\"reject_edge_matches\":true,\"threshold\":5,"
			"\"ref_beta\":0.01,\"ref_beta_init\":0.25,"
			"\"ref_beta_tau\":2.0,\"min_confidence\":0.55,"
			"\"flux_floor\":20,\"min_valid_subaps\":6,"
			"\"rolling_shutter\":true,\"max_row_shear\":2.5,"
			"\"init_y\":16,\"init_x\":16,\"reacquire_after\":10}");
		struct aylp_device self_s = {0};
		self_s.params = params;
		CHECK(wfs_com_init(&self_s) == 0,
			"production 32x32 config should initialize");
		struct aylp_wfs_com_data *d = self_s.device_data;
		struct aylp_state st = {0};
		int rejected = 0;
		double first = 0.0, worst = 0.0;
		for (int f = 0; f < N; f++) {
			st.matrix_uchar = img;
			st.header.status = 0;
			wfs_com_proc(&self_s, &st);
			if (st.header.status & AYLP_FRAME_REJECTED) { rejected++; continue; }
			double x = (d->last_x+1.0)/2.0*(IMG-1);
			if (!first) first = x;
			if (fabs(x-first) > worst) worst = fabs(x-first);
		}
		printf("%d/%d frames rejected; peak excursion %.3f px\n",
			rejected, N, worst);
		// Only the warm-up frame may be rejected: a motionless, bright,
		// well-matched beam is the easiest case there is.
		CHECK(rejected <= 1, "a static beam should not be rejected after "
			"warm-up (%d of %d frames rejected)", rejected, N);
		CHECK(!(st.header.status & AYLP_BEAM_LOST),
			"a static beam should never be reported lost");
		// Bounded, converging settle -- not the unbounded ramp that the old
		// per-subaperture registration produced (which ran to the search
		// boundary, ~2 px, and then broke the loop entirely).
		CHECK(worst < 1.0, "self-generated drift on an identical frame "
			"should stay sub-pixel (peak excursion %.3f px)", worst);
		wfs_com_fini(&self_s);
		json_object_put(params);
		gsl_matrix_uchar_free(img);
	}

	if (failures) {
		printf("\n%d CHECK(S) FAILED\n", failures);
		return 1;
	}
	printf("\nall checks passed\n");
	return 0;
}
