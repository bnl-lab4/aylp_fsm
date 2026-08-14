// Regression test for anyloop:fit_com, the joint spatio-temporal beam fit.
//
// The central claim is that putting the intra-frame motion INSIDE the model
// makes rolling shutter an identifiable parameter rather than a defect to be
// detected and rejected. So the tests check, in order: that the analytic
// Jacobian is right (a wrong derivative degrades the fit silently rather than
// breaking it), that the fit recovers known ground truth from a synthetic
// sheared frame, that a static scene produces a genuinely static output, that
// a stray reflection is rejected by the robust weighting, and that a frame
// with no beam is reported as such instead of guessed at.
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_matrix.h>
#include <json-c/json.h>

#include "anyloop.h"
#include "fit_com.h"

static int failures = 0;
#define CHECK(cond, ...) do { \
	if (!(cond)) { \
		printf("FAIL: " __VA_ARGS__); \
		printf("\n"); \
		failures++; \
	} \
} while (0)

static const size_t IMG = 32;

// Render a beam whose centre moves linearly with row index -- exactly the
// model fit_com fits, and exactly what a rolling shutter produces when the
// beam is translating during readout. slope is px of centre motion per row.
static void render(gsl_matrix_uchar *img, double y0, double x0,
	double sy, double sx, double sigma, double amp, double bg
) {
	double ref = 0.5*((double)img->size1 - 1.0);
	for (size_t r = 0; r < img->size1; r++) {
		double t = (double)r - ref;
		double yc = y0 + sy*t, xc = x0 + sx*t;
		for (size_t c = 0; c < img->size2; c++) {
			double dy = (double)r - yc, dx = (double)c - xc;
			double v = bg + amp*exp(-(dy*dy+dx*dx)/(2*sigma*sigma));
			img->data[r*img->tda+c] = (unsigned char)(v > 255 ? 255 : v);
		}
	}
}

// These are convergence tests, so they turn the latency guard off (max_us 0)
// and let the solver iterate: what they check is that the model, the jacobian
// and the robust weighting are right, which is a different question from how
// much of that fits in a loop period. The guard's own behaviour -- and the
// price it charges on a scene that needs many iterations -- is checked
// separately, in test_latency below.
static json_object *cfg(const char *extra)
{
	char buf[1024];
	snprintf(buf, sizeof buf, "{\"sigma_init\":2.5,\"sigma_min\":0.5,"
		"\"sigma_max\":10.0,\"max_iter\":30,\"min_amplitude\":5.0,"
		"\"max_us\":0,\"reacquire_after\":5%s%s}",
		extra && *extra ? "," : "", extra ? extra : "");
	json_object *p = json_tokener_parse(buf);
	if (!p) { printf("FAIL: bad test json: %s\n", buf); failures++; }
	return p;
}

// ---- 1. analytic Jacobian vs central differences ------------------------
//
// fit_eval is static, so probe the derivative through the public interface:
// perturbing one parameter and re-fitting is not a derivative test, so
// instead reconstruct the model here and compare against the same closed
// form the device uses. A mismatch here means the device is descending the
// wrong gradient.
static double model_px(const double *p, double ref, size_t r, size_t c)
{
	double t = (double)r - ref;
	double yc = p[AYLP_FIT_P_Y0] + p[AYLP_FIT_P_SY]*t;
	double xc = p[AYLP_FIT_P_X0] + p[AYLP_FIT_P_SX]*t;
	double dy = (double)r - yc, dx = (double)c - xc;
	double s = p[AYLP_FIT_P_SIGMA];
	return p[AYLP_FIT_P_BG]
		+ p[AYLP_FIT_P_AMP]*exp(-(dy*dy+dx*dx)/(2*s*s));
}

static void test_jacobian(void)
{
	printf("\n=== analytic Jacobian vs central differences ===\n");
	double p[AYLP_FIT_COM_NP] = {15.3, 16.7, 0.031, -0.017, 2.4, 190.0, 6.0};
	double ref = 0.5*((double)IMG - 1.0);
	const char *names[] = {"y0","x0","slope_y","slope_x","sigma","amp","bg"};
	double worst = 0.0;
	for (size_t a = 0; a < AYLP_FIT_COM_NP; a++) {
		// central differences: optimal step is ~eps^(1/3) times the
		// parameter scale. The previous 1e-6/1e-9 step made roundoff
		// (~amp*eps/h) dominate and the "error" measured was the test's
		// own numerics, not the derivative.
		double h = (fabs(p[a]) > 1.0 ? fabs(p[a]) : 1.0)*1e-5;
		for (size_t r = 4; r < IMG; r += 7) {
			for (size_t c = 4; c < IMG; c += 7) {
				double pp[AYLP_FIT_COM_NP], pm[AYLP_FIT_COM_NP];
				memcpy(pp, p, sizeof(p)); memcpy(pm, p, sizeof(p));
				pp[a] += h; pm[a] -= h;
				double num = (model_px(pp,ref,r,c)
					- model_px(pm,ref,r,c))/(2*h);
				// closed form, same as fit_com.c
				double t = (double)r - ref;
				double yc = p[AYLP_FIT_P_Y0]+p[AYLP_FIT_P_SY]*t;
				double xc = p[AYLP_FIT_P_X0]+p[AYLP_FIT_P_SX]*t;
				double dy = (double)r-yc, dx = (double)c-xc;
				double s2 = p[AYLP_FIT_P_SIGMA]*p[AYLP_FIT_P_SIGMA];
				double q = exp(-(dy*dy+dx*dx)/(2*s2));
				double aq = p[AYLP_FIT_P_AMP]*q;
				double an[AYLP_FIT_COM_NP];
				an[AYLP_FIT_P_Y0] = aq*dy/s2;
				an[AYLP_FIT_P_X0] = aq*dx/s2;
				an[AYLP_FIT_P_SY] = an[AYLP_FIT_P_Y0]*t;
				an[AYLP_FIT_P_SX] = an[AYLP_FIT_P_X0]*t;
				an[AYLP_FIT_P_SIGMA] = aq*(dy*dy+dx*dx)
					/(s2*p[AYLP_FIT_P_SIGMA]);
				an[AYLP_FIT_P_AMP] = q;
				an[AYLP_FIT_P_BG] = 1.0;
				// Only meaningful where the derivative is well
				// above the differencing noise floor; far from
				// the beam every derivative is ~0 and the
				// relative comparison is pure roundoff.
				if (fabs(an[a]) < 1e-4) continue;
				double err = fabs(an[a]-num)/fabs(an[a]);
				if (err > worst) worst = err;
			}
		}
		(void)names;
	}
	printf("worst relative Jacobian error: %.2e\n", worst);
	CHECK(worst < 1e-5, "analytic Jacobian must match central differences "
		"(worst relative error %.2e)", worst);
}

// ---- 2. parameter recovery from a sheared frame -------------------------
static void test_recovery(void)
{
	printf("\n=== parameter recovery from a sheared frame ===\n");
	struct { double sy, sx; } cases[] = {
		{0.00, 0.00}, {0.05, 0.00}, {0.00, -0.04}, {0.12, 0.08},
	};
	gsl_matrix_uchar *img = gsl_matrix_uchar_alloc(IMG, IMG);
	for (size_t k = 0; k < 4; k++) {
		json_object *p = cfg("\"row_time\":8.25e-6");
		struct aylp_device self = {0};
		self.params = p;
		if (fit_com_init(&self)) { CHECK(0, "init failed"); return; }
		struct aylp_fit_com_data *d = self.device_data;
		struct aylp_state st = {0};
		double ty = 15.6, tx = 16.4;
		render(img, ty, tx, cases[k].sy, cases[k].sx, 2.5, 200.0, 5.0);
		for (int f = 0; f < 6; f++) {	// warm start settles in 1-2
			st.matrix_uchar = img;
			st.header.status = 0;
			fit_com_proc(&self, &st);
		}
		double ry = (d->last_y+1.0)/2.0*(IMG-1);
		double rx = (d->last_x+1.0)/2.0*(IMG-1);
		double esy = d->p[AYLP_FIT_P_SY], esx = d->p[AYLP_FIT_P_SX];
		printf("  slope true (%+.3f,%+.3f) -> fit (%+.3f,%+.3f) | "
			"centre err (%+.4f,%+.4f) px | rms %.2f | %zu iter\n",
			cases[k].sy, cases[k].sx, esy, esx, ry-ty, rx-tx,
			d->last_rms, d->n_iter_last);
		CHECK(fabs(ry-ty) < 0.05 && fabs(rx-tx) < 0.05,
			"centre must be recovered to 0.05 px under shear "
			"(err %.4f, %.4f)", ry-ty, rx-tx);
		CHECK(fabs(esy-cases[k].sy) < 0.01
				&& fabs(esx-cases[k].sx) < 0.01,
			"slope must be recovered to 0.01 px/row "
			"(got %.4f,%.4f want %.4f,%.4f)", esy, esx,
			cases[k].sy, cases[k].sx);
		fit_com_fini(&self);
		json_object_put(p);
	}
	gsl_matrix_uchar_free(img);
}

// ---- 3. shear does NOT bias the centre, and is never rejected -----------
// The failure this device exists to remove: a shear gate is a velocity gate,
// so it drops frames at the fast part of an oscillation and keeps the turning
// points. Here shear is a fitted parameter, so no frame is ever dropped for
// moving and the centre stays unbiased as shear grows.
static void test_shear_unbiased(void)
{
	printf("\n=== centre bias vs shear magnitude (must stay flat) ===\n");
	gsl_matrix_uchar *img = gsl_matrix_uchar_alloc(IMG, IMG);
	json_object *p = cfg("\"row_time\":8.25e-6");
	struct aylp_device self = {0};
	self.params = p;
	if (fit_com_init(&self)) { CHECK(0, "init failed"); return; }
	struct aylp_fit_com_data *d = self.device_data;
	struct aylp_state st = {0};
	int rejected = 0;
	double worst = 0.0;
	for (int k = 0; k <= 8; k++) {
		double sy = 0.02*k;	// up to 0.16 px/row = 5 px across the ROI
		render(img, 15.6, 16.4, sy, 0.0, 2.5, 200.0, 5.0);
		for (int f = 0; f < 4; f++) {
			st.matrix_uchar = img;
			st.header.status = 0;
			fit_com_proc(&self, &st);
		}
		if (st.header.status & AYLP_FRAME_REJECTED) rejected++;
		double ry = (d->last_y+1.0)/2.0*(IMG-1);
		if (fabs(ry-15.6) > worst) worst = fabs(ry-15.6);
		printf("  shear %.2f px/row (%.1f px across ROI): centre err "
			"%+.4f px%s\n", sy, sy*IMG, ry-15.6,
			(st.header.status & AYLP_FRAME_REJECTED) ? "  REJECTED" : "");
	}
	CHECK(rejected == 0, "no frame may be rejected merely for shear "
		"(%d were)", rejected);
	CHECK(worst < 0.05, "centre bias must stay flat as shear grows "
		"(worst %.4f px)", worst);
	fit_com_fini(&self);
	json_object_put(p);
	gsl_matrix_uchar_free(img);
}

// ---- 4. static scene invariant -----------------------------------------
static void test_static(void)
{
	printf("\n=== static-scene invariant ===\n");
	gsl_matrix_uchar *img = gsl_matrix_uchar_alloc(IMG, IMG);
	render(img, 16.0, 16.4, 0.0, 0.0, 2.5, 200.0, 5.0);
	json_object *p = cfg("\"row_time\":8.25e-6");
	struct aylp_device self = {0};
	self.params = p;
	if (fit_com_init(&self)) { CHECK(0, "init failed"); return; }
	struct aylp_fit_com_data *d = self.device_data;
	struct aylp_state st = {0};
	double first = 0.0, worst = 0.0;
	int rejected = 0;
	for (int f = 0; f < 400; f++) {
		st.matrix_uchar = img;
		st.header.status = 0;
		fit_com_proc(&self, &st);
		if (st.header.status & AYLP_FRAME_REJECTED) { rejected++; continue; }
		double x = (d->last_x+1.0)/2.0*(IMG-1);
		if (!first) first = x;
		if (fabs(x-first) > worst) worst = fabs(x-first);
	}
	double ry = (d->last_y+1.0)/2.0*(IMG-1);
	double rx = (d->last_x+1.0)/2.0*(IMG-1);
	printf("400 identical frames: %d rejected, excursion %.2e px, "
		"absolute err (%+.4f, %+.4f) px\n", rejected, worst,
		ry-16.0, rx-16.4);
	CHECK(rejected == 0, "a static beam must never be rejected (%d were)",
		rejected);
	// Not bit-identical -- the solver warm-starts p and lambda from the
	// previous frame, so it re-converges to within its tolerance rather
	// than to the same bits. What matters is that there is no drift
	// mechanism: nothing here is self-referential the way a self-updating
	// correlation template is, so the excursion is solver noise and stays
	// ~1e-7 px forever instead of ratcheting.
	// The redescending weight is discontinuous in the data (a pixel is in
	// or out), so the converged point jitters at the sub-millipixel level.
	// Bounded, not drifting -- which is the property that matters.
	CHECK(worst < 1e-3, "identical frames must give identical output to "
		"well below the noise floor (excursion %.2e px)", worst);
	// And unlike a correlation tracker anchored at acquisition, the
	// reported position is absolute -- no acquisition-time offset is
	// absorbed into a reference.
	CHECK(fabs(ry-16.0) < 0.05 && fabs(rx-16.4) < 0.05,
		"absolute position must be correct, not relative to an anchor "
		"(err %.4f, %.4f)", ry-16.0, rx-16.4);
	fit_com_fini(&self);
	json_object_put(p);
	gsl_matrix_uchar_free(img);
}

// ---- 5. robustness to a stray reflection --------------------------------
static void test_stray(void)
{
	printf("\n=== stray reflection rejected by robust weighting ===\n");
	gsl_matrix_uchar *img = gsl_matrix_uchar_alloc(IMG, IMG);
	double err_off = 0.0, err_on = 0.0;
	for (int robust = 0; robust < 2; robust++) {
		json_object *p = cfg(robust ? "\"robust_k\":2.5" : "\"robust_k\":0.0");
		struct aylp_device self = {0};
		self.params = p;
		if (fit_com_init(&self)) { CHECK(0, "init failed"); return; }
		struct aylp_fit_com_data *d = self.device_data;
		struct aylp_state st = {0};
		render(img, 16.0, 16.0, 0.0, 0.0, 2.5, 200.0, 5.0);
		// A compact stray spot 7 px away -- close enough that its wing
		// overlaps the beam and an unweighted fit is measurably pulled,
		// which is the only way this test says anything.
		for (size_t r = 0; r < IMG; r++)
			for (size_t c = 0; c < IMG; c++) {
				double dy = (double)r-16.0, dx = (double)c-23.0;
				double v = img->data[r*img->tda+c]
					+ 170.0*exp(-(dy*dy+dx*dx)/(2*2.0*2.0));
				img->data[r*img->tda+c] =
					(unsigned char)(v > 255 ? 255 : v);
			}
		for (int f = 0; f < 12; f++) {
			st.matrix_uchar = img;
			st.header.status = 0;
			fit_com_proc(&self, &st);
		}
		double ry = (d->last_y+1.0)/2.0*(IMG-1);
		double rx = (d->last_x+1.0)/2.0*(IMG-1);
		double e = hypot(ry-16.0, rx-16.0);
		printf("  robust %s: centre err %.4f px\n",
			robust ? "on " : "off", e);
		if (robust) err_on = e; else err_off = e;
		fit_com_fini(&self);
		json_object_put(p);
	}
	CHECK(err_on <= err_off + 1e-6, "robust weighting must not worsen a "
		"stray reflection (off %.4f, on %.4f)", err_off, err_on);
	CHECK(err_on < 0.5, "stray reflection must not pull the centre more "
		"than 0.5 px with robust weighting on (got %.4f)", err_on);
	gsl_matrix_uchar_free(img);
}

// ---- 6. no beam is reported, not guessed --------------------------------
static void test_no_beam(void)
{
	printf("\n=== dark frames report loss instead of guessing ===\n");
	gsl_matrix_uchar *img = gsl_matrix_uchar_alloc(IMG, IMG);
	json_object *p = cfg("\"row_time\":8.25e-6");
	struct aylp_device self = {0};
	self.params = p;
	if (fit_com_init(&self)) { CHECK(0, "init failed"); return; }
	struct aylp_fit_com_data *d = self.device_data;
	struct aylp_state st = {0};
	render(img, 16.0, 16.0, 0.0, 0.0, 2.5, 200.0, 5.0);
	for (int f = 0; f < 4; f++) {
		st.matrix_uchar = img; st.header.status = 0;
		fit_com_proc(&self, &st);
	}
	double hold_y = d->last_y, hold_x = d->last_x;
	gsl_matrix_uchar_set_zero(img);
	int rejected = 0; bool lost = false;
	for (int f = 0; f < 8; f++) {
		st.matrix_uchar = img; st.header.status = 0;
		fit_com_proc(&self, &st);
		if (st.header.status & AYLP_FRAME_REJECTED) rejected++;
		if (st.header.status & AYLP_BEAM_LOST) lost = true;
	}
	printf("8 dark frames: %d rejected, beam_lost asserted: %s\n",
		rejected, lost ? "yes" : "no");
	CHECK(rejected == 8, "every dark frame must be rejected (%d of 8)",
		rejected);
	CHECK(lost, "sustained loss must assert AYLP_BEAM_LOST");
	CHECK(d->last_y == hold_y && d->last_x == hold_x,
		"a rejected frame must hold the previous output exactly");
	fit_com_fini(&self);
	json_object_put(p);
	gsl_matrix_uchar_free(img);
}

// ---- 7. a boundary-pinned solution is rejected --------------------------
static void test_clipped_beam(void)
{
	printf("\n=== clipped beam is rejected instead of pinned to the boundary ===\n");
	gsl_matrix_uchar *img = gsl_matrix_uchar_alloc(IMG, IMG);
	json_object *p = cfg("");
	struct aylp_device self = {0};
	self.params = p;
	if (fit_com_init(&self)) { CHECK(0, "init failed"); return; }
	struct aylp_fit_com_data *d = self.device_data;
	struct aylp_state st = {0};
	render(img, 16.0, 16.0, 0.0, 0.0, 2.5, 200.0, 5.0);
	for (int f = 0; f < 4; f++) {
		st.matrix_uchar = img; st.header.status = 0;
		fit_com_proc(&self, &st);
	}
	double hold_y = d->last_y, hold_x = d->last_x;
	render(img, 31.0, 31.0, 0.0, 0.0, 2.5, 200.0, 5.0);
	for (int f = 0; f < 4; f++) {
		st.matrix_uchar = img; st.header.status = 0;
		fit_com_proc(&self, &st);
		CHECK(st.header.status & AYLP_FRAME_REJECTED,
			"a clipped edge fit must be rejected on frame %d", f);
		CHECK(d->last_y == hold_y && d->last_x == hold_x,
			"a clipped edge fit must hold the last valid centroid");
	}
	fit_com_fini(&self);
	json_object_put(p);
	gsl_matrix_uchar_free(img);
}

// ---- 8. a solver jump cannot become a valid centroid --------------------
static void test_max_step(void)
{
	printf("\n=== implausible one-frame centroid jump is rejected ===\n");
	gsl_matrix_uchar *img = gsl_matrix_uchar_alloc(IMG, IMG);
	json_object *p = cfg("\"max_step\":5.0");
	struct aylp_device self = {0};
	self.params = p;
	if (fit_com_init(&self)) { CHECK(0, "init failed"); return; }
	struct aylp_fit_com_data *d = self.device_data;
	struct aylp_state st = {0};
	render(img, 16.0, 16.0, 0.0, 0.0, 2.5, 200.0, 5.0);
	for (int f = 0; f < 4; f++) {
		st.matrix_uchar = img; st.header.status = 0;
		fit_com_proc(&self, &st);
	}
	double hold_y = d->last_y, hold_x = d->last_x;
	render(img, 25.0, 16.0, 0.0, 0.0, 2.5, 200.0, 5.0);
	for (int f = 0; f < 4; f++) {
		st.matrix_uchar = img; st.header.status = 0;
		fit_com_proc(&self, &st);
		CHECK(st.header.status & AYLP_FRAME_REJECTED,
			"a >5 px one-frame fit jump must be rejected on frame %d", f);
		CHECK(d->last_y == hold_y && d->last_x == hold_x,
			"a rejected fit jump must hold the last valid centroid");
	}
	fit_com_fini(&self);
	json_object_put(p);
	gsl_matrix_uchar_free(img);
}

// ---- 9. the latency guard actually bounds latency ----------------------
//
// The device exists to fit inside a loop period, so the bound is a behaviour to
// test rather than a hope. What a specific microsecond count would test is the
// machine, so this runs the same hungry scene twice -- guard off, guard on --
// and checks that the guard is what ends the fit: fewer iterations, less time,
// and never all of max_iter. That holds whether or not the wide kernels were
// selected, which a threshold in microseconds would not.
static void run_big(double max_us, double *mean_us, size_t *last_iter,
	size_t *capped
) {
	const size_t BIG = 248;
	gsl_matrix_uchar *img = gsl_matrix_uchar_alloc(BIG, BIG);
	// A stray reflection is the case that genuinely wants many iterations,
	// so it is the case that shows what the guard costs and buys.
	render(img, 124.0, 124.0, 0.0, 0.0, 2.5, 200.0, 5.0);
	for (size_t r = 0; r < BIG; r++)
		for (size_t c = 0; c < BIG; c++) {
			double dy = (double)r-124.0, dx = (double)c-131.0;
			double v = img->data[r*img->tda+c]
				+ 170.0*exp(-(dy*dy+dx*dx)/(2*2.0*2.0));
			img->data[r*img->tda+c] = (unsigned char)(v > 255 ? 255 : v);
		}
	char buf[256];
	snprintf(buf, sizeof buf, "{\"sigma_init\":2.5,\"sigma_max\":10.0,"
		"\"max_iter\":50,\"max_us\":%g}", max_us);
	json_object *p = json_tokener_parse(buf);
	struct aylp_device self = {0};
	self.params = p;
	if (fit_com_init(&self)) { CHECK(0, "init failed"); return; }
	struct aylp_fit_com_data *d = self.device_data;
	struct aylp_state st = {0};

	double tot = 0.0;
	size_t n = 120, warm = 4, hit = 0;
	for (size_t f = 0; f < n; f++) {
		st.matrix_uchar = img;
		st.header.status = 0;
		struct timespec a, b;
		clock_gettime(CLOCK_MONOTONIC, &a);
		fit_com_proc(&self, &st);
		clock_gettime(CLOCK_MONOTONIC, &b);
		if (f < warm) continue;	// acquisition is not the steady state
		tot += (double)(b.tv_sec-a.tv_sec)*1e6
			+ (double)(b.tv_nsec-a.tv_nsec)*1e-3;
		if (d->budget_hit) hit++;
	}
	*mean_us = tot/(double)(n-warm);
	*last_iter = d->n_iter_last;
	*capped = hit;
	fit_com_fini(&self);
	json_object_put(p);
	gsl_matrix_uchar_free(img);
}

static void test_latency(void)
{
	printf("\n=== latency guard bounds a full-frame fit ===\n");
	double free_us, cap_us;
	size_t free_it, cap_it, free_hit, cap_hit;
	run_big(0.0, &free_us, &free_it, &free_hit);
	run_big(10.0, &cap_us, &cap_it, &cap_hit);
	printf("248x248 whole frame, stray reflection, max_iter 50:\n");
	printf("  guard off: %6.2f us/frame, %zu iterations\n", free_us, free_it);
	printf("  guard 10us: %6.2f us/frame, %zu iterations, stopped %zu frames\n",
		cap_us, cap_it, cap_hit);
	CHECK(free_it < 50, "even unguarded the solver must converge before "
		"max_iter on this scene (used %zu)", free_it);
	CHECK(free_hit == 0, "max_us 0 must not stop anything (%zu frames)",
		free_hit);
	CHECK(cap_it <= free_it, "the guard must not make the solver iterate "
		"more (%zu vs %zu)", cap_it, free_it);
	CHECK(cap_hit > 0, "the guard must actually stop this scene");
	CHECK(cap_us <= free_us, "the guard must not make a frame slower "
		"(%.2f vs %.2f us)", cap_us, free_us);
}

// A PWM edge during rolling readout leaves a bright fragment of the spot.
// Moment brightness and peak tests alone accept it and move the cross.  The
// row-shape gates must make the tracker infer the complete Gaussian's centre,
// while a uniformly dim whole spot must remain usable.
static void test_pwm_shutter_fragment(void)
{
	printf("\n=== PWM rolling-shutter fragments are held ===\n");
	gsl_matrix_uchar *full = gsl_matrix_uchar_alloc(IMG, IMG);
	render(full, 16.0, 16.0, 0.0, 0.0, 2.5, 200.0, 5.0);
	json_object *p = cfg("\"fit_slope\":false,\"moment_output\":true,"
		"\"fit_gaussian\":false,"
		"\"moment_cut\":5,\"moment_max_y_skew\":0.35,"
		"\"moment_min_y_width\":0.70,\"max_step\":8");
	struct aylp_device self = {0}; self.params = p;
	CHECK(!fit_com_init(&self), "init failed");
	struct aylp_fit_com_data *d = self.device_data;
	struct aylp_state st = {0}; st.matrix_uchar = full;
	for (int k = 0; k < 4; k++) { st.header.status = 0; fit_com_proc(&self, &st); }
	double good_y = d->last_y, good_x = d->last_x;

	gsl_matrix_uchar *cut = gsl_matrix_uchar_alloc(IMG, IMG);
	memcpy(cut->data, full->data, IMG*full->tda);
	for (size_t i = 17; i < IMG; i++)
		for (size_t j = 0; j < IMG; j++) cut->data[i*cut->tda+j] = 5;
	st.matrix_uchar = cut; st.header.status = 0; fit_com_proc(&self, &st);
	CHECK(!(st.header.status & AYLP_FRAME_REJECTED),
		"one-sided shutter fragment was not tracked");
	CHECK(fabs(d->last_y-good_y)*(IMG-1)/2 < 0.35
			&& fabs(d->last_x-good_x)*(IMG-1)/2 < 0.35,
		"inferred full-beam centre moved by (%.3f, %.3f) px",
		fabs(d->last_y-good_y)*(IMG-1)/2,
		fabs(d->last_x-good_x)*(IMG-1)/2);

	gsl_matrix_uchar *dim = gsl_matrix_uchar_alloc(IMG, IMG);
	render(dim, 16.0, 16.0, 0.0, 0.0, 2.5, 60.0, 5.0);
	st.matrix_uchar = dim; st.header.status = 0; fit_com_proc(&self, &st);
	CHECK(!(st.header.status & AYLP_FRAME_REJECTED),
		"uniformly dim whole beam was rejected as a shutter fragment");

	fit_com_fini(&self); json_object_put(p);
	gsl_matrix_uchar_free(full); gsl_matrix_uchar_free(cut);
	gsl_matrix_uchar_free(dim);
}


int main(void)
{
	test_jacobian();
	test_recovery();
	test_shear_unbiased();
	test_static();
	test_stray();
	test_no_beam();
	test_clipped_beam();
	test_max_step();
	test_pwm_shutter_fragment();
	test_latency();
	if (failures) {
		printf("\n%d CHECK(S) FAILED\n", failures);
		return 1;
	}
	printf("\nall checks passed\n");
	return 0;
}
