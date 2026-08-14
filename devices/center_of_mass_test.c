/** Unit tests for center_of_mass's partial-beam gate.
*
* The case these exist for is a rolling shutter reading out across a chopped
* source: some rows catch the source's on-window and some do not, so a frame can
* hold a full-brightness fragment of the beam. `min_peak` cannot see that -- the
* lit rows are at normal brightness -- and the centroid of the fragment is
* pulled toward those rows, which is a systematic error along the shutter axis
* rather than noise that averages away.
*
* `ref_cut` learns what a whole beam's row profile looks like and divides each
* frame by it. A whole beam gives the same ratio on every row, so the frame's
* brightness cancels; a cut beam gives two groups of rows, lit and not. The
* tests below pin both halves of that: cuts are caught (including on a beam too
* tight for any shape-based test), and frames that are merely dim are not.
*
* Run: ninja -C build test  (or build/test_center_of_mass) */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "anyloop.h"
#include "logging.h"
#include "center_of_mass.h"

static unsigned n_fail = 0;
static const char *current_test = "";

static int double_asc(const void *a, const void *b)
{
	double x=*(const double *)a, y=*(const double *)b;
	return x<y?-1:x>y;
}

#define CHECK(cond, ...) do { \
	if (!(cond)) { \
		n_fail++; \
		fprintf(stderr, "FAIL %s:%d [%s]: ", \
			__FILE__, __LINE__, current_test); \
		fprintf(stderr, __VA_ARGS__); \
		fprintf(stderr, "\n"); \
	} \
} while (0)

#define CHECK_NEAR(got, want, tol) \
	CHECK(fabs((got) - (want)) <= (tol), \
		#got " = %.6f, want %.6f (tol %g)", (double)(got), \
		(double)(want), (double)(tol))


// ---------------------------------------------------------------- fixtures

#define IMG 32

/** A Gaussian beam centred at (cy, cx), peaking at `peak` counts over a
* background of 2, with every row outside [cut_lo, cut_hi] forced to background
* -- which is what a rolling shutter does to a chopped source. */
static gsl_matrix_uchar *make_beam(
	double cy, double cx, double sigma, unsigned char peak,
	double cut_lo, double cut_hi
) {
	const unsigned char bg = 2;
	gsl_matrix_uchar *m = gsl_matrix_uchar_alloc(IMG, IMG);
	gsl_matrix_uchar_set_all(m, bg);
	for (size_t i = 0; i < IMG; i++) {
		if ((double)i < cut_lo || (double)i > cut_hi)
			continue;
		for (size_t j = 0; j < IMG; j++) {
			double dy = (double)i - cy, dx = (double)j - cx;
			double v = peak * exp(-(dy*dy + dx*dx)
				/ (2.0*sigma*sigma));
			double t = v + bg;
			m->data[i*m->tda + j] = t > 255.0 ? 255
				: (unsigned char)(t + 0.5);
		}
	}
	return m;
}

/** The whole beam, uncut, at the reference brightness. */
static gsl_matrix_uchar *full_beam(double cy, double cx)
{
	return make_beam(cy, cx, 2.0, 200, -1e9, 1e9);
}

static int make_com(struct aylp_device *dev, const char *params_json)
{
	memset(dev, 0, sizeof(*dev));
	dev->uri = "anyloop:center_of_mass";
	dev->params = json_tokener_parse(params_json);
	if (!dev->params) {
		fprintf(stderr, "test bug: bad json: %s\n", params_json);
		exit(2);
	}
	return center_of_mass_init(dev);
}

static void free_com(struct aylp_device *dev)
{
	if (dev->fini)
		dev->fini(dev);
	json_object_put(dev->params);
}

/** Push one frame through the device and return the emitted [y, x]. */
static void run_frame(
	struct aylp_device *dev, gsl_matrix_uchar *img, double *out_y,
	double *out_x
) {
	struct aylp_state state = {0};
	state.header.type = AYLP_T_MATRIX_UCHAR;
	state.matrix_uchar = img;
	int err = dev->proc(dev, &state);
	CHECK(!err, "proc returned %d", err);
	*out_y = state.vector->data[0];
	*out_x = state.vector->data[1];
}

/** Run the warmup so the reference is learned and the gate goes live. */
static void warm_up(struct aylp_device *dev, gsl_matrix_uchar *img, size_t n)
{
	double oy, ox;
	for (size_t k = 0; k < n; k++)
		run_frame(dev, img, &oy, &ox);
	struct aylp_center_of_mass_data *data = dev->device_data;
	CHECK(data->ref_ready, "reference not ready after %zu warmup frames", n);
}

/** Image-normalized coordinate of image row/column `p`. */
static double norm(double p)
{
	return -1.0 + 2.0*p/(IMG - 1);
}

#define WARMUP 5
// ref_cut 0.5 with the defaults. Measured on this fixture: a whole beam reads
// 1.00, a beam dimmed 4x still reads 0.61, and any cut through the significant
// rows reads 0.00 -- see test_dim_beam_is_not_a_cut
#define BASE "\"region_height\":12,\"region_width\":12,\"threshold\":10," \
	"\"track\":true,\"init_y\":16,\"init_x\":16," \
	"\"ref_warmup\":5,\"ref_cut\":0.5,"


// ------------------------------------------------------------------- tests

/** An uncut beam passes, and the warmup itself does not reject anything. */
static void test_full_beam_passes(void)
{
	current_test = "full_beam_passes";
	struct aylp_device dev;
	CHECK(!make_com(&dev, "{" BASE "}"), "init failed");
	gsl_matrix_uchar *img = full_beam(16.0, 16.0);
	warm_up(&dev, img, WARMUP);
	double oy, ox;
	run_frame(&dev, img, &oy, &ox);
	CHECK_NEAR(oy, norm(16.0), 5e-3);
	CHECK_NEAR(ox, norm(16.0), 5e-3);
	struct aylp_center_of_mass_data *data = dev.device_data;
	CHECK(data->n_held == 0, "uncut beam was rejected");
	gsl_matrix_uchar_free(img);
	free_com(&dev);
}

/** A beam that has moved still passes: the gate must track a real excursion,
* not freeze the loop the moment the beam is not where it was. */
static void test_moving_beam_passes(void)
{
	current_test = "moving_beam_passes";
	struct aylp_device dev;
	CHECK(!make_com(&dev, "{" BASE "}"), "init failed");
	gsl_matrix_uchar *img = full_beam(16.0, 16.0);
	warm_up(&dev, img, WARMUP);
	gsl_matrix_uchar_free(img);
	double oy, ox;
	for (double cy = 16.0; cy <= 19.0; cy += 0.5) {
		gsl_matrix_uchar *m = full_beam(cy, 16.0);
		run_frame(&dev, m, &oy, &ox);
		CHECK_NEAR(oy, norm(cy), 5e-3);
		gsl_matrix_uchar_free(m);
	}
	struct aylp_center_of_mass_data *data = dev.device_data;
	CHECK(data->n_held == 0, "a moving beam was rejected %zu times",
		data->n_held);
	free_com(&dev);
}

/** The defect: the cut frame passes every brightness test there is. Without the
* profile gate it is accepted and drags the centroid a pixel off a beam that
* never moved. */
static void test_brightness_tests_cannot_see_a_cut(void)
{
	current_test = "brightness_tests_cannot_see_a_cut";
	struct aylp_device dev;
	CHECK(!make_com(&dev, "{\"region_height\":12,\"region_width\":12,"
		"\"threshold\":10,\"track\":true,\"init_y\":16,\"init_x\":16,"
		"\"min_peak\":100}"), "init failed");
	double oy, ox;
	gsl_matrix_uchar *full = full_beam(16.0, 16.0);
	run_frame(&dev, full, &oy, &ox);
	CHECK_NEAR(oy, norm(16.0), 5e-3);

	gsl_matrix_uchar *cut = make_beam(16.0, 16.0, 2.0, 200, -1e9, 16.0);
	run_frame(&dev, cut, &oy, &ox);
	struct aylp_center_of_mass_data *data = dev.device_data;
	CHECK(data->n_held == 0,
		"min_peak was expected to let the cut frame through");
	CHECK(oy < norm(16.0) - 2.0/(IMG - 1)*0.7,
		"cut frame biased the centroid only to %.6f, from %.6f", oy,
		norm(16.0));

	gsl_matrix_uchar_free(full);
	gsl_matrix_uchar_free(cut);
	free_com(&dev);
}

/** ref_cut rejects that frame and holds the last good centroid. */
static void test_cut_is_rejected_and_held(void)
{
	current_test = "cut_is_rejected_and_held";
	struct aylp_device dev;
	CHECK(!make_com(&dev, "{" BASE "\"min_peak\":100}"), "init failed");
	gsl_matrix_uchar *full = full_beam(16.0, 16.0);
	warm_up(&dev, full, WARMUP);
	double oy, ox;
	run_frame(&dev, full, &oy, &ox);
	double good_y = oy, good_x = ox;

	gsl_matrix_uchar *cut = make_beam(16.0, 16.0, 2.0, 200, -1e9, 16.0);
	run_frame(&dev, cut, &oy, &ox);
	// held exactly, not merely close: nothing downstream should be able to
	// tell a held frame from a repeated measurement
	CHECK(oy == good_y && ox == good_x,
		"cut frame moved the output to (%.6f,%.6f), want held "
		"(%.6f,%.6f)", oy, ox, good_y, good_x);

	// the window did not chase the fragment either
	struct aylp_center_of_mass_data *data = dev.device_data;
	CHECK(data->win_y == 16, "window centre drifted to row %zu on a "
		"rejected frame", data->win_y);

	run_frame(&dev, full, &oy, &ox);
	CHECK_NEAR(oy, good_y, 1e-12);
	CHECK(data->n_held == 1, "n_held = %zu, want 1", data->n_held);
	CHECK(data->n_episodes == 1, "n_episodes = %zu, want 1",
		data->n_episodes);

	gsl_matrix_uchar_free(full);
	gsl_matrix_uchar_free(cut);
	free_com(&dev);
}

/** Cuts from either side, and a cut leaving a thin slice, are all rejected. The
* shutter phase decides which, and the bias flips sign with it, so a gate that
* caught only one would leave a residual that still beat against the chop. */
static void test_cut_from_either_side(void)
{
	current_test = "cut_from_either_side";
	struct aylp_device dev;
	CHECK(!make_com(&dev, "{" BASE "\"reacquire_after\":1000}"),
		"init failed");
	gsl_matrix_uchar *full = full_beam(16.0, 16.0);
	warm_up(&dev, full, WARMUP);
	double oy, ox;
	run_frame(&dev, full, &oy, &ox);
	double good_y = oy, good_x = ox;

	struct { const char *name; double lo, hi; } cuts[] = {
		{"bottom blanked", -1e9, 16.0},
		{"top blanked",    16.0,  1e9},
		{"thin slice",     15.0, 17.0},
		{"off-centre cut", -1e9, 17.0},
		{"one row short",   14.0, 1e9},
	};
	for (size_t k = 0; k < sizeof cuts / sizeof *cuts; k++) {
		gsl_matrix_uchar *cut = make_beam(16.0, 16.0, 2.0, 200,
			cuts[k].lo, cuts[k].hi);
		run_frame(&dev, cut, &oy, &ox);
		CHECK(oy == good_y && ox == good_x,
			"%s leaked into the output: (%.6f,%.6f)",
			cuts[k].name, oy, ox);
		gsl_matrix_uchar_free(cut);
		run_frame(&dev, full, &oy, &ox);
	}
	struct aylp_center_of_mass_data *data = dev.device_data;
	CHECK(data->n_held == 5, "n_held = %zu, want 5", data->n_held);

	gsl_matrix_uchar_free(full);
	free_com(&dev);
}

/** The property the whole design is for: a beam that is merely DIM is not a cut
* beam. Dividing by the reference makes the frame's own brightness cancel, so
* the source can drift a long way without the gate firing.
*
* Measured on this fixture, dimmest significant row as a fraction of the typical
* one: 1.00 at full brightness, 0.98 at 0.9x, 0.93 at 0.6x, 0.81 at 0.4x, 0.61
* at 0.25x. Every cut reads 0.00. The residual slope with brightness is the
* `threshold` subtraction, which is not linear in pixel value and so does not
* cancel perfectly -- it is the one part of this that is not exactly
* scale-free. */
static void test_dim_beam_is_not_a_cut(void)
{
	current_test = "dim_beam_is_not_a_cut";
	struct aylp_device dev;
	CHECK(!make_com(&dev, "{" BASE "}"), "init failed");
	gsl_matrix_uchar *bright = full_beam(16.0, 16.0);
	warm_up(&dev, bright, WARMUP);
	double oy, ox;

	// a quarter of the brightness the reference was learned at
	unsigned char levels[] = {180, 160, 120, 80, 50};
	for (size_t k = 0; k < sizeof levels / sizeof *levels; k++) {
		gsl_matrix_uchar *dim = make_beam(16.0, 16.0, 2.0, levels[k],
			-1e9, 1e9);
		run_frame(&dev, dim, &oy, &ox);
		CHECK_NEAR(oy, norm(16.0), 1.5e-2);
		gsl_matrix_uchar_free(dim);
	}
	struct aylp_center_of_mass_data *data = dev.device_data;
	CHECK(data->n_held == 0,
		"a whole but dimmed beam was rejected %zu times", data->n_held);

	gsl_matrix_uchar_free(bright);
	free_com(&dev);
}

/** A beam too tight for any shape-based test still works, because the reference
* is the beam's own profile and carries no assumption about how wide it is. A
* sigma-1 beam's own row-to-row change is 0.61 of its peak row, which no
* steepness threshold could tell from a cut. */
static void test_tight_beam(void)
{
	current_test = "tight_beam";
	struct aylp_device dev;
	// a tight beam has fewer rows above the floor, so lower it
	CHECK(!make_com(&dev, "{\"region_height\":12,\"region_width\":12,"
		"\"threshold\":10,\"track\":true,\"init_y\":16,\"init_x\":16,"
		"\"ref_warmup\":5,\"ref_cut\":0.5,\"ref_floor\":0.1}"),
		"init failed");
	gsl_matrix_uchar *full = make_beam(16.0, 16.0, 1.0, 200, -1e9, 1e9);
	warm_up(&dev, full, WARMUP);
	double oy, ox;
	run_frame(&dev, full, &oy, &ox);
	double good_y = oy, good_x = ox;
	struct aylp_center_of_mass_data *data = dev.device_data;
	CHECK(data->n_held == 0, "a whole tight beam was rejected");

	gsl_matrix_uchar *cut = make_beam(16.0, 16.0, 1.0, 200, -1e9, 16.0);
	run_frame(&dev, cut, &oy, &ox);
	CHECK(oy == good_y && ox == good_x,
		"a cut tight beam leaked into the output: (%.6f,%.6f)", oy, ox);
	CHECK(data->n_held == 1, "n_held = %zu, want 1", data->n_held);

	gsl_matrix_uchar_free(full);
	gsl_matrix_uchar_free(cut);
	free_com(&dev);
}

/** The bootstrap takes the row-wise maximum, so it learns the uncut profile
* even when the source is already chopping while it runs. A mean would learn a
* blend, set the reference low, and then wave the cut frames through. */
static void test_bootstrap_survives_chopping(void)
{
	current_test = "bootstrap_survives_chopping";
	struct aylp_device dev;
	CHECK(!make_com(&dev, "{\"region_height\":12,\"region_width\":12,"
		"\"threshold\":10,\"track\":true,\"init_y\":16,\"init_x\":16,"
		"\"ref_warmup\":8,\"ref_cut\":0.5,\"reacquire_after\":1000}"),
		"init failed");
	gsl_matrix_uchar *full = full_beam(16.0, 16.0);
	gsl_matrix_uchar *cut = make_beam(16.0, 16.0, 2.0, 200, -1e9, 16.0);
	double oy, ox;
	// warm up on an alternating chopped stream, which is what a cold start
	// into a live chop actually looks like
	for (int k = 0; k < 8; k++)
		run_frame(&dev, k % 2 ? cut : full, &oy, &ox);
	struct aylp_center_of_mass_data *data = dev.device_data;
	CHECK(data->ref_ready, "reference not ready after the warmup");

	// the reference must describe the WHOLE beam, so cuts are now caught
	run_frame(&dev, full, &oy, &ox);
	double good_y = oy;
	size_t held = data->n_held;
	run_frame(&dev, cut, &oy, &ox);
	CHECK(data->n_held == held + 1,
		"a cut frame was accepted after a chopped warmup");
	CHECK(oy == good_y, "cut frame leaked into the output: %.6f", oy);

	gsl_matrix_uchar_free(full);
	gsl_matrix_uchar_free(cut);
	free_com(&dev);
}

/** Until the reference exists there is nothing to compare against, so the gate
* stays off and frames go through on the brightness tests alone. Documented
* rather than incidental: at loop rate the warmup is milliseconds. */
static void test_gate_off_during_warmup(void)
{
	current_test = "gate_off_during_warmup";
	struct aylp_device dev;
	CHECK(!make_com(&dev, "{" BASE "}"), "init failed");
	double oy, ox;
	gsl_matrix_uchar *cut = make_beam(16.0, 16.0, 2.0, 200, -1e9, 16.0);
	run_frame(&dev, cut, &oy, &ox);
	struct aylp_center_of_mass_data *data = dev.device_data;
	CHECK(!data->ref_ready, "reference ready after one frame");
	CHECK(data->n_held == 0, "gate fired before the reference existed");
	gsl_matrix_uchar_free(cut);
	free_com(&dev);
}

/** ref_cut 0 (and absent) is exactly the old behaviour, and allocates nothing.
*/
static void test_gate_off_is_old_behaviour(void)
{
	current_test = "gate_off_is_old_behaviour";
	struct aylp_device dev;
	CHECK(!make_com(&dev, "{\"region_height\":12,\"region_width\":12,"
		"\"threshold\":10,\"track\":true,\"init_y\":16,\"init_x\":16}"),
		"init failed");
	struct aylp_center_of_mass_data *data = dev.device_data;
	CHECK(!data->ref, "reference allocated with the gate off");
	double oy, ox;
	gsl_matrix_uchar *cut = make_beam(16.0, 16.0, 2.0, 200, -1e9, 16.0);
	for (int k = 0; k < 10; k++)
		run_frame(&dev, cut, &oy, &ox);
	CHECK(data->n_held == 0, "cut frames rejected with the gate off");
	gsl_matrix_uchar_free(cut);
	free_com(&dev);
}

/** Losing the beam long enough to re-acquire throws the reference away: it
* describes a beam at the old window position, which is what we just stopped
* believing in. */
static void test_reference_reset_on_reacquire(void)
{
	current_test = "reference_reset_on_reacquire";
	struct aylp_device dev;
	CHECK(!make_com(&dev, "{" BASE "\"reacquire_after\":3}"),
		"init failed");
	gsl_matrix_uchar *full = full_beam(16.0, 16.0);
	warm_up(&dev, full, WARMUP);
	struct aylp_center_of_mass_data *data = dev.device_data;

	double oy, ox;
	// everything under threshold: no beam at all, for long enough to
	// trigger the re-acquire
	gsl_matrix_uchar *dark = make_beam(16.0, 16.0, 2.0, 0, -1e9, 1e9);
	for (int k = 0; k < 4; k++)
		run_frame(&dev, dark, &oy, &ox);
	CHECK(!data->ref_ready, "reference survived a re-acquire");
	CHECK(data->ref_seen == 0, "ref_seen = %zu after reset", data->ref_seen);

	// and it relearns. Not in WARMUP frames flat, though: re-acquiring off a
	// dark frame parks the window on the brightest speck of nothing, so the
	// window has to find the beam again (another reacquire_after frames)
	// before there is a profile to learn at all
	for (int k = 0; k < 20; k++)
		run_frame(&dev, full, &oy, &ox);
	CHECK(data->ref_ready, "reference did not relearn after re-acquiring");
	CHECK_NEAR(oy, norm(16.0), 5e-3);

	gsl_matrix_uchar_free(full);
	gsl_matrix_uchar_free(dark);
	free_com(&dev);
}

/** Nothing above threshold is still "no beam", and still holds. */
static void test_dark_frame_holds(void)
{
	current_test = "dark_frame_holds";
	struct aylp_device dev;
	CHECK(!make_com(&dev, "{" BASE "\"reacquire_after\":1000}"),
		"init failed");
	gsl_matrix_uchar *full = full_beam(16.0, 16.0);
	warm_up(&dev, full, WARMUP);
	double oy, ox;
	run_frame(&dev, full, &oy, &ox);
	double good_y = oy, good_x = ox;

	gsl_matrix_uchar *dark = make_beam(16.0, 16.0, 2.0, 0, -1e9, 1e9);
	run_frame(&dev, dark, &oy, &ox);
	CHECK(oy == good_y && ox == good_x, "dark frame did not hold");

	gsl_matrix_uchar_free(full);
	gsl_matrix_uchar_free(dark);
	free_com(&dev);
}

/** Auto settings must learn the background gate and seed/reacquire without
 * requiring ROI- or camera-specific constants. */
static void test_auto_tracking_settings(void)
{
	current_test="auto_tracking_settings";
	struct aylp_device dev;
	CHECK(!make_com(&dev,"{\"region_height\":12,\"region_width\":12,"
		"\"track\":true,\"threshold\":\"auto\","
		"\"min_peak\":\"auto\",\"init_y\":\"auto\","
		"\"init_x\":\"auto\",\"reacquire_after\":\"auto\"}"),
		"auto init failed");
	gsl_matrix_uchar *full=full_beam(20.0,21.0),*dark=
		make_beam(20.0,21.0,2.0,0,-1e9,1e9);
	double oy,ox;
	run_frame(&dev,full,&oy,&ox);
	struct aylp_center_of_mass_data *d=dev.device_data;
	CHECK(d->acquired,"auto init did not find the beam");
	CHECK(d->threshold<20,"auto threshold %u includes too much beam signal",
		d->threshold);
	CHECK(d->min_peak> d->threshold,"auto min_peak did not clear threshold");
	for(size_t k=0;k<30;k++)run_frame(&dev,dark,&oy,&ox);
	CHECK(d->reacquire_after==60,"auto reacquire_after %zu, want 60",
		d->reacquire_after);
	gsl_matrix_uchar_free(full);gsl_matrix_uchar_free(dark);free_com(&dev);
}

/** Config errors are caught at init rather than at 3 kHz. */
static void test_init_validation(void)
{
	struct aylp_device dev;

	current_test = "ref_cut_out_of_range";
	CHECK(make_com(&dev, "{" BASE "\"ref_cut\":1.5}"),
		"init accepted ref_cut > 1");
	free_com(&dev);
	CHECK(make_com(&dev, "{" BASE "\"ref_cut\":-0.1}"),
		"init accepted a negative ref_cut");
	free_com(&dev);

	current_test = "ref_knobs_out_of_range";
	CHECK(make_com(&dev, "{" BASE "\"ref_rate\":0}"),
		"init accepted ref_rate 0");
	free_com(&dev);
	CHECK(make_com(&dev, "{" BASE "\"ref_rate\":2}"),
		"init accepted ref_rate > 1");
	free_com(&dev);
	CHECK(make_com(&dev, "{" BASE "\"ref_floor\":1}"),
		"init accepted ref_floor 1");
	free_com(&dev);
	CHECK(make_com(&dev, "{" BASE "\"ref_warmup\":0}"),
		"init accepted ref_warmup 0");
	free_com(&dev);

	current_test = "ref_cut_needs_rows";
	CHECK(make_com(&dev, "{\"region_height\":3,\"region_width\":12,"
		"\"track\":true,\"ref_cut\":0.5}"),
		"init accepted a window too short to have a row profile");
	free_com(&dev);

	current_test = "ref_cut_ignored_without_track";
	CHECK(!make_com(&dev, "{\"region_height\":12,\"region_width\":12,"
		"\"ref_cut\":0.5}"), "init failed");
	struct aylp_center_of_mass_data *data = dev.device_data;
	CHECK(data->ref_cut == 0.0, "ref_cut = %G outside track mode, want 0",
		data->ref_cut);
	free_com(&dev);
}

/** The chop case end to end: the cut edge walks through the beam the way a chop
* beating against the frame rate does, while the source brightness also drifts.
* The emitted centroid must stay put throughout. */
static void test_chopped_sequence(void)
{
	current_test = "chopped_sequence";
	struct aylp_device dev;
	CHECK(!make_com(&dev, "{" BASE "\"min_peak\":50,"
		"\"reacquire_after\":1000}"), "init failed");
	gsl_matrix_uchar *full = full_beam(16.0, 16.0);
	warm_up(&dev, full, WARMUP);
	double oy, ox;
	run_frame(&dev, full, &oy, &ox);
	double good_y = oy, good_x = ox;
	const double tol = 0.1 * 2.0/(IMG - 1);	// a tenth of a pixel

	size_t n_cut = 0, n_caught = 0;
	unsigned char level = 200;
	for (double edge = 14.0; edge <= 18.0; edge += 0.5) {
		// the source dims as the sequence runs, to keep the gate
		// honest about telling "dim" from "cut"
		level = (unsigned char)(level * 0.95);
		gsl_matrix_uchar *cut = make_beam(16.0, 16.0, 2.0, level,
			-1e9, edge);
		struct aylp_center_of_mass_data *d = dev.device_data;
		size_t before = d->n_held;
		run_frame(&dev, cut, &oy, &ox);
		n_caught += d->n_held > before;
		n_cut++;
		CHECK(fabs(oy - good_y) <= tol && fabs(ox - good_x) <= tol,
			"cut at row %.1f (peak %u) moved the output to "
			"(%.6f,%.6f), more than 0.1 px off", edge, level,
			oy, ox);
		gsl_matrix_uchar_free(cut);

		gsl_matrix_uchar *whole = make_beam(16.0, 16.0, 2.0, level,
			-1e9, 1e9);
		run_frame(&dev, whole, &oy, &ox);
		CHECK(fabs(oy - good_y) <= tol,
			"whole beam at peak %u was mishandled: %.6f", level, oy);
		gsl_matrix_uchar_free(whole);
	}
	CHECK(n_caught == n_cut, "only %zu of %zu cuts were rejected",
		n_caught, n_cut);

	gsl_matrix_uchar_free(full);
	free_com(&dev);
}

/** Registration reports translation, not brightness-induced centroid motion.
* The deliberately asymmetric two-lobed target is not Gaussian and its lobes
* change relative brightness between frames. */
static gsl_matrix_uchar *make_pattern_size(
	size_t height, size_t width, double dy, double dx, double gain
)
{
	gsl_matrix_uchar *m=gsl_matrix_uchar_alloc(height,width);
	gsl_matrix_uchar_set_all(m,3);
	double cy=0.5*(height-1),cx=0.5*(width-1);
	for(size_t i=0;i<height;i++) for(size_t j=0;j<width;j++) {
		double y=i-dy,x=j-dx;
		double a=150*exp(-((y-(cy-3))*(y-(cy-3))
			+(x-(cx-2))*(x-(cx-2)))/8.0);
		double b=70*exp(-((y-(cy+4))*(y-(cy+4))
			+(x-(cx+3))*(x-(cx+3)))/18.0);
		// Lobes change in opposite directions: this is internal deformation,
		// not merely a global gain change, and should be downweighted.
		double v=3+gain*a+(2.0-gain)*b;
		m->data[i*m->tda+j]=(unsigned char)(v>255?255:v);
	}
	return m;
}

static gsl_matrix_uchar *make_pattern(double dy, double dx, double gain)
{
	gsl_matrix_uchar *m=gsl_matrix_uchar_alloc(IMG,IMG);
	gsl_matrix_uchar_set_all(m,3);
	for(size_t i=0;i<IMG;i++) for(size_t j=0;j<IMG;j++) {
		double y=i-dy,x=j-dx;
		double a=150*exp(-((y-13)*(y-13)+(x-14)*(x-14))/8.0);
		double b=70*exp(-((y-20)*(y-20)+(x-19)*(x-19))/18.0);
		double v=3+gain*a+(2.0-gain)*b;
		m->data[i*m->tda+j]=(unsigned char)(v>255?255:v);
	}
	return m;
}

static void test_registration_translation_and_brightness(void)
{
	current_test="registration_translation_and_brightness";
	struct aylp_device dev;
	CHECK(!make_com(&dev,"{\"region_height\":28,\"region_width\":28,"
		"\"track\":true,\"init_y\":16,\"init_x\":16,"
		"\"registration\":true}"),"init failed");
	double y0,x0,y1,x1,y2,x2;
	gsl_matrix_uchar *a=make_pattern(0,0,1.0);
	gsl_matrix_uchar *b=make_pattern(1.25,-0.75,0.90);
	gsl_matrix_uchar *c=make_pattern(1.25,-0.75,1.10);
	run_frame(&dev,a,&y0,&x0); run_frame(&dev,b,&y1,&x1);
	run_frame(&dev,c,&y2,&x2);
	double px=2.0/(IMG-1);
	CHECK_NEAR(y1-y0,1.25*px,0.16*px);
	CHECK_NEAR(x1-x0,-0.75*px,0.16*px);
	CHECK_NEAR(y2,y1,0.15*px);
	CHECK_NEAR(x2,x1,0.15*px);
	gsl_matrix_uchar_free(a); gsl_matrix_uchar_free(b);
	gsl_matrix_uchar_free(c); free_com(&dev);
}

/** Exercise the same registration contract at the two camera ROI sizes used
 * by the steering configurations.  Accuracy is in pixels, so changing the
 * normalization denominator must not change the answer. */
static void test_registration_roi_sizes(void)
{
	const size_t sizes[] = {112, 384};
	for (size_t k=0; k<sizeof sizes/sizeof sizes[0]; k++) {
		size_t n=sizes[k];
		char json[192];
		snprintf(json,sizeof json,"{\"region_height\":%zu,"
			"\"region_width\":%zu,\"track\":true,\"init_y\":%zu,"
			"\"init_x\":%zu,\"registration\":true}",n,n,n/2,n/2);
		struct aylp_device dev;
		CHECK(!make_com(&dev,json),"%zux%zu init failed",n,n);
		double y0,x0,y1,x1;
		gsl_matrix_uchar *a=make_pattern_size(n,n,0,0,1.0);
		gsl_matrix_uchar *b=make_pattern_size(n,n,1.25,-0.75,0.9);
		run_frame(&dev,a,&y0,&x0); run_frame(&dev,b,&y1,&x1);
		double px=2.0/(n-1);
		CHECK_NEAR(y1-y0,1.25*px,0.16*px);
		CHECK_NEAR(x1-x0,-0.75*px,0.16*px);
		struct aylp_center_of_mass_data *d=dev.device_data;
		CHECK(d->registration_sample_capacity==(n-4)*(n-4),
			"%zux%zu registration storage is not ROI-derived",n,n);
		gsl_matrix_uchar_free(a); gsl_matrix_uchar_free(b); free_com(&dev);
	}
}

static gsl_matrix_uchar *make_rolling_pattern(
	size_t n, double dy, double dx, double cut, bool keep_below
)
{
	gsl_matrix_uchar *m=make_pattern_size(n,n,dy,dx,1.0);
	if(!isfinite(cut))return m;
	for(size_t y=0;y<n;y++) {
		bool dark=keep_below ? (double)y<cut : (double)y>cut;
		if(dark)memset(m->data+y*m->tda,3,n);
	}
	return m;
}

/** A rolling shutter may expose only one side of an otherwise valid pattern.
 * The flux centroid then follows the exposure edge instead of the target.
 * Registration is allowed to accept the common rigid structure or hold its
 * last result, but it must never make the known position estimate worse.
 * Run at both production ROI sizes so neither the target nor the acceptance
 * criterion encodes a preferred frame geometry. */
static void test_registration_rolling_shutter_only_improves(void)
{
	current_test="registration_rolling_shutter_only_improves";
	const size_t sizes[]={112,384};
	for(size_t z=0;z<sizeof sizes/sizeof sizes[0];z++) {
		size_t n=sizes[z];
		char plain_json[160],reg_json[192];
		snprintf(plain_json,sizeof plain_json,"{\"region_height\":%zu,"
			"\"region_width\":%zu,\"track\":true,\"init_y\":%zu,"
			"\"init_x\":%zu}",n,n,n/2,n/2);
		snprintf(reg_json,sizeof reg_json,"{\"region_height\":%zu,"
			"\"region_width\":%zu,\"track\":true,\"init_y\":%zu,"
			"\"init_x\":%zu,\"ref_cut\":0.5,\"ref_warmup\":2,"
			"\"registration\":true}",n,n,n/2,n/2);
		struct aylp_device plain,reg;
		CHECK(!make_com(&plain,plain_json),"%zux%zu plain init failed",n,n);
		CHECK(!make_com(&reg,reg_json),"%zux%zu registration init failed",n,n);
		const double dy[]={0,.5,.5,.5,.5,1.0};
		const double dx[]={0,-.3,-.3,-.3,-.3,-.6};
		const double cut[]={NAN,NAN,.5*n,.5*n,.5*n,NAN};
		const bool below[]={false,false,true,false,true,false};
		double ref_y=0,ref_x=0,sum_plain=0,sum_reg=0;
		for(size_t k=0;k<sizeof dy/sizeof dy[0];k++) {
			gsl_matrix_uchar *img=make_rolling_pattern(
				n,dy[k],dx[k],cut[k],below[k]);
			double py,px,ry,rx;
			run_frame(&plain,img,&py,&px);run_frame(&reg,img,&ry,&rx);
			py=(py+1)*(n-1)/2;px=(px+1)*(n-1)/2;
			ry=(ry+1)*(n-1)/2;rx=(rx+1)*(n-1)/2;
			if(!k){ref_y=ry;ref_x=rx;}
			double want_y=ref_y+dy[k],want_x=ref_x+dx[k];
			double ep=hypot(py-want_y,px-want_x);
			double er=hypot(ry-want_y,rx-want_x);
			CHECK(er<=ep+0.10,
				"%zux%zu frame %zu: registration error %.3f px worsened COM %.3f px",
				n,n,k,er,ep);
			sum_plain+=ep*ep;sum_reg+=er*er;
			gsl_matrix_uchar_free(img);
		}
		CHECK(sum_reg<=0.25*sum_plain,
			"%zux%zu rolling registration RMS did not materially improve COM "
			"(squared errors %.4f vs %.4f)",n,n,sum_reg,sum_plain);
		free_com(&plain);free_com(&reg);
	}
}

/** Long rolling-keyframe chains must not integrate subpixel bias.  The target
 * stays at a known position while its two lobes exchange brightness through
 * many keyframe renewals; this is the failure mode that looked like tens of
 * pixels of motion in the repository captures. */
static void test_registration_long_term_shape_drift(void)
{
	current_test="registration_long_term_shape_drift";
	const size_t sizes[]={112,384};
	for(size_t z=0;z<sizeof sizes/sizeof sizes[0];z++) {
		size_t n=sizes[z];char json[192];
		snprintf(json,sizeof json,"{\"region_height\":%zu,"
			"\"region_width\":%zu,\"track\":true,\"init_y\":%zu,"
			"\"init_x\":%zu,\"registration\":true}",n,n,n/2,n/2);
		struct aylp_device dev;
		CHECK(!make_com(&dev,json),"%zux%zu long-term init failed",n,n);
		double y0=0,x0=0,se=0,last=0;
		for(size_t k=0;k<256;k++) {
			double gain=1.0+0.3*sin(2*M_PI*k/37.0);
			gsl_matrix_uchar *img=make_pattern_size(n,n,0,0,gain);
			double y,x;run_frame(&dev,img,&y,&x);
			y=(y+1)*(n-1)/2;x=(x+1)*(n-1)/2;
			if(!k){y0=y;x0=x;}
			double e=hypot(y-y0,x-x0);se+=e*e;last=e;
			gsl_matrix_uchar_free(img);
		}
		CHECK(sqrt(se/256)<0.30,
			"%zux%zu changing-shape RMS drift %.3f px",n,n,sqrt(se/256));
		CHECK(last<0.30,"%zux%zu terminal drift %.3f px",n,n,last);
		free_com(&dev);
	}
}

static double test_bilinear_bg(
	const gsl_matrix_uchar *src, double y, double x, double bg
)
{
	if(y<0||x<0||y+1>=src->size1||x+1>=src->size2)return bg;
	size_t y0=(size_t)y,x0=(size_t)x;
	double fy=y-y0,fx=x-x0;
	const unsigned char *p=src->data+y0*src->tda+x0;
	return (1-fy)*((1-fx)*p[0]+fx*p[1])
		+fy*((1-fx)*p[src->tda]+fx*p[src->tda+1]);
}

/** Inject a closed, known subpixel trajectory into one real illuminated frame.
 * Real sensor noise, beam shape, reflections and saturation are retained; only
 * translation, global gain and background are controlled. */
static void test_registration_injected_real_frame(
	const gsl_matrix_uchar *base, const char *path
)
{
	size_t h=base->size1,w=base->size2;
	char json[256];
	snprintf(json,sizeof json,"{\"region_height\":%zu,"
		"\"region_width\":%zu,\"track\":true,\"init_y\":%zu,"
		"\"init_x\":%zu,\"registration\":true}",h,w,h/2,w/2);
	struct aylp_device dev;
	CHECK(!make_com(&dev,json),"injected tracker init failed for %s",path);
	gsl_matrix_uchar *img=gsl_matrix_uchar_alloc(h,w);
	double border=0,bn=0;
	for(size_t y=0;y<h;y++)for(size_t x=0;x<w;x++)
		if(!y||!x||y+1==h||x+1==w){border+=base->data[y*base->tda+x];bn++;}
	double bg=bn?border/bn:0,oy0=0,ox0=0,se_y=0,se_x=0;
	double syy=0,sxx=0,smy=0,smx=0,max_e=0,last_e=0;
	const size_t frames=257;
	for(size_t k=0;k<frames;k++) {
		double phase=2*M_PI*k/(frames-1);
		double dy=4.0*sin(phase),dx=3.0*sin(2*phase);
		double gain=0.85+0.05*cos(3*phase);
		double off=2.0+1.0*sin(5*phase);
		for(size_t y=0;y<h;y++)for(size_t x=0;x<w;x++) {
			double v=test_bilinear_bg(base,y-dy,x-dx,bg);
			v=off+bg+gain*(v-bg);
			img->data[y*img->tda+x]=(unsigned char)(v<0?0:v>255?255:v+0.5);
		}
		double oy,ox;run_frame(&dev,img,&oy,&ox);
		oy=(oy+1)*(h-1)/2;ox=(ox+1)*(w-1)/2;
		if(!k){oy0=oy;ox0=ox;}
		double my=oy-oy0,mx=ox-ox0,ey=my-dy,ex=mx-dx;
		se_y+=ey*ey;se_x+=ex*ex;syy+=dy*dy;sxx+=dx*dx;
		smy+=dy*my;smx+=dx*mx;
		last_e=hypot(ey,ex);if(last_e>max_e)max_e=last_e;
	}
	double rms_y=sqrt(se_y/frames),rms_x=sqrt(se_x/frames);
	double gain_y=syy?smy/syy:0,gain_x=sxx?smx/sxx:0;
	fprintf(stderr,"injected replay: %s: %zux%zu, RMS %.4f/%.4f px, "
		"gain %.4f/%.4f, max %.4f px, terminal %.4f px\n",
		path,h,w,rms_y,rms_x,gain_y,gain_x,max_e,last_e);
	CHECK(rms_y<0.35&&rms_x<0.35,
		"injected RMS %.3f/%.3f px on %s",rms_y,rms_x,path);
	CHECK(fabs(gain_y-1)<0.05&&fabs(gain_x-1)<0.05,
		"injected gain %.3f/%.3f on %s",gain_y,gain_x,path);
	CHECK(last_e<0.35,"injected closed-loop terminal error %.3f px on %s",
		last_e,path);
	gsl_matrix_uchar_free(img);free_com(&dev);
}

static bool test_pair_delta(const gsl_matrix_uchar *a,const gsl_matrix_uchar *b,
	double *dy,double *dx)
{
	char json[192];size_t h=a->size1,w=a->size2;
	snprintf(json,sizeof json,"{\"region_height\":%zu,\"region_width\":%zu,"
		"\"track\":true,\"init_y\":%zu,\"init_x\":%zu,"
		"\"registration\":true}",h,w,h/2,w/2);
	struct aylp_device dev;if(make_com(&dev,json))return false;
	double ay,ax,by,bx;run_frame(&dev,(gsl_matrix_uchar*)a,&ay,&ax);
	size_t held=((struct aylp_center_of_mass_data*)dev.device_data)->n_held;
	run_frame(&dev,(gsl_matrix_uchar*)b,&by,&bx);
	bool ok=((struct aylp_center_of_mass_data*)dev.device_data)->n_held==held;
	*dy=(by-ay)*(h-1)/2;*dx=(bx-ax)*(w-1)/2;free_com(&dev);return ok;
}

/** Ground-truth-free checks on three actual consecutive camera frames. A
 * valid translation estimate must invert and compose even when brightness and
 * beam shape evolve between exposures. */
static void test_registration_cycle(const gsl_matrix_uchar *a,
	const gsl_matrix_uchar *b,const gsl_matrix_uchar *c,const char *path)
{
	double aby,abx,bay,bax,bcy,bcx,cby,cbx,acy,acx;
	bool ok=test_pair_delta(a,b,&aby,&abx)&&test_pair_delta(b,a,&bay,&bax)
		&&test_pair_delta(b,c,&bcy,&bcx)&&test_pair_delta(c,b,&cby,&cbx)
		&&test_pair_delta(a,c,&acy,&acx);
	CHECK(ok,"cycle registration rejected a selected bright pair on %s",path);
	if(!ok)return;
	double inv_ab=hypot(aby+bay,abx+bax);
	double inv_bc=hypot(bcy+cby,bcx+cbx);
	double compose=hypot(aby+bcy-acy,abx+bcx-acx);
	CHECK(inv_ab<0.25&&inv_bc<0.25,
		"forward/backward closure %.4f/%.4f px on %s",inv_ab,inv_bc,path);
	CHECK(compose<0.75,"three-frame composition error %.4f px on %s",compose,path);
	fprintf(stderr,"cycle replay: %s: inverse %.4f/%.4f px, composition %.4f px\n",
		path,inv_ab,inv_bc,compose);
}

/** Optional real-capture replay used during tracker validation.  Keeping it in
 * this binary makes the 40-byte AYLP parsing and the production device path
 * reproducible without making the repository's large captures test fixtures. */
static void replay_registration_capture(const char *path)
{
	current_test="registration_capture_replay";
	FILE *fp=fopen(path,"rb");
	CHECK(fp,"cannot open registration capture %s",path);
	if(!fp)return;
	unsigned char header[40];
	CHECK(fread(header,1,sizeof header,fp)==sizeof header,
		"short AYLP header in %s",path);
	uint64_t h=0,w=0;
	memcpy(&h,header+8,sizeof h); memcpy(&w,header+16,sizeof w);
	CHECK(header[6]==AYLP_T_MATRIX_UCHAR && h>=7 && w>=7,
		"%s is not a usable MATRIX_UCHAR capture",path);
	if(header[6]!=AYLP_T_MATRIX_UCHAR || h<7 || w<7){fclose(fp);return;}
	char json[256];
	snprintf(json,sizeof json,"{\"region_height\":%llu,"
		"\"region_width\":%llu,\"track\":true,\"init_y\":%llu,"
		"\"init_x\":%llu,\"threshold\":1,\"min_peak\":100,"
		"\"reacquire_after\":600,\"registration\":true}",
		(unsigned long long)h,(unsigned long long)w,
		(unsigned long long)(h/2),(unsigned long long)(w/2));
	struct aylp_device dev;
	CHECK(!make_com(&dev,json),"capture tracker init failed for %s",path);
	char plain_json[256];
	snprintf(plain_json,sizeof plain_json,"{\"region_height\":%llu,"
		"\"region_width\":%llu,\"track\":true,\"init_y\":%llu,"
		"\"init_x\":%llu,\"threshold\":1,\"min_peak\":100,"
		"\"reacquire_after\":600}",
		(unsigned long long)h,(unsigned long long)w,
		(unsigned long long)(h/2),(unsigned long long)(w/2));
	struct aylp_device plain;
	CHECK(!make_com(&plain,plain_json),"plain capture tracker init failed for %s",path);
	gsl_matrix_uchar *img=gsl_matrix_uchar_alloc(h,w);
	gsl_matrix_uchar *injected_base=(h<=384&&w<=384)
		?gsl_matrix_uchar_alloc(h,w):0;
	gsl_matrix_uchar *cycle_prev2=injected_base?gsl_matrix_uchar_alloc(h,w):0;
	gsl_matrix_uchar *cycle_prev1=injected_base?gsl_matrix_uchar_alloc(h,w):0;
	gsl_matrix_uchar *cycle_best[3]={0};
	for(size_t k=0;k<3&&injected_base;k++)cycle_best[k]=gsl_matrix_uchar_alloc(h,w);
	double prev_score2=0,prev_score1=0,best_cycle_score=-1;
	bool prev_bright2=false,prev_bright1=false,have_cycle=false;
	bool have_injected_base=false;
	double injected_score=-1;
	double py=0,px=0,max_step=0,oy=0,ox=0;
	double mean_y=0,mean_x=0,m2_y=0,m2_x=0,last_ay=0,last_ax=0;
	double diff_y=0,diff_x=0;
	double delta_y=0,delta_x=0,delta_m2_y=0,delta_m2_x=0;
	double first_delta_y=0,first_delta_x=0,last_delta_y=0,last_delta_x=0;
	size_t frames=0,bright=0,bright_accepted=0,recoveries=0;
	size_t accepted_n=0,diff_n=0;
	double bright_us[2000]; size_t bright_us_n=0;
	double plain_mean_y=0,plain_mean_x=0,plain_m2_y=0,plain_m2_x=0;
	double plain_diff_y=0,plain_diff_x=0,plain_last_y=0,plain_last_x=0;
	size_t plain_n=0,plain_diff_n=0,plain_bright_accepted=0;
	bool plain_previous=false;
	bool was_bright=false,previous_accepted=false;
	do {
		if(fread(img->data,1,h*w,fp)!=h*w)break;
		unsigned char peak=0;double frame_score=0;
		for(size_t p=0;p<h*w;p++) {
			if(img->data[p]>peak)peak=img->data[p];
			frame_score+=img->data[p];
		}
		bool is_bright=peak>=100;
		if(frames>=2&&is_bright&&prev_bright1&&prev_bright2) {
			double cycle_score=fmin(frame_score,fmin(prev_score1,prev_score2));
			if(cycle_score>best_cycle_score) {
				gsl_matrix_uchar_memcpy(cycle_best[0],cycle_prev2);
				gsl_matrix_uchar_memcpy(cycle_best[1],cycle_prev1);
				gsl_matrix_uchar_memcpy(cycle_best[2],img);
				best_cycle_score=cycle_score;have_cycle=true;
			}
		}
		if(is_bright&&injected_base&&frame_score>injected_score) {
			gsl_matrix_uchar_memcpy(injected_base,img);
			have_injected_base=true;
			injected_score=frame_score;
		}
		if(is_bright&&!was_bright)recoveries++;
		was_bright=is_bright;
		size_t held_before=((struct aylp_center_of_mass_data*)dev.device_data)->n_held;
		size_t plain_held_before=((struct aylp_center_of_mass_data*)plain.device_data)->n_held;
		double plain_oy,plain_ox;
		run_frame(&plain,img,&plain_oy,&plain_ox);
		bool plain_accepted=((struct aylp_center_of_mass_data*)plain.device_data)->n_held
			==plain_held_before;
		if(is_bright&&plain_accepted) {
			double pay=(plain_oy+1)*(h-1)/2,pax=(plain_ox+1)*(w-1)/2;
			plain_bright_accepted++;plain_n++;
			double qy=pay-plain_mean_y,qx=pax-plain_mean_x;
			plain_mean_y+=qy/plain_n;plain_mean_x+=qx/plain_n;
			plain_m2_y+=qy*(pay-plain_mean_y);plain_m2_x+=qx*(pax-plain_mean_x);
			if(plain_previous){plain_diff_y+=(pay-plain_last_y)*(pay-plain_last_y);
				plain_diff_x+=(pax-plain_last_x)*(pax-plain_last_x);plain_diff_n++;}
			plain_last_y=pay;plain_last_x=pax;plain_previous=true;
		} else plain_previous=false;
		struct timespec proc_t0,proc_t1;
		clock_gettime(CLOCK_MONOTONIC_RAW,&proc_t0);
		run_frame(&dev,img,&oy,&ox);
		clock_gettime(CLOCK_MONOTONIC_RAW,&proc_t1);
		if(is_bright && bright_us_n<sizeof bright_us/sizeof bright_us[0])
			bright_us[bright_us_n++]=1e6*((proc_t1.tv_sec-proc_t0.tv_sec)
				+1e-9*(proc_t1.tv_nsec-proc_t0.tv_nsec));
		bool accepted=((struct aylp_center_of_mass_data*)dev.device_data)->n_held==held_before;
		if(is_bright){bright++;if(accepted)bright_accepted++;}
		if(accepted) {
			double ay=(oy+1)*(h-1)/2,ax=(ox+1)*(w-1)/2;
			// Independent whole-frame, border-subtracted intensity centroid.
			// It is shape-sensitive, but registration drifting tens of pixels
			// relative to it is a useful diagnostic even when neither is truth.
			double border=0,bn=0,rw=0,ry=0,rx=0;
			for(size_t y=0;y<h;y++)for(size_t x=0;x<w;x++)
				if(!y||!x||y+1==h||x+1==w){border+=img->data[y*img->tda+x];bn++;}
			double bg=bn?border/bn:0;
			for(size_t y=0;y<h;y++)for(size_t x=0;x<w;x++){
				double v=img->data[y*img->tda+x]-bg;
				if(v<=2)continue;
				rw+=v;ry+=v*y;rx+=v*x;
			}
			double ddy=rw?ay-ry/rw:0,ddx=rw?ax-rx/rw:0;
			accepted_n++;
			if(accepted_n==1){first_delta_y=ddy;first_delta_x=ddx;}
			last_delta_y=ddy;last_delta_x=ddx;
			double qdy=ddy-delta_y,qdx=ddx-delta_x;
			delta_y+=qdy/accepted_n;delta_x+=qdx/accepted_n;
			delta_m2_y+=qdy*(ddy-delta_y);delta_m2_x+=qdx*(ddx-delta_x);
			double qy=ay-mean_y,qx=ax-mean_x;
			mean_y+=qy/accepted_n;mean_x+=qx/accepted_n;
			m2_y+=qy*(ay-mean_y);m2_x+=qx*(ax-mean_x);
			if(previous_accepted){
				diff_y+=(ay-last_ay)*(ay-last_ay);
				diff_x+=(ax-last_ax)*(ax-last_ax);diff_n++;
			}
			last_ay=ay;last_ax=ax;
		}
		if(frames&&accepted&&previous_accepted)max_step=fmax(max_step,hypot((oy-py)*(h-1)/2,
			(ox-px)*(w-1)/2));
		py=oy;px=ox;frames++;
		previous_accepted=accepted;
		if(cycle_prev1) {
			gsl_matrix_uchar_memcpy(cycle_prev2,cycle_prev1);
			gsl_matrix_uchar_memcpy(cycle_prev1,img);
			prev_score2=prev_score1;prev_score1=frame_score;
			prev_bright2=prev_bright1;prev_bright1=is_bright;
		}
	} while(frames<2000 && fread(header,1,sizeof header,fp)==sizeof header);
	struct aylp_center_of_mass_data *d=dev.device_data;
	CHECK(frames>=2,"%s has fewer than two complete frames",path);
	if(bright) {
		CHECK(d->registration_ready,"registration never became ready on %s",path);
		// Peak alone cannot distinguish a complete exposure from a bright rolling
		// fragment.  Repository surveys independently classify as many as 24% of
		// peak-bright frames as unusable, so require registration to retain the
		// clear majority without forcing it to accept those malformed frames.
		CHECK(bright_accepted*10>=bright*7,
			"registration accepted only %zu of %zu illuminated frames on %s",
			bright_accepted,bright,path);
	} else {
		CHECK(d->n_held==frames,
			"beamless capture emitted %zu updates instead of holding all frames: %s",
			frames-d->n_held,path);
	}
	qsort(bright_us,bright_us_n,sizeof *bright_us,double_asc);
	fprintf(stderr,"registration replay: %s: %llux%llu, %zu frames, "
		"%zu samples, cond %.1f, max step %.4f px, bright %zu/%zu accepted, "
		"cost p50 %.1f us / p99 %.1f us / max %.1f us, "
		"dark transitions %zu, held %zu, RMS %.4f/%.4f px, "
		"HF %.4f/%.4f px, raw disagreement %.4f/%.4f px, "
		"raw drift %.4f/%.4f px\n",path,
		(unsigned long long)h,(unsigned long long)w,frames,
		d->registration_n_samples,d->registration_condition,max_step,bright_accepted,bright,
		bright_us_n?bright_us[bright_us_n/2]:0,
		bright_us_n?bright_us[(bright_us_n-1)*99/100]:0,
		bright_us_n?bright_us[bright_us_n-1]:0,
		recoveries,d->n_held,
		accepted_n?sqrt(m2_y/accepted_n):0,
		accepted_n?sqrt(m2_x/accepted_n):0,
		diff_n?sqrt(diff_y/(2*diff_n)):0,
		diff_n?sqrt(diff_x/(2*diff_n)):0,
		accepted_n?sqrt(delta_m2_y/accepted_n):0,
		accepted_n?sqrt(delta_m2_x/accepted_n):0,
		last_delta_y-first_delta_y,last_delta_x-first_delta_x);
	fprintf(stderr,"plain replay: %s: bright %zu/%zu accepted, RMS %.4f/%.4f px, "
		"HF %.4f/%.4f px\n",path,plain_bright_accepted,bright,
		plain_n?sqrt(plain_m2_y/plain_n):0,plain_n?sqrt(plain_m2_x/plain_n):0,
		plain_diff_n?sqrt(plain_diff_y/(2*plain_diff_n)):0,
		plain_diff_n?sqrt(plain_diff_x/(2*plain_diff_n)):0);
	if(have_injected_base)test_registration_injected_real_frame(injected_base,path);
	if(have_cycle)test_registration_cycle(cycle_best[0],cycle_best[1],cycle_best[2],path);
	for(size_t k=0;k<3;k++)if(cycle_best[k])gsl_matrix_uchar_free(cycle_best[k]);
	if(cycle_prev2)gsl_matrix_uchar_free(cycle_prev2);
	if(cycle_prev1)gsl_matrix_uchar_free(cycle_prev1);
	if(injected_base)gsl_matrix_uchar_free(injected_base);
	gsl_matrix_uchar_free(img); free_com(&plain); free_com(&dev); fclose(fp);
}

struct movement_stats {
	size_t records,valid;
	double sum[2],sum2[2],first[2],second[2];
	size_t first_n,second_n;
};

/** Validate the only recorded actuator-step COM pair without pretending the
 * beamless 0 V file contains a position reference. The test pins flag handling
 * and checks that the driven acquisition is finite, continuous and stationary
 * enough to be a usable physical measurement. */
static bool read_movement_com(const char *path,struct movement_stats *s)
{
	FILE *fp=fopen(path,"rb");
	CHECK(fp,"cannot open movement COM %s",path);if(!fp)return false;
	struct aylp_header h;double v[2];
	while(fread(&h,1,sizeof h,fp)==sizeof h) {
		if(h.magic!=AYLP_MAGIC||h.type!=AYLP_T_VECTOR
				||h.log_dim.y!=2||h.log_dim.x!=1
				||fread(v,sizeof *v,2,fp)!=2) {
			CHECK(false,"invalid movement COM record in %s",path);
			fclose(fp);return false;
		}
		s->records++;
		if(h.status&AYLP_NO_SIGNAL)continue;
		CHECK(isfinite(v[0])&&isfinite(v[1]),"non-finite COM in %s",path);
		if(!isfinite(v[0])||!isfinite(v[1]))continue;
		for(size_t a=0;a<2;a++){s->sum[a]+=v[a];s->sum2[a]+=v[a]*v[a];}
		s->valid++;
		// These recordings contain 500 samples; splitting by record number keeps
		// the drift check independent of missing/flagged values.
		if(s->records<=250){s->first[0]+=v[0];s->first[1]+=v[1];s->first_n++;}
		else{s->second[0]+=v[0];s->second[1]+=v[1];s->second_n++;}
	}
	CHECK(feof(fp),"truncated movement COM %s",path);fclose(fp);return true;
}

static void test_recorded_movement_com(const char *zero_path,const char *step_path)
{
	current_test="recorded_movement_com";
	struct movement_stats z={0},d={0};
	if(!read_movement_com(zero_path,&z)||!read_movement_com(step_path,&d))return;
	CHECK(z.records==500&&z.valid==0,
		"0 V recording should contain 500 explicitly beamless samples, got %zu/%zu valid",
		z.valid,z.records);
	CHECK(d.records==500&&d.valid==500,
		"+0.5 V recording should contain 500 valid samples, got %zu/%zu",d.valid,d.records);
	if(!d.valid||!d.first_n||!d.second_n)return;
	double mean[2],sd[2],drift[2];
	for(size_t a=0;a<2;a++) {
		mean[a]=d.sum[a]/d.valid;
		sd[a]=sqrt(fmax(0,d.sum2[a]/d.valid-mean[a]*mean[a]));
		drift[a]=d.second[a]/d.second_n-d.first[a]/d.first_n;
		CHECK(sd[a]<0.06,"driven COM axis %zu is not stable: SD %.4f",a,sd[a]);
		CHECK(fabs(drift[a])<0.08,"driven COM axis %zu drifts %.4f",a,drift[a]);
	}
	CHECK(hypot(mean[0],mean[1])>0.20,
		"driven COM has no resolved offset: mean %.4f/%.4f",mean[0],mean[1]);
	fprintf(stderr,"movement replay: 0 V %zu/%zu valid; +0.5 V mean %.4f/%.4f, "
		"SD %.4f/%.4f, half-drift %.4f/%.4f\n",z.valid,z.records,
		mean[0],mean[1],sd[0],sd[1],drift[0],drift[1]);
}


int main(int argc, char **argv)
{
	// the gate tests deliberately trip the beam-lost path, and its log
	// lines are expected output rather than failures
	log_init(LOG_ERROR);

	test_full_beam_passes();
	test_moving_beam_passes();
	test_brightness_tests_cannot_see_a_cut();
	test_cut_is_rejected_and_held();
	test_cut_from_either_side();
	test_dim_beam_is_not_a_cut();
	test_tight_beam();
	test_bootstrap_survives_chopping();
	test_gate_off_during_warmup();
	test_gate_off_is_old_behaviour();
	test_reference_reset_on_reacquire();
	test_dark_frame_holds();
	test_auto_tracking_settings();
	test_init_validation();
	test_chopped_sequence();
	test_registration_translation_and_brightness();
	test_registration_roi_sizes();
	test_registration_rolling_shutter_only_improves();
	test_registration_long_term_shape_drift();
	for(int i=1;i<argc;i++) {
		if(!strcmp(argv[i],"--movement-check")&&i+2<argc) {
			test_recorded_movement_com(argv[i+1],argv[i+2]);i+=2;
		} else replay_registration_capture(argv[i]);
	}

	if (n_fail)
		fprintf(stderr, "%u check(s) failed\n", n_fail);
	else
		fprintf(stderr, "all center_of_mass gate checks passed\n");
	return n_fail ? 1 : 0;
}
