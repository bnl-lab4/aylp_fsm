/** Unit tests for gain_test, against a simulated plant of known gain.
*
* The point of these is that the number the device prints is the number the
* plant actually has. A staircase gain measurement has three ways to lie --
* averaging before the plant has settled, fitting levels where the sensor has
* run out of range, and mistaking beam drift for gain -- and each is checked
* here against a plant whose true gain is known exactly.
*
* The simulated plant is a one-pole lag with a pure delay, driven by whatever
* command the device publishes, which is how the device is used in a real
* pipeline (sensor -> gain_test -> DAC). Wall-clock dwells are scaled right
* down; the device times its levels off CLOCK_MONOTONIC, so the test sleeps a
* realistic loop period per iteration rather than spinning.
*
* Run: ninja -C build test  (or build/test_gain_test) */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "anyloop.h"
#include "logging.h"
#include "gain_test.h"

static unsigned n_fail = 0;
static const char *current_test = "";

#define CHECK(cond, ...) do { \
	if (!(cond)) { \
		n_fail++; \
		fprintf(stderr, "FAIL %s:%d [%s]: ", \
			__FILE__, __LINE__, current_test); \
		fprintf(stderr, __VA_ARGS__); \
		fprintf(stderr, "\n"); \
	} \
} while (0)

// a loop period in the region the rig actually runs at (~3.8 kHz)
#define LOOP_NS 260000

struct plant {
	double gain;		// response units per command unit
	double pole;		// one-pole lag coefficient, 0 = instant
	double drift;		// response units added per iteration
	double clip;		// |response| saturates here (0 = no limit)
	double noise;		// uniform noise half-width
	double y;		// state
	unsigned long n;
	unsigned seed;
};

static double plant_step(struct plant *p, double u)
{
	double target = p->gain * u + p->drift * (double)p->n++;
	p->y = p->pole * p->y + (1.0 - p->pole) * target;
	double out = p->y;
	if (p->clip > 0.0) {
		// a sensor that has run out of range: it stops moving AND stops
		// showing noise, exactly as a tracked center_of_mass holding its
		// last coordinate does
		if (out > p->clip) return p->clip;
		if (out < -p->clip) return -p->clip;
	}
	if (p->noise > 0.0) {
		p->seed = p->seed * 1103515245u + 12345u;
		out += p->noise * (2.0 * ((p->seed >> 16) / 65536.0) - 1.0);
	}
	return out;
}

/** Run the device against the plant until it raises AYLP_DONE, or until
 * max_iter iterations have gone by. Returns iterations run. */
static size_t run(struct aylp_device *dev, struct plant *p, size_t max_iter)
{
	struct timespec sl = { 0, LOOP_NS };
	double u = 0.0;
	size_t i = 0;
	for (; i < max_iter; i++) {
		gsl_vector *in = gsl_vector_alloc(2);
		gsl_vector_set(in, 0, plant_step(p, u));
		gsl_vector_set(in, 1, 0.0);
		struct aylp_state state = {0};
		state.header.type = AYLP_T_VECTOR;
		state.header.units = AYLP_U_MINMAX;
		state.vector = in;
		int err = dev->proc(dev, &state);
		CHECK(!err, "proc returned %d at iteration %zu", err, i);
		u = state.vector->data[0];
		gsl_vector_free(in);
		if (state.header.status & AYLP_DONE) { i++; break; }
		nanosleep(&sl, NULL);
	}
	return i;
}

static int make(struct aylp_device *dev, const char *params_json)
{
	memset(dev, 0, sizeof(*dev));
	dev->uri = "anyloop:gain_test";
	dev->params = json_tokener_parse(params_json);
	if (!dev->params) {
		fprintf(stderr, "test bug: bad json: %s\n", params_json);
		exit(2);
	}
	return gain_test_init(dev);
}

static void unmake(struct aylp_device *dev)
{
	if (dev->fini) dev->fini(dev);
	json_object_put(dev->params);
}

// 21 levels from -0.12 to +0.12, the shipped sweep, with the dwells scaled
// down so the whole thing runs in about a second
#define SWEEP21 "{\"index_cmd\":0,\"index_err\":0,\"out_size\":2," \
	"\"low\":-0.12,\"high\":0.12,\"step\":0.012," \
	"\"dwell\":0.04,\"settle_frac\":0.5,\"warmup\":0.05,\"ramp\":0.02," \
	"\"pixel_scale\":31.5,\"output_file\":\"\""


/** The headline claim: the fitted slope is the plant's gain. */
static void test_recovers_gain(void)
{
	current_test = "recovers_gain";
	struct aylp_device dev;
	CHECK(!make(&dev, SWEEP21 "}"), "init failed");
	struct plant p = { .gain = 4.05, .pole = 0.7, .noise = 0.002,
		.seed = 1 };
	run(&dev, &p, 200000);
	struct aylp_gain_test_data *d = dev.device_data;
	CHECK(d->stage == GT_DONE, "did not finish (stage %d)", d->stage);
	CHECK(d->n_done == 21, "measured %zu levels, want 21", d->n_done);
	CHECK(d->have_fit, "no fit was produced");
	CHECK(fabs(d->slope - 4.05) < 0.02, "slope %g, want 4.05", d->slope);
	CHECK(d->r2 > 0.999, "R2 %g on a linear plant", d->r2);
	// with a linear plant the largest residual is noise, not curvature
	CHECK(d->nonlin < 0.01, "max residual %g on a linear plant", d->nonlin);
	unmake(&dev);
}

/** The command really is stepped over the whole requested range, in the
 * requested increments -- the fit is only as good as the staircase under it. */
static void test_staircase_shape(void)
{
	current_test = "staircase_shape";
	struct aylp_device dev;
	CHECK(!make(&dev, SWEEP21 "}"), "init failed");
	struct plant p = { .gain = 4.05, .pole = 0.7, .seed = 1 };
	run(&dev, &p, 200000);
	struct aylp_gain_test_data *d = dev.device_data;
	CHECK(fabs(d->cmds[0] + 0.12) < 1e-9, "first level %g, want -0.12",
		d->cmds[0]);
	CHECK(fabs(d->cmds[d->n_done-1] - 0.12) < 1e-9,
		"last level %g, want +0.12", d->cmds[d->n_done-1]);
	for (size_t i = 1; i < d->n_done; i++)
		if (fabs(d->cmds[i] - d->cmds[i-1] - 0.012) > 1e-9) {
			CHECK(0, "level %zu steps by %g, want 0.012", i,
				d->cmds[i] - d->cmds[i-1]);
			break;
		}
	// and each level is averaged over a real number of settled samples
	for (size_t i = 0; i < d->n_done; i++)
		if (d->counts[i] < 8) {
			CHECK(0, "level %zu averaged only %zu samples", i,
				d->counts[i]);
			break;
		}
	unmake(&dev);
}

/** A plant that compresses at large command is reported as compressing: the
 * full-range slope is pulled down, the small-signal fit is not, and the
 * nonlinearity figure is well clear of the noise. This is the case the rig
 * actually has -- the FSM loses a few percent of gain above ~50 mV. */
static void test_compression_is_visible(void)
{
	current_test = "compression_is_visible";
	struct aylp_device dev;
	CHECK(!make(&dev, SWEEP21 ",\"fit_range\":0.05}"), "init failed");
	// y = K*(u - c*u^3): 4.05 at zero, ~7% down at the ends
	struct plant p = { .gain = 4.05, .pole = 0.7, .noise = 0.0005,
		.seed = 3 };
	struct timespec sl = { 0, LOOP_NS };
	double u = 0.0;
	for (size_t i = 0; i < 200000; i++) {
		double eff = u - 5.0 * u*u*u;
		gsl_vector *in = gsl_vector_alloc(2);
		gsl_vector_set(in, 0, plant_step(&p, eff));
		gsl_vector_set(in, 1, 0.0);
		struct aylp_state state = {0};
		state.header.type = AYLP_T_VECTOR;
		state.header.units = AYLP_U_MINMAX;
		state.vector = in;
		dev.proc(&dev, &state);
		u = state.vector->data[0];
		gsl_vector_free(in);
		if (state.header.status & AYLP_DONE) break;
		nanosleep(&sl, NULL);
	}
	struct aylp_gain_test_data *d = dev.device_data;
	CHECK(d->have_fit, "no fit was produced");
	// full-range slope sits below the small-signal one
	CHECK(d->slope < 4.0, "full-range slope %g, want < 4.0 on a "
		"compressing plant", d->slope);
	CHECK(d->n_small > 0, "small-signal fit did not run");
	CHECK(d->slope_small > d->slope + 0.05,
		"small-signal slope %g is not above the full-range %g",
		d->slope_small, d->slope);
	CHECK(fabs(d->slope_small - 4.05) < 0.12,
		"small-signal slope %g, want about 4.05", d->slope_small);
	CHECK(d->nonlin > 0.005, "curvature %g was not reported", d->nonlin);
	unmake(&dev);
}

/** Levels where the sensor has run out of range are excluded, not fitted, and
 * two of them in a row end the sweep instead of walking the beam further off.
 * Without this a saturating sensor quietly flattens the slope. */
static void test_excludes_saturated_levels(void)
{
	current_test = "excludes_saturated_levels";
	struct aylp_device dev;
	// resp_max 0.5 with a gain of 4.05 puts the limit at +-0.123 command
	// units, so the sweep should die a couple of levels from the top
	CHECK(!make(&dev, SWEEP21 ",\"resp_max\":0.35}"), "init failed");
	struct plant p = { .gain = 4.05, .pole = 0.7, .noise = 0.002,
		.clip = 0.4, .seed = 5 };
	run(&dev, &p, 200000);
	struct aylp_gain_test_data *d = dev.device_data;
	CHECK(d->stage == GT_DONE, "did not finish (stage %d)", d->stage);
	CHECK(d->n_done < 21, "sweep ran to the end past a saturated sensor "
		"(%zu levels)", d->n_done);
	size_t bad = 0;
	for (size_t i = 0; i < d->n_done; i++) if (!d->valid[i]) bad++;
	CHECK(bad >= 2, "only %zu levels marked unusable", bad);
	CHECK(d->have_fit, "no fit from the usable levels");
	CHECK(fabs(d->slope - 4.05) < 0.05,
		"slope %g from the usable levels, want 4.05", d->slope);
	unmake(&dev);
}

/** A beam that drifts during the sweep biases the plain slope. An up-and-down
 * sweep repeats every command at a later time, which is what makes gain and
 * drift separable at all; the drift-carrying fit then recovers the true gain,
 * and the gap between the two slopes is what says to distrust the plain one. */
static void test_drift_is_separable(void)
{
	current_test = "drift_is_separable";
	struct aylp_device dev;
	CHECK(!make(&dev, SWEEP21 ",\"updown\":true}"), "init failed");
	// drift of ~1e-5 units per iteration is ~0.04 units/s at this loop rate
	struct plant p = { .gain = 4.05, .pole = 0.7, .noise = 0.001,
		.drift = 1e-5, .seed = 7 };
	run(&dev, &p, 400000);
	struct aylp_gain_test_data *d = dev.device_data;
	CHECK(d->have_fit, "no fit was produced");
	CHECK(fabs(d->slope_drift - 4.05) < 0.05,
		"drift-corrected slope %g, want 4.05", d->slope_drift);
	CHECK(d->drift_rate > 0.0, "drift rate %g, want positive",
		d->drift_rate);
	// the drift shows up as an offset between the two branches, which is
	// exactly what a real hysteresis measurement has to be read against
	CHECK(d->hysteresis > 0.01,
		"drift of %g units/s left no up/down offset (%g)",
		d->drift_rate, d->hysteresis);
	unmake(&dev);
}

/** A single rising staircase cannot separate drift from gain -- command and
 * time are the same axis -- so the device must decline to report a drift
 * number rather than inventing a split. */
static void test_monotonic_sweep_declines_drift_fit(void)
{
	current_test = "monotonic_sweep_declines_drift_fit";
	struct aylp_device dev;
	CHECK(!make(&dev, SWEEP21 "}"), "init failed");
	struct plant p = { .gain = 4.05, .pole = 0.7, .noise = 0.001,
		.drift = 1e-5, .seed = 7 };
	run(&dev, &p, 200000);
	struct aylp_gain_test_data *d = dev.device_data;
	CHECK(d->have_fit, "no fit was produced");
	CHECK(d->slope > 4.10, "test bug: drift did not bias the plain slope "
		"(%g)", d->slope);
	CHECK(d->slope_drift == 0.0 && d->drift_rate == 0.0,
		"reported a drift split (%g, %g) from a monotonic sweep",
		d->slope_drift, d->drift_rate);
	unmake(&dev);
}

/** An up-and-down sweep measures both branches and reports the offset between
 * them, which on a real actuator is hysteresis. */
static void test_updown_measures_hysteresis(void)
{
	current_test = "updown_measures_hysteresis";
	struct aylp_device dev;
	CHECK(!make(&dev, SWEEP21 ",\"updown\":true}"), "init failed");
	struct plant p = { .gain = 4.05, .pole = 0.7, .noise = 0.001,
		.seed = 11 };
	run(&dev, &p, 400000);
	struct aylp_gain_test_data *d = dev.device_data;
	CHECK(d->n_done == 42, "measured %zu levels, want 42", d->n_done);
	CHECK(d->branch[0] == GT_BRANCH_UP, "first level is not the up branch");
	CHECK(d->branch[d->n_done-1] == GT_BRANCH_DOWN,
		"last level is not the down branch");
	CHECK(fabs(d->cmds[d->n_done-1] + 0.12) < 1e-9,
		"the down branch ended at %g, want -0.12", d->cmds[d->n_done-1]);
	// this plant has no hysteresis, so the branches must agree
	CHECK(fabs(d->hysteresis) < 0.01,
		"hysteresis %g reported on a memoryless plant", d->hysteresis);
	CHECK(fabs(d->slope - 4.05) < 0.02, "slope %g, want 4.05", d->slope);
	unmake(&dev);
}

/** Configuration errors are refused at init rather than producing a run that
 * cannot mean anything. */
static void test_rejects_bad_config(void)
{
	current_test = "rejects_bad_config";
	struct aylp_device dev;
	CHECK(make(&dev, "{\"low\":0.1,\"high\":-0.1,\"step\":0.01}"),
		"accepted low > high");
	unmake(&dev);
	CHECK(make(&dev, "{\"low\":-0.1,\"high\":0.1,\"step\":0}"),
		"accepted step 0");
	unmake(&dev);
	CHECK(make(&dev, "{\"low\":-0.1,\"high\":0.1,\"step\":0.2}"),
		"accepted a two-level sweep, which is a step and not a line");
	unmake(&dev);
	CHECK(make(&dev, "{\"index_cmd\":5,\"out_size\":2}"),
		"accepted index_cmd outside the output vector");
	unmake(&dev);
}


int main(void)
{
	// the saturation test deliberately trips warnings; they are expected
	// output rather than failures
	log_init(LOG_ERROR);

	test_recovers_gain();
	test_staircase_shape();
	test_compression_is_visible();
	test_excludes_saturated_levels();
	test_drift_is_separable();
	test_monotonic_sweep_declines_drift_fit();
	test_updown_measures_hysteresis();
	test_rejects_bad_config();

	if (n_fail) fprintf(stderr, "%u check(s) failed\n", n_fail);
	else fprintf(stderr, "all gain_test checks passed\n");
	return n_fail ? 1 : 0;
}
