/** Unit tests for pid's park-on-lost-signal behaviour.
*
* Holding the last command is right for a brief dropout and wrong for a beam
* that is actually gone: the error is frozen, so the integrator winds on it at a
* constant rate and a lost beam becomes a large excursion instead of a stalled
* one. Measured on this rig, that was 0.44 command units over 8.5 s.
*
* `park_after` counts consecutive iterations carrying AYLP_NO_SIGNAL and, past
* the threshold, drives every output element to `park_value` and empties the
* controller's memory so recovery starts from rest.
*
* Run: ninja -C build test  (or build/test_pid_park) */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "anyloop.h"
#include "logging.h"
#include "pid.h"

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


// ---------------------------------------------------------------- fixtures

static int make_pid(struct aylp_device *dev, const char *params_json)
{
	memset(dev, 0, sizeof(*dev));
	dev->uri = "anyloop:pid";
	dev->params = json_tokener_parse(params_json);
	if (!dev->params) {
		fprintf(stderr, "test bug: bad json: %s\n", params_json);
		exit(2);
	}
	return pid_init(dev);
}

static void free_pid(struct aylp_device *dev)
{
	if (dev->fini)
		dev->fini(dev);
	json_object_put(dev->params);
}

/** Run one iteration with error `e` on both elements and the given status;
* returns the two output elements. The input vector is rebuilt every call
* because pid replaces state->vector with its own result. */
static void step(
	struct aylp_device *dev, double e, aylp_status status,
	double *o0, double *o1
) {
	gsl_vector *in = gsl_vector_alloc(2);
	gsl_vector_set_all(in, e);
	struct aylp_state state = {0};
	state.header.type = AYLP_T_VECTOR;
	state.header.units = AYLP_U_MINMAX;
	state.header.status = status;
	state.vector = in;
	int err = dev->proc(dev, &state);
	CHECK(!err, "proc returned %d", err);
	// the pipeline type must survive a park, or every device after pid
	// fails on the first parked iteration
	CHECK(state.header.type != AYLP_T_UNCHANGED,
		"pid nulled the pipeline type");
	*o0 = state.vector->data[0 * state.vector->stride];
	*o1 = state.vector->data[1 * state.vector->stride];
	gsl_vector_free(in);
}

#define PARK30 "{\"type\":\"vector\",\"p\":5,\"i\":5,\"clamp\":10," \
	"\"park_after\":30,\"park_value\":0}"


// ------------------------------------------------------------------- tests

/** With a signal present, pid controls normally and park never fires. */
static void test_controls_with_signal(void)
{
	current_test = "controls_with_signal";
	struct aylp_device dev;
	CHECK(!make_pid(&dev, PARK30), "init failed");
	double a, b;
	for (int k = 0; k < 100; k++)
		step(&dev, 0.1, 0, &a, &b);
	// p=5 on an error of 0.1 is -0.5 before the integrator adds anything
	CHECK(a < -0.4, "output %g does not look like active control", a);
	struct aylp_pid_data *d = dev.device_data;
	CHECK(!d->parked, "parked while the signal was present");
	CHECK(d->no_signal == 0, "no_signal = %zu with a signal present",
		d->no_signal);
	free_pid(&dev);
}

/** Short dropouts do NOT park: the whole point of the threshold is that an
* ordinary 1-3 frame dropout is ridden out by holding, which is what the loop
* wants. Only a real loss parks. */
static void test_short_dropout_does_not_park(void)
{
	current_test = "short_dropout_does_not_park";
	struct aylp_device dev;
	CHECK(!make_pid(&dev, PARK30), "init failed");
	double a, b;
	for (int k = 0; k < 50; k++)
		step(&dev, 0.1, 0, &a, &b);
	double before = a;
	// 29 consecutive no-signal iterations, one short of the threshold
	for (int k = 0; k < 29; k++)
		step(&dev, 0.1, AYLP_NO_SIGNAL, &a, &b);
	struct aylp_pid_data *d = dev.device_data;
	CHECK(!d->parked, "parked after 29 iterations with park_after 30");
	CHECK(a != 0.0, "output was zeroed before the threshold");
	CHECK(fabs(a - before) < 0.5, "output ran away during a short dropout: "
		"%g -> %g", before, a);
	free_pid(&dev);
}

/** The 30th consecutive no-signal iteration parks both elements at exactly
* park_value -- not near it. */
static void test_parks_on_threshold(void)
{
	current_test = "parks_on_threshold";
	struct aylp_device dev;
	CHECK(!make_pid(&dev, PARK30), "init failed");
	double a, b;
	for (int k = 0; k < 50; k++)
		step(&dev, 0.1, 0, &a, &b);
	CHECK(a != 0.0, "test bug: output was already zero before the park");

	for (int k = 0; k < 29; k++)
		step(&dev, 0.1, AYLP_NO_SIGNAL, &a, &b);
	CHECK(a != 0.0, "parked one iteration early");
	step(&dev, 0.1, AYLP_NO_SIGNAL, &a, &b);	// the 30th
	CHECK(a == 0.0 && b == 0.0,
		"30th no-signal iteration gave [%g, %g], want exactly [0, 0]",
		a, b);
	struct aylp_pid_data *d = dev.device_data;
	CHECK(d->parked, "parked flag not set");

	// and it stays parked for as long as the signal is missing
	for (int k = 0; k < 500; k++) {
		step(&dev, 0.1, AYLP_NO_SIGNAL, &a, &b);
		if (a != 0.0 || b != 0.0) {
			CHECK(0, "output left the park at iteration %d: "
				"[%g, %g]", k, a, b);
			break;
		}
	}
	free_pid(&dev);
}

/** Parking empties the controller's memory, so when the beam comes back the
* first good frame is not handed a stale integral. This is the half that makes
* parking safe rather than merely quiet. */
static void test_recovery_starts_from_rest(void)
{
	current_test = "recovery_starts_from_rest";
	struct aylp_device dev;
	CHECK(!make_pid(&dev, PARK30), "init failed");
	double a, b;
	// wind up a substantial integral on a real error
	for (int k = 0; k < 2000; k++)
		step(&dev, 0.2, 0, &a, &b);
	double wound = a;
	CHECK(fabs(wound) > 1.0, "test bug: integral did not wind up (%g)",
		wound);

	for (int k = 0; k < 30; k++)
		step(&dev, 0.2, AYLP_NO_SIGNAL, &a, &b);
	CHECK(a == 0.0, "did not park");

	// signal returns, and the error is now zero. If the accumulator had
	// survived, the output would jump straight back to the wound-up value
	step(&dev, 0.0, 0, &a, &b);
	struct aylp_pid_data *d = dev.device_data;
	CHECK(!d->parked, "still parked after the signal returned");
	CHECK(fabs(a) < 0.01,
		"output was %g on the first good frame; a stale integral of "
		"about %g came back with it", a, wound);
	free_pid(&dev);
}

/** A dropout that recovers before the threshold resets the counter, so dropouts
* do not accumulate toward a park across unrelated episodes. */
static void test_counter_resets(void)
{
	current_test = "counter_resets";
	struct aylp_device dev;
	CHECK(!make_pid(&dev, PARK30), "init failed");
	double a, b;
	struct aylp_pid_data *d = dev.device_data;
	for (int episode = 0; episode < 10; episode++) {
		for (int k = 0; k < 20; k++)
			step(&dev, 0.1, AYLP_NO_SIGNAL, &a, &b);
		CHECK(d->no_signal == 20, "no_signal = %zu, want 20",
			d->no_signal);
		step(&dev, 0.1, 0, &a, &b);
		CHECK(d->no_signal == 0, "counter not reset by a good frame");
		CHECK(!d->parked, "parked after 10 episodes of 20 (< 30 each)");
	}
	free_pid(&dev);
}

/** park_after 0 disables the feature entirely -- the behaviour this device has
* always had, so existing configs are untouched. */
static void test_disabled_by_default(void)
{
	current_test = "disabled_by_default";
	struct aylp_device dev;
	CHECK(!make_pid(&dev, "{\"type\":\"vector\",\"p\":5,\"i\":5,"
		"\"clamp\":10}"), "init failed");
	double a, b;
	for (int k = 0; k < 50; k++)
		step(&dev, 0.1, 0, &a, &b);
	for (int k = 0; k < 500; k++)
		step(&dev, 0.1, AYLP_NO_SIGNAL, &a, &b);
	struct aylp_pid_data *d = dev.device_data;
	CHECK(d->park_after == 0, "park_after defaulted to %zu, want 0",
		d->park_after);
	CHECK(!d->parked, "parked with the feature disabled");
	CHECK(a != 0.0, "output zeroed with the feature disabled");
	free_pid(&dev);
}

/** park_value is honoured, so a rig whose neutral is not zero command can say
* so. */
static void test_nonzero_park_value(void)
{
	current_test = "nonzero_park_value";
	struct aylp_device dev;
	CHECK(!make_pid(&dev, "{\"type\":\"vector\",\"p\":5,\"i\":5,"
		"\"clamp\":10,\"park_after\":5,\"park_value\":-1.25}"),
		"init failed");
	double a, b;
	for (int k = 0; k < 20; k++)
		step(&dev, 0.1, 0, &a, &b);
	for (int k = 0; k < 5; k++)
		step(&dev, 0.1, AYLP_NO_SIGNAL, &a, &b);
	CHECK(a == -1.25 && b == -1.25,
		"parked at [%g, %g], want [-1.25, -1.25]", a, b);
	free_pid(&dev);
}


int main(void)
{
	// the tests deliberately trip the park, and its log lines are expected
	// output rather than failures
	log_init(LOG_ERROR);

	test_controls_with_signal();
	test_short_dropout_does_not_park();
	test_parks_on_threshold();
	test_recovery_starts_from_rest();
	test_counter_resets();
	test_disabled_by_default();
	test_nonzero_park_value();

	if (n_fail)
		fprintf(stderr, "%u check(s) failed\n", n_fail);
	else
		fprintf(stderr, "all pid park checks passed\n");
	return n_fail ? 1 : 0;
}
