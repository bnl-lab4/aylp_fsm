/** Unit tests for anyloop:fsm, the composite steering stage.
*
* Everything here runs on test_source, clamp and file_sink, so there is no
* camera, no DAC and no config file involved: the substages are real devices
* driven exactly as they would be on the bench, and the "mirror" is a file we
* read back afterwards.
*
* What is worth testing about a container:
*
*  - the stages run in ROLE order, not in the order the json keys happen to be
*    written in (a config that lists "actuator" first must not put the DAC
*    before the sensor);
*  - the composite's own type_in/type_out describe the chain, so the outer
*    pipeline is still typechecked, and a chain that cannot work is refused at
*    init rather than at 3 kHz;
*  - the mirror is parked exactly once on the way out, with the shape the
*    actuator has been seeing.
*
* Run: ninja -C build test  (or build/test_fsm) */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "anyloop.h"
#include "logging.h"
#include "fsm.h"

static unsigned n_fail = 0;
static const char *current_test = "";
static char tmpdir[] = "/tmp/aylp_fsm_test_XXXXXX";

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

/** Write the path of a scratch file in our own temp directory into `buf`.
* The caller owns the buffer: two live paths at once is normal here (a tap and
* an actuator record), and a shared static one would silently alias them. */
static const char *scratch(char *buf, size_t n, const char *name)
{
	snprintf(buf, n, "%s/%s", tmpdir, name);
	return buf;
}

static int make_fsm(struct aylp_device *dev, const char *params_json)
{
	memset(dev, 0, sizeof(*dev));
	dev->uri = "anyloop:fsm";
	dev->params = json_tokener_parse(params_json);
	if (!dev->params) {
		fprintf(stderr, "test bug: bad json: %s\n", params_json);
		exit(2);
	}
	return fsm_init(dev);
}

static void free_fsm(struct aylp_device *dev)
{
	if (dev->fini) dev->fini(dev);
	json_object_put(dev->params);
}

/** Run the composite `n` times on a fresh state, returning the last state. */
static struct aylp_state run(struct aylp_device *dev, size_t n)
{
	struct aylp_state state = {0};
	state.header.magic = AYLP_MAGIC;
	state.header.version = AYLP_SCHEMA_VERSION;
	for (size_t i = 0; i < n; i++) {
		int err = dev->proc(dev, &state);
		CHECK(!err, "proc returned %d on iteration %zu", err, i);
	}
	return state;
}

/** Read back an AYLP file of `width`-element double records. Returns the number
* of records and writes their data into `out` (at most `max` records' worth). */
static size_t read_records(const char *path, size_t width, double *out,
	size_t max
){
	FILE *fp = fopen(path, "rb");
	if (!fp) {
		fprintf(stderr, "test bug: can't open %s\n", path);
		exit(2);
	}
	struct aylp_header header;
	size_t n = 0;
	while (fread(&header, sizeof header, 1, fp) == 1) {
		double vals[16];
		if (width > 16) exit(2);
		if (fread(vals, sizeof(double), width, fp) != width) break;
		if (n < max) memcpy(out + n*width, vals, width*sizeof(double));
		n++;
	}
	fclose(fp);
	return n;
}


// ------------------------------------------------------------------- tests

/** The stages run in role order even when the json lists them backwards, and
* the clamp between the source and the sink really is applied. */
static void test_role_order(void)
{
	current_test = "role_order";
	char pathbuf[256];
	const char *path = scratch(pathbuf, sizeof pathbuf, "order.aylp");
	char params[512];
	snprintf(params, sizeof params,
		"{"
		"\"actuator\": {\"uri\": \"anyloop:file_sink\","
			" \"params\": {\"filename\": \"%s\"}},"
		"\"clamp\": {\"min\": -0.4, \"max\": 0.4},"
		"\"sensor\": {\"uri\": \"anyloop:test_source\", \"params\": {"
			"\"type\": \"vector\", \"kind\": \"constant\","
			" \"size1\": 2, \"offset\": 0.9, \"amplitude\": 0}}"
		"}", path
	);
	struct aylp_device dev;
	CHECK(!make_fsm(&dev, params), "init failed");

	// test_source accepts anything and emits a minmax vector, so that is
	// what the composite reports outward
	CHECK(dev.type_in == AYLP_T_ANY, "type_in = 0x%hhX", dev.type_in);
	CHECK(dev.units_in == AYLP_U_ANY, "units_in = 0x%hhX", dev.units_in);
	CHECK(dev.type_out == AYLP_T_VECTOR, "type_out = 0x%hhX", dev.type_out);
	CHECK(dev.units_out == AYLP_U_MINMAX, "units_out = 0x%hhX",
		dev.units_out);

	struct aylp_state state = run(&dev, 3);
	CHECK(state.header.type == AYLP_T_VECTOR, "type 0x%hhX",
		state.header.type);
	CHECK(state.vector && state.vector->size == 2, "no command vector");
	if (state.vector && state.vector->size == 2) {
		// 0.9 out of the source, clamped to 0.4 before the sink: if the
		// stages had run in json order the sink would have recorded the
		// unclamped value and this would read 0.9
		CHECK(gsl_vector_get(state.vector, 0) == 0.4, "cmd0 = %g",
			gsl_vector_get(state.vector, 0));
		CHECK(gsl_vector_get(state.vector, 1) == 0.4, "cmd1 = %g",
			gsl_vector_get(state.vector, 1));
	}
	free_fsm(&dev);

	double recs[8*2];
	size_t n = read_records(path, 2, recs, 8);
	// three iterations plus the park write
	CHECK(n == 4, "recorded %zu records, expected 4", n);
	if (n == 4) {
		for (size_t i = 0; i < 3; i++)
			CHECK(recs[i*2] == 0.4 && recs[i*2+1] == 0.4,
				"record %zu = [%g, %g], expected clamped 0.4",
				i, recs[i*2], recs[i*2+1]);
		CHECK(recs[6] == 0.0 && recs[7] == 0.0,
			"park record = [%g, %g], expected zeros",
			recs[6], recs[7]);
	}
	unlink(path);
}


/** park_value moves the park command; park_on_exit false removes it. */
static void test_park_options(void)
{
	current_test = "park_options";
	const char *fmt =
		"{"
		"\"sensor\": {\"uri\": \"anyloop:test_source\", \"params\": {"
			"\"type\": \"vector\", \"kind\": \"constant\","
			" \"size1\": 2, \"offset\": 0.5, \"amplitude\": 0}},"
		"\"actuator\": {\"uri\": \"anyloop:file_sink\","
			" \"params\": {\"filename\": \"%s\"}},"
		"%s"
		"}";
	char params[512];
	char pathbuf[256];
	struct aylp_device dev;
	double recs[8*2];

	const char *path = scratch(pathbuf, sizeof pathbuf, "park_value.aylp");
	snprintf(params, sizeof params, fmt, path, "\"park_value\": -0.25");
	CHECK(!make_fsm(&dev, params), "init failed (park_value)");
	run(&dev, 2);
	free_fsm(&dev);
	size_t n = read_records(path, 2, recs, 8);
	CHECK(n == 3, "park_value: %zu records, expected 3", n);
	if (n == 3) CHECK(recs[4] == -0.25 && recs[5] == -0.25,
		"park record = [%g, %g], expected -0.25", recs[4], recs[5]);
	unlink(path);

	path = scratch(pathbuf, sizeof pathbuf, "park_off.aylp");
	snprintf(params, sizeof params, fmt, path, "\"park_on_exit\": false");
	CHECK(!make_fsm(&dev, params), "init failed (park_on_exit false)");
	run(&dev, 2);
	free_fsm(&dev);
	n = read_records(path, 2, recs, 8);
	CHECK(n == 2, "park_on_exit false: %zu records, expected 2", n);
	unlink(path);
}


/** A chain whose stages don't fit together is refused at init. */
static void test_typecheck(void)
{
	current_test = "typecheck";
	struct aylp_device dev;
	// a matrix_uchar source into a stage that clamps vectors
	CHECK(make_fsm(&dev,
		"{"
		"\"sensor\": {\"uri\": \"anyloop:test_source\", \"params\": {"
			"\"type\": \"matrix_uchar\", \"kind\": \"constant\","
			" \"size1\": 2, \"size2\": 2, \"offset\": 3}},"
		"\"clamp\": {\"type\": \"vector\", \"min\": -1, \"max\": 1}"
		"}"
	), "a vector clamp after a matrix_uchar source should not initialize");
	free_fsm(&dev);
}


/** Unknown substage devices fail the init of the whole composite. */
static void test_bad_substage(void)
{
	current_test = "bad_substage";
	struct aylp_device dev;
	CHECK(make_fsm(&dev, "{\"controller\": \"anyloop:not_a_device\"}"),
		"an unknown substage uri should not initialize");
	free_fsm(&dev);

	CHECK(make_fsm(&dev, "{}"),
		"a composite with no substages should not initialize");
	free_fsm(&dev);
}


/** A role can hold an array, which is how logging taps get inside: the sink
* between the sensor and the controller sees the measurement, not the command.
* Also covers bare device names (no scheme) and the clamp type default. */
static void test_array_role(void)
{
	current_test = "array_role";
	char tapbuf[128], outbuf[128];
	const char *tap = scratch(tapbuf, sizeof tapbuf, "tap.aylp");
	const char *out = scratch(outbuf, sizeof outbuf, "array_out.aylp");
	char params[1024];
	snprintf(params, sizeof params,
		"{"
		"\"sensor\": ["
			"{\"uri\": \"test_source\", \"params\": {"
				"\"type\": \"vector\", \"kind\": \"constant\","
				" \"size1\": 2, \"offset\": 0.8,"
				" \"amplitude\": 0}},"
			"{\"uri\": \"file_sink\", \"params\":"
				" {\"filename\": \"%s\"}}"
		"],"
		"\"clamp\": {\"min\": -0.2, \"max\": 0.2},"
		"\"actuator\": {\"uri\": \"file_sink\","
			" \"params\": {\"filename\": \"%s\"}}"
		"}", tap, out
	);
	struct aylp_device dev;
	CHECK(!make_fsm(&dev, params), "init failed");
	run(&dev, 2);
	free_fsm(&dev);

	double recs[8*2];
	size_t n = read_records(tap, 2, recs, 8);
	CHECK(n == 2, "tap: %zu records, expected 2 (no park write)", n);
	if (n == 2) CHECK(recs[0] == 0.8 && recs[1] == 0.8,
		"tap saw [%g, %g], expected the unclamped measurement 0.8",
		recs[0], recs[1]);
	n = read_records(out, 2, recs, 8);
	CHECK(n == 3, "actuator: %zu records, expected 2 plus the park", n);
	if (n == 3) {
		CHECK(recs[0] == 0.2, "actuator saw %g, expected clamped 0.2",
			recs[0]);
		CHECK(recs[4] == 0.0 && recs[5] == 0.0,
			"park record = [%g, %g], expected zeros",
			recs[4], recs[5]);
	}
	unlink(tap);
	unlink(out);
}


int main(void)
{
	log_init(LOG_WARN);
	if (!mkdtemp(tmpdir)) {
		fprintf(stderr, "test bug: mkdtemp failed\n");
		return 2;
	}
	test_role_order();
	test_park_options();
	test_typecheck();
	test_bad_substage();
	test_array_role();
	rmdir(tmpdir);
	if (n_fail) {
		fprintf(stderr, "%u check(s) failed\n", n_fail);
		return 1;
	}
	printf("all fsm tests passed\n");
	return 0;
}
