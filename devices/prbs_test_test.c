/** Unit tests for prbs_test, against a simulated plant of known delay.
*
* A delay measurement is worth exactly as much as its ability to return the
* delay you put in, so that is what these check: a plant with a programmed
* whole-frame delay and a one-pole rise is driven by the device's own PRBS, and
* the reported onset, peak and phase-slope numbers are compared against the
* delay the plant was built with. The awkward cases are covered too -- a plant
* whose gain is negative (the optical path inverts on one axis here), a
* stimulus buried in noise, and a plant that isn't moving at all, which must be
* reported as a failure rather than as a fast loop.
*
* The device counts frames rather than seconds for everything except its
* warmup, so these run at whatever rate the machine manages; a realistic loop
* period is slept anyway, since the ms and Hz figures come from wall clock.
*
* Run: ninja -C build test  (or build/test_prbs_test) */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "anyloop.h"
#include "logging.h"
#include "prbs_test.h"

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

#define LOOP_NS 120000
#define MAXDELAY 16

/** Pure delay of `delay` frames followed by a one-pole rise, which is what the
 * loop's feedback path looks like: transport (DAC write, actuator dead time,
 * exposure, readout, centroid) and then the actuator's own motion. */
struct plant {
	double gain;
	unsigned delay;		// whole frames of pure transport delay
	double pole;		// one-pole lag, 0 = the plant is instant
	double noise;
	int dead;		// plant does not respond at all
	double hist[MAXDELAY];
	size_t n;
	double y;
	unsigned seed;
	// A chopped source, modelled the way a centroid stage behaves on one:
	// for `chop_dark` frames out of every `chop_period` the sensor cannot
	// fit the beam, so it re-publishes its previous coordinate and raises
	// AYLP_FRAME_REJECTED. 0 = the source is CW.
	size_t chop_period, chop_dark;
	size_t frames;
	double held;
};

static double plant_step(struct plant *p, double u)
{
	double old = p->hist[p->n % MAXDELAY];
	p->hist[(p->n + p->delay) % MAXDELAY] = u;
	p->n++;
	if (p->dead) old = 0.0;
	p->y = p->pole * p->y + (1.0 - p->pole) * p->gain * old;
	double out = p->y;
	if (p->noise > 0.0) {
		// two draws, so the noise is not flat-spectrum trivial
		p->seed = p->seed * 1103515245u + 12345u;
		double a = (p->seed >> 16) / 65536.0;
		p->seed = p->seed * 1103515245u + 12345u;
		double b = (p->seed >> 16) / 65536.0;
		out += p->noise * (a + b - 1.0);
	}
	return out;
}

static size_t run(struct aylp_device *dev, struct plant *p, size_t max_iter)
{
	struct timespec sl = { 0, LOOP_NS };
	double u = 0.0;
	size_t i = 0;
	for (; i < max_iter; i++) {
		double y = plant_step(p, u);
		bool dark = p->chop_period
			&& p->frames % p->chop_period < p->chop_dark;
		if (dark) y = p->held;
		else p->held = y;
		p->frames++;
		gsl_vector *in = gsl_vector_alloc(2);
		gsl_vector_set(in, 0, y);
		gsl_vector_set(in, 1, 0.0);
		struct aylp_state state = {0};
		state.header.type = AYLP_T_VECTOR;
		state.header.units = AYLP_U_MINMAX;
		if (dark) state.header.status |= AYLP_FRAME_REJECTED;
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
	dev->uri = "anyloop:prbs_test";
	dev->params = json_tokener_parse(params_json);
	if (!dev->params) {
		fprintf(stderr, "test bug: bad json: %s\n", params_json);
		exit(2);
	}
	return prbs_test_init(dev);
}

static void unmake(struct aylp_device *dev)
{
	if (dev->fini) dev->fini(dev);
	json_object_put(dev->params);
}

// order 7 = 127 chips, one frame per chip, 16 bursts: ~5000 iterations
#define BURSTS "{\"index_cmd\":0,\"index_err\":0,\"out_size\":2," \
	"\"amplitude\":0.05,\"order\":7,\"chip_frames\":1,\"n_bursts\":16," \
	"\"quiet_frames\":96,\"warmup\":0.02,\"max_lag\":32,\"neg_lags\":16," \
	"\"phase_f_lo\":100,\"phase_f_hi\":2000,\"output_file\":\"\""


/** The headline claim: the reported delay is the plant's delay. A plant delay
 * of 3 frames is a command-to-response lag of 4 (the harness applies the
 * command the iteration after proc published it, as the real pipeline does --
 * see test_pure_delay_is_exact). Onset reads slightly early of that, the peak
 * sits past it by the rise, and neither is below the one frame that is
 * physically possible. */
static void test_recovers_delay(void)
{
	current_test = "recovers_delay";
	struct aylp_device dev;
	CHECK(!make(&dev, BURSTS "}"), "init failed");
	struct plant p = { .gain = 4.0, .delay = 3, .pole = 0.5,
		.noise = 0.01, .seed = 1 };
	run(&dev, &p, 100000);
	struct aylp_prbs_test_data *d = dev.device_data;
	CHECK(d->stage == PT_DONE, "did not finish (stage %d)", d->stage);
	CHECK(d->have_result, "no result produced");
	CHECK(d->burst == 16, "recorded %zu bursts, want 16", d->burst);
	CHECK(d->lag_onset >= 1.0, "onset %g is below one frame", d->lag_onset);
	CHECK(d->lag_onset > 2.5 && d->lag_onset < 4.1,
		"onset %g frames, want just under the 4-frame lag",
		d->lag_onset);
	CHECK(d->lag_peak >= d->lag_onset - 0.01,
		"peak %g is before onset %g", d->lag_peak, d->lag_onset);
	CHECK(d->lag_peak > 3.9 && d->lag_peak < 7.0,
		"peak %g frames, want the 4-frame lag plus a one-pole rise",
		d->lag_peak);
	CHECK(d->rho_peak > 0.5, "peak correlation %g is weak", d->rho_peak);
	// the floor here is not just sensor noise: the linear (rather than
	// circular) autocorrelation of a finite m-sequence has sidelobes of
	// order 1/sqrt(n_chips), a few percent of the peak
	CHECK(fabs(d->rho_peak) > 8.0 * d->rho_noise,
		"peak %g is not clear of the %g noise floor", d->rho_peak,
		d->rho_noise);
	unmake(&dev);
}

/** Against a plant that is nothing but a pure delay, the reported lag is that
 * delay exactly -- no hidden frame either way. This is the check that keeps
 * this device's numbers comparable with bode_plot's `fit_frames` (which has to
 * add one back, because its correlator references the previous iteration's
 * excitation phase) and directly usable as an fsp delay/delay_frac pair.
 *
 * The harness applies the command one iteration after proc published it, which
 * is what the real pipeline does too, so a plant delay of D frames is a
 * command-to-response lag of D+1. */
static void test_pure_delay_is_exact(void)
{
	current_test = "pure_delay_is_exact";
	for (unsigned D = 1; D <= 4; D++) {
		struct aylp_device dev;
		CHECK(!make(&dev, "{\"index_cmd\":0,\"index_err\":0,"
			"\"out_size\":2,\"amplitude\":0.05,\"order\":7,"
			"\"n_bursts\":4,\"quiet_frames\":128,\"warmup\":0,"
			"\"max_lag\":32,\"neg_lags\":16,\"phase_f_lo\":50,"
			"\"phase_f_hi\":800,\"output_file\":\"\"}"),
			"init failed");
		// pole 0 and no noise: the plant IS the delay
		struct plant p = { .gain = 4.0, .delay = D, .seed = 21 };
		run(&dev, &p, 100000);
		struct aylp_prbs_test_data *d = dev.device_data;
		double want = D + 1.0;
		CHECK(d->have_result, "no result for a pure %u-frame delay", D);
		CHECK(fabs(d->lag_peak - want) < 0.05,
			"peak %g for a lag of %g frames", d->lag_peak, want);
		double frames = d->tau_phase_ms / (1e3 / d->fs);
		CHECK(fabs(frames - want) < 0.05,
			"phase slope %g frames for a lag of %g", frames, want);
		// onset reads early: the m-sequence autocorrelation is a
		// triangle a chip wide either side, so the correlogram starts
		// rising a chip before the true lag
		CHECK(d->lag_onset < want && d->lag_onset > want - 1.2,
			"onset %g for a lag of %g frames -- expected slightly "
			"early, within one chip", d->lag_onset, want);
		unmake(&dev);
	}
}

/** The correlogram is causal: the acausal lags hold nothing but noise. If they
 * did hold signal, the command and the response would be misaligned and every
 * delay in the report would be off by that much. */
static void test_no_acausal_response(void)
{
	current_test = "no_acausal_response";
	struct aylp_device dev;
	CHECK(!make(&dev, BURSTS "}"), "init failed");
	struct plant p = { .gain = 4.0, .delay = 3, .pole = 0.5,
		.noise = 0.005, .seed = 2 };
	run(&dev, &p, 100000);
	struct aylp_prbs_test_data *d = dev.device_data;
	CHECK(d->have_result, "no result produced");
	double worst = 0.0;
	for (size_t i = 0; i < d->neg_lags; i++)
		if (fabs(d->rho[i]) > worst) worst = fabs(d->rho[i]);
	CHECK(worst < 0.2 * fabs(d->rho_peak),
		"an acausal lag reached %g against a peak of %g", worst,
		d->rho_peak);
	// lag 0 is likewise impossible: this iteration's sensor value was
	// captured before this iteration's command went out
	CHECK(fabs(d->rho[d->neg_lags]) < 0.5 * fabs(d->rho_peak),
		"lag 0 carries %g of the peak %g", d->rho[d->neg_lags],
		d->rho_peak);
	unmake(&dev);
}

/** The phase-slope group delay agrees with the correlogram. This is the number
 * meant to be compared with a bode sweep's tau and put into an fsp
 * delay/delay_frac pair, so it has to be right, and sub-frame. */
static void test_phase_slope_agrees(void)
{
	current_test = "phase_slope_agrees";
	struct aylp_device dev;
	CHECK(!make(&dev, BURSTS "}"), "init failed");
	struct plant p = { .gain = 4.0, .delay = 3, .pole = 0.5,
		.noise = 0.005, .seed = 4 };
	run(&dev, &p, 100000);
	struct aylp_prbs_test_data *d = dev.device_data;
	CHECK(d->have_result, "no result produced");
	CHECK(d->phase_n > 4, "phase fit used %zu bins", d->phase_n);
	double dt_ms = 1e3 / d->fs;
	double frames = d->tau_phase_ms / dt_ms;
	// a pole of 0.5 adds about one frame of group delay on top of the 3
	CHECK(frames > 2.5 && frames < 6.0,
		"phase-slope delay %g frames (%g ms at %g Hz), want about 4",
		frames, d->tau_phase_ms, d->fs);
	CHECK(fabs(frames - d->lag_peak) < 2.0,
		"phase-slope %g and peak %g disagree", frames, d->lag_peak);
	unmake(&dev);
}

/** One axis of this rig inverts optically, so the plant's gain is negative.
 * The delay must come out the same; only the sign of the correlation changes. */
static void test_negative_gain(void)
{
	current_test = "negative_gain";
	struct aylp_device dev;
	CHECK(!make(&dev, BURSTS "}"), "init failed");
	struct plant p = { .gain = -4.0, .delay = 3, .pole = 0.5,
		.noise = 0.005, .seed = 6 };
	run(&dev, &p, 100000);
	struct aylp_prbs_test_data *d = dev.device_data;
	CHECK(d->have_result, "no result produced");
	CHECK(d->rho_peak < 0.0, "peak correlation %g, want negative",
		d->rho_peak);
	CHECK(d->lag_onset > 2.5 && d->lag_onset < 4.1,
		"onset %g frames, want just under the 4-frame lag",
		d->lag_onset);
	CHECK(d->lag_peak > 3.0 && d->lag_peak < 7.0, "peak %g frames",
		d->lag_peak);
	double frames = d->tau_phase_ms / (1e3 / d->fs);
	CHECK(frames > 2.5 && frames < 6.0,
		"phase-slope delay %g frames on an inverting plant", frames);
	unmake(&dev);
}

/** Averaging bursts is what buys the SNR: the same measurement under noise
 * several times the response still lands on the right delay. */
static void test_survives_noise(void)
{
	current_test = "survives_noise";
	struct aylp_device dev;
	CHECK(!make(&dev, "{\"index_cmd\":0,\"index_err\":0,\"out_size\":2,"
		"\"amplitude\":0.05,\"order\":7,\"n_bursts\":48,"
		"\"quiet_frames\":96,\"warmup\":0.02,\"max_lag\":32,"
		"\"neg_lags\":16,\"output_file\":\"\"}"), "init failed");
	// response to a chip is ~0.1 units; the noise here is 5x that
	struct plant p = { .gain = 4.0, .delay = 3, .pole = 0.5,
		.noise = 0.5, .seed = 8 };
	run(&dev, &p, 200000);
	struct aylp_prbs_test_data *d = dev.device_data;
	CHECK(d->have_result, "no result under noise");
	CHECK(fabs(d->lag_onset - 3.0) <= 1.5,
		"onset %g frames under noise, want about 3", d->lag_onset);
	CHECK(d->lag_peak > 2.5 && d->lag_peak < 8.0,
		"peak %g frames under noise", d->lag_peak);
	// per-burst estimates are noisy by construction; their median is not
	CHECK(d->n_burst_ok > 0, "no single burst produced an estimate");
	CHECK(fabs(d->lag_peak_med - d->lag_peak) < 2.0,
		"per-burst median %g disagrees with the ensemble peak %g",
		d->lag_peak_med, d->lag_peak);
	unmake(&dev);
}

/** A plant that is not moving must be reported as a failure. The dangerous
 * outcome here is a confident small delay read off the noise. */
static void test_dead_plant_is_reported(void)
{
	current_test = "dead_plant_is_reported";
	struct aylp_device dev;
	CHECK(!make(&dev, BURSTS "}"), "init failed");
	struct plant p = { .gain = 4.0, .delay = 3, .pole = 0.5,
		.noise = 0.05, .dead = 1, .seed = 9 };
	run(&dev, &p, 100000);
	struct aylp_prbs_test_data *d = dev.device_data;
	CHECK(!d->have_result,
		"reported a delay of %g frames from a plant that never moved",
		d->lag_peak);
	unmake(&dev);
}

/** Finding a correlation lobe is necessary but no longer sufficient for a
 * bench calibration. Strict runs reject poor phase fits and inadequate live
 * coverage while preserving the candidate diagnostics for troubleshooting. */
static void test_quality_gates(void)
{
	current_test = "quality_gates_pass";
	struct aylp_device dev;
	CHECK(!make(&dev, BURSTS ",\"min_peak_snr\":8,"
		"\"max_phase_resid_deg\":20,\"min_phase_bins\":4,"
		"\"min_burst_frac\":0.5,\"max_peak_mad\":0.5,"
		"\"min_live_frac\":0.7,\"max_holes\":0,"
		"\"max_phase_peak_delta\":2}"), "init failed");
	struct plant p = { .gain = 4.0, .delay = 3, .pole = 0.5,
		.noise = 0.005, .seed = 19 };
	run(&dev, &p, 100000);
	struct aylp_prbs_test_data *d = dev.device_data;
	CHECK(d->have_result, "strong run found no correlation");
	CHECK(d->quality_ok, "strong run failed quality gates: %s",
		d->quality_reason);
	unmake(&dev);

	current_test = "quality_gates_reject_phase";
	CHECK(!make(&dev, BURSTS ",\"max_phase_resid_deg\":0.000001}"),
		"init failed");
	p = (struct plant){ .gain = 4.0, .delay = 3, .pole = 0.5,
		.noise = 0.01, .seed = 20 };
	run(&dev, &p, 100000);
	d = dev.device_data;
	CHECK(d->have_result, "candidate correlation was lost");
	CHECK(!d->quality_ok, "accepted a deliberately impossible phase gate");
	CHECK(strstr(d->quality_reason, "phase residual") != NULL,
		"quality failure did not identify phase residual: %s",
		d->quality_reason);
	unmake(&dev);

	current_test = "quality_gates_reject_live_fraction";
	CHECK(!make(&dev, BURSTS ",\"min_live_frac\":0.9}"), "init failed");
	p = (struct plant){ .gain = 4.0, .delay = 3, .pole = 0.5,
		.noise = 0.01, .seed = 22, .chop_period = 51, .chop_dark = 13 };
	run(&dev, &p, 100000);
	d = dev.device_data;
	CHECK(d->have_result, "candidate correlation was lost through chop");
	CHECK(!d->quality_ok, "accepted insufficient live-frame coverage");
	CHECK(strstr(d->quality_reason, "live frames") != NULL,
		"quality failure did not identify live coverage: %s",
		d->quality_reason);
	unmake(&dev);
}

/** Every tap mask in the table really does give a maximum-length sequence: a
 * short one would put a periodic ridge in the autocorrelation and could be
 * mistaken for a delay. init walks the whole cycle to check, so a bad mask
 * fails here rather than in a measurement. */
static void test_all_orders_are_maximal(void)
{
	current_test = "all_orders_are_maximal";
	for (unsigned order = 5; order <= 16; order++) {
		char params[256];
		snprintf(params, sizeof params,
			"{\"index_cmd\":0,\"index_err\":0,\"out_size\":2,"
			"\"order\":%u,\"n_bursts\":1,\"quiet_frames\":128,"
			"\"output_file\":\"\"}", order);
		struct aylp_device dev;
		int err = make(&dev, params);
		CHECK(!err, "order %u rejected by init", order);
		if (!err) {
			struct aylp_prbs_test_data *d = dev.device_data;
			CHECK(d->n_chips == ((size_t)1 << order) - 1,
				"order %u gave %zu chips", order, d->n_chips);
			// a maximum-length sequence is balanced to one chip
			long sum = 0;
			for (size_t i = 0; i < d->n_chips; i++)
				sum += d->chips[i];
			CHECK(labs(sum) == 1, "order %u sequence sums to %ld, "
				"want +-1", order, sum);
		}
		unmake(&dev);
	}
}

/** The reason the masking exists. A chopped source makes the sensor hold its
 * last coordinate through the dark frames, so a quarter of the stream carries
 * values several frames old.
 *
 * What that costs, measured over this simulation at 25% dark in 13-frame
 * streaks: the PEAK barely moves (+0.013 fr) but the phase slope -- the number
 * this device exists to hand to an fsp delay/delay_frac pair -- reads short,
 * and the correlation drops from 0.865 to 0.836. At 51% dark it is 0.66 frames
 * short at rho 0.73. Short, not long: held samples have a spread of ages, so
 * they decorrelate the response rather than shifting it coherently, and what a
 * decorrelated high-frequency end does to a phase-slope fit is flatten it.
 * (Below ~10% dark the sign is not even stable, for the same reason: at that
 * point it is all decorrelation and no shift.)
 *
 * How far short depends on where [phase_f_lo, phase_f_hi] sits in the spectrum,
 * which moves with the loop rate the machine happens to manage -- 0.09 frames
 * at this harness's rate, 0.24 at a faster one. So the thresholds below check
 * a direction and a floor, not a value.
 *
 * So the test is two-sided. Masked, the answer has to match what the same
 * plant reads on a CW source; unmasked, it has to be visibly worse, or the
 * masking would be a no-op nobody could check. */
static void test_chop_does_not_bias_lag(void)
{
	current_test = "chop_does_not_bias_lag";
	// 13 dark frames in every 51 = 25% held, the duty this rig measured.
	// 51 is coprime with the 223-frame window, so the chop walks across
	// the burst instead of dimming the same positions every time.
	struct plant base = { .gain = 4.0, .delay = 3, .pole = 0.5,
		.noise = 0.01, .seed = 11 };
	struct plant clean = base, chopped = base, biased = base;
	chopped.chop_period = biased.chop_period = 51;
	chopped.chop_dark = biased.chop_dark = 13;

	struct aylp_device dev;
	CHECK(!make(&dev, BURSTS "}"), "init failed");
	run(&dev, &clean, 100000);
	struct aylp_prbs_test_data *d = dev.device_data;
	CHECK(d->have_result, "no result on a CW source");
	double cw_peak = d->lag_peak, cw_rho = d->rho_peak;
	double cw_cent = d->lag_centroid;
	double cw_tau = d->tau_phase_ms / (1e3 / d->fs);
	CHECK(d->n_live == d->n_frames,
		"a CW source reported %zu of %zu frames live",
		d->n_live, d->n_frames);
	unmake(&dev);

	CHECK(!make(&dev, BURSTS "}"), "init failed");
	run(&dev, &chopped, 100000);
	d = dev.device_data;
	CHECK(d->have_result, "no result through the chop");
	// the duty cycle is measured, not assumed: 38/51 frames are live
	double live = (double)d->n_live / (double)d->n_frames;
	CHECK(fabs(live - 38.0/51.0) < 0.02,
		"tracked live fraction %g, want the %g the harness chopped",
		live, 38.0/51.0);
	CHECK(!d->n_holes, "%zu holes from a chop that walks the window",
		d->n_holes);
	// the chop should leave no trace at all: the ensemble mean is built
	// from live samples only, and 16 bursts of a walking chop cover every
	// window position many times over
	CHECK(fabs(d->lag_peak - cw_peak) < 0.05,
		"masked peak %g frames against %g on a CW source",
		d->lag_peak, cw_peak);
	double tau = d->tau_phase_ms / (1e3 / d->fs);
	CHECK(fabs(tau - cw_tau) < 0.15,
		"masked phase-slope tau %g frames against %g on a CW source",
		tau, cw_tau);
	CHECK(fabs(d->rho_peak - cw_rho) < 0.02,
		"masked correlation %g against %g on a CW source -- the chop "
		"should not cost peak height either", d->rho_peak, cw_rho);
	unmake(&dev);

	// same data correlated the old way, so the bias being guarded against
	// is exercised rather than assumed
	CHECK(!make(&dev, BURSTS ",\"use_rejected\":true}"), "init failed");
	run(&dev, &biased, 100000);
	d = dev.device_data;
	CHECK(d->have_result, "no result with use_rejected");
	// The two checks that do not depend on the loop rate: lags are counted
	// in frames and rho is a correlation, so neither moves with whatever
	// speed the machine managed. The phase slope shifts as well, but by an
	// amount set by where the fit band lands in the spectrum, so asserting
	// on it here would be asserting on the test machine.
	CHECK(d->lag_centroid - cw_cent > 0.01,
		"correlating held samples put the lobe centroid at %g frames "
		"against the CW %g; the bias this guards against did not "
		"reproduce, so the guard is untested",
		d->lag_centroid, cw_cent);
	CHECK(cw_rho - d->rho_peak > 0.01,
		"correlating held samples gave rho %g against the CW %g; the "
		"decorrelation this guards against did not reproduce",
		d->rho_peak, cw_rho);
	unmake(&dev);
}

/** A chop locked to the burst period is the case the ensemble mean cannot
 * repair: the same window positions are dark in every burst. Those become
 * holes, the run says so, and the correlogram skips them rather than
 * correlating an interpolated guess. */
static void test_chop_locked_to_window_is_reported(void)
{
	current_test = "chop_locked_to_window_is_reported";
	struct aylp_device dev;
	CHECK(!make(&dev, BURSTS "}"), "init failed");
	// window is 127 + 96 = 223 frames; a chop of exactly that period dims
	// the same 13 positions in every burst
	struct plant p = { .gain = 4.0, .delay = 3, .pole = 0.5,
		.noise = 0.01, .seed = 12, .chop_period = 223, .chop_dark = 13 };
	run(&dev, &p, 100000);
	struct aylp_prbs_test_data *d = dev.device_data;
	CHECK(d->n_holes == 13,
		"%zu window positions with no live sample, want the 13 dimmed "
		"in every burst", d->n_holes);
	CHECK(d->worst_cnt == 0, "worst position had %zu live bursts",
		d->worst_cnt);
	CHECK(d->have_result, "no result with 13 of 223 positions dark");
	CHECK(d->lag_peak > 2.5 && d->lag_peak < 7.0,
		"peak %g frames with a locked chop", d->lag_peak);
	unmake(&dev);
}

/** Configuration that cannot produce a measurement is refused at init. */
static void test_rejects_bad_config(void)
{
	current_test = "rejects_bad_config";
	struct aylp_device dev;
	CHECK(make(&dev, "{\"order\":4}"), "accepted order 4");
	unmake(&dev);
	CHECK(make(&dev, "{\"order\":17}"), "accepted order 17");
	unmake(&dev);
	CHECK(make(&dev, "{\"amplitude\":0}"), "accepted amplitude 0");
	unmake(&dev);
	// lags longer than the window would correlate off the end of it
	CHECK(make(&dev, "{\"order\":5,\"quiet_frames\":0,\"max_lag\":64,"
		"\"neg_lags\":32}"), "accepted lags longer than the window");
	unmake(&dev);
	CHECK(make(&dev, "{\"index_cmd\":9,\"out_size\":2}"),
		"accepted index_cmd outside the output vector");
	unmake(&dev);
}


int main(void)
{
	// the dead-plant test deliberately trips an error log; it is expected
	// output rather than a failure
	log_init(LOG_FATAL);

	test_recovers_delay();
	test_pure_delay_is_exact();
	test_no_acausal_response();
	test_phase_slope_agrees();
	test_negative_gain();
	test_survives_noise();
	test_dead_plant_is_reported();
	test_quality_gates();
	test_chop_does_not_bias_lag();
	test_chop_locked_to_window_is_reported();
	test_all_orders_are_maximal();
	test_rejects_bad_config();

	if (n_fail) fprintf(stderr, "%u check(s) failed\n", n_fail);
	else fprintf(stderr, "all prbs_test checks passed\n");
	return n_fail ? 1 : 0;
}
