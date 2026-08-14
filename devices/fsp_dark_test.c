/** Closed-loop integration test for FSP rejected-frame prediction.
 *
 * This calls the real fsp_proc (not contrib/fsp_replay) with the intended
 * 2310 Hz / 512-tap / broad_lp=7 steering parameters.  The simulated sensor
 * behaves like fit_com: on a rejected frame it republishes the last live
 * coordinate while the physical disturbance and delayed actuator continue.
 *
 * Besides the 10-off-in-50 chop, random blank frames are injected
 * only in the nominally live part of the duty cycle.  Results compare the old
 * baseline-mask/held-command path with direct dark-frame prediction.
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_vector.h>
#include <json-c/json.h>

#include "anyloop.h"
#include "fsp.h"

#define FS 2310.0
#define DELAY 3
#define PERIOD 50
#define DUTY_DARK 10
#define TRAIN_FRAMES (25 * 2310)
#define SCORE_FRAMES (20 * 2310)
#define PIXEL_SCALE 31.5
#define COMMAND_LIMIT 5.0

static unsigned n_fail;

#define CHECK(c, ...) do { if (!(c)) { \
	n_fail++; fprintf(stderr, "FAIL %s:%d: ", __FILE__, __LINE__); \
	fprintf(stderr, __VA_ARGS__); fputc('\n', stderr); \
} } while (0)

struct result {
	double rms, live_rms, dark_rms, duty_rms, random_rms;
	double disturbance_rms, peak, cmd_rms, cmd_peak, wnorm;
	size_t n_live, n_dark, n_duty_dark, n_random_dark, n_predicted;
	size_t n_clamped;
};

struct replay_data {
	size_t n;
	double *v[2];
	unsigned char *dark;
};

static unsigned long rng;

static double uniform01(void)
{
	rng = rng * 6364136223846793005UL + 1442695040888963407UL;
	return (double)(rng >> 11) / 9007199254740992.0;
}

static double gaussian(void)
{
	double u1 = uniform01();
	double u2 = uniform01();
	if (u1 < 1e-15) u1 = 1e-15;
	return sqrt(-2.0 * log(u1)) * cos(2.0*M_PI*u2);
}

static int make_fsp_config(struct aylp_device *dev, bool predict,
	bool scheduled_mu)
{
	char json[2048];
	snprintf(json, sizeof json,
		"{\"type\":\"vector\",\"units\":\"V\",\"fs\":2310,"
		"\"delay\":3,\"delay_frac\":0,\"broad_order\":512,"
		"\"broad_mu\":0.03,%s\"broad_lp\":7,"
		"\"broad_freeze_closed\":false,\"dark_predict\":%s,"
		"\"dark_predict_max\":16,\"dark_bank_train\":4,"
		"\"dark_train_min_real\":0.5,\"start_delay\":0,\"ramp\":0,"
		"\"guard_ratio\":0,\"adapt_period\":0,"
		"\"clamp_min\":-5.0,\"clamp_max\":5.0,"
		"\"y\":{\"K\":1,\"r\":1,\"freqs\":[20],"
		"\"zeta\":[0.1],\"q\":[1e-5]},"
		"\"x\":{\"K\":1,\"r\":1,\"freqs\":[20],"
		"\"zeta\":[0.1],\"q\":[1e-5]}}",
		scheduled_mu ? "\"broad_mu_init\":0.5,\"broad_mu_tau\":60," : "",
		predict ? "true" : "false");
	memset(dev, 0, sizeof *dev);
	dev->uri = "anyloop:fsp";
	dev->params = json_tokener_parse(json);
	if (!dev->params) return -1;
	return fsp_init(dev);
}

static int make_fsp(struct aylp_device *dev, bool predict)
{
	return make_fsp_config(dev, predict, false);
}

static void free_fsp(struct aylp_device *dev)
{
	if (dev->fini) dev->fini(dev);
	json_object_put(dev->params);
}

static struct result run_case(bool predict, double random_blank_prob,
	size_t duty_dark_frames, double target_rms_px, unsigned long seed)
{
	struct result out = {0};
	struct aylp_device dev;
	CHECK(!make_fsp(&dev, predict), "fsp init failed");
	if (!dev.device_data) return out;

	const size_t total = TRAIN_FRAMES + SCORE_FRAMES;
	double delayed[2][DELAY] = {{0}};
	double last_meas[2] = {0};
	double ar[2] = {0};
	double e2 = 0, phi2 = 0, live2 = 0, dark2 = 0;
	double duty2 = 0, random2 = 0;
	double u2 = 0;
	/* RMS of the six mutually non-coherent tones plus stationary AR(1). */
	const double base_rms = sqrt((0.13*0.13 + 0.11*0.11
		+ 0.12*0.12 + 0.07*0.07 + 0.055*0.055
		+ 0.035*0.035)/2.0
		+ (0.0015*0.0015)/(1.0 - 0.997*0.997));
	const double jitter_scale = target_rms_px / (PIXEL_SCALE*base_rms);
	rng = seed;

	for (size_t k = 0; k < total; k++) {
		double t = (double)k / FS;
		double phi[2];
		for (int a = 0; a < 2; a++) {
			ar[a] = 0.997*ar[a] + 0.0015*gaussian();
			double p = a ? 0.37 : 1.11;
			phi[a] = jitter_scale*(0.13*sin(2*M_PI*9.2*t + p)
				+ 0.11*sin(2*M_PI*11.8*t + 0.7*p)
				+ 0.12*sin(2*M_PI*21.2*t + 1.3*p)
				+ 0.07*sin(2*M_PI*31.0*t + 0.4*p)
				+ 0.055*sin(2*M_PI*57.5*t + 1.8*p)
				+ 0.035*sin(2*M_PI*83.7*t + 0.2*p)
				+ ar[a]);
		}

		size_t phase = k % PERIOD;
		bool duty_dark = duty_dark_frames > 0
			&& phase >= PERIOD - duty_dark_frames;
		bool random_dark = !duty_dark
			&& uniform01() < random_blank_prob;
		bool dark = duty_dark || random_dark;
		double truth[2];
		for (int a = 0; a < 2; a++) {
			truth[a] = phi[a] + delayed[a][k % DELAY];
			if (!dark) last_meas[a] = truth[a];
		}

		double in_data[2] = {last_meas[0], last_meas[1]};
		gsl_vector in = {.size = 2, .stride = 1, .data = in_data,
			.block = NULL, .owner = 0};
		struct aylp_state state = {0};
		state.header.type = AYLP_T_VECTOR;
		state.header.units = AYLP_U_MINMAX;
		state.header.status = dark ? AYLP_FRAME_REJECTED : 0;
		state.vector = &in;
		int err = dev.proc(&dev, &state);
		CHECK(!err, "fsp proc failed at frame %zu", k);
		if (err) break;

		for (int a = 0; a < 2; a++) {
			double u = state.vector->data[a * state.vector->stride];
			CHECK(isfinite(u), "non-finite command at frame %zu", k);
			CHECK(fabs(u) <= COMMAND_LIMIT + 1e-6,
				"command escaped clamp at frame %zu: %g", k, u);
			delayed[a][k % DELAY] = u;
			if (k < TRAIN_FRAMES) continue;
			double ae = fabs(truth[a]);
			phi2 += phi[a]*phi[a];
			e2 += truth[a]*truth[a];
			u2 += u*u;
			if (ae > out.peak) out.peak = ae;
			if (fabs(u) > out.cmd_peak) out.cmd_peak = fabs(u);
			out.n_clamped += fabs(u) >= COMMAND_LIMIT - 1e-6;
			if (dark) dark2 += truth[a]*truth[a];
			else live2 += truth[a]*truth[a];
			if (duty_dark) duty2 += truth[a]*truth[a];
			else if (random_dark) random2 += truth[a]*truth[a];
		}
		if (k >= TRAIN_FRAMES) {
			if (dark) out.n_dark++; else out.n_live++;
			out.n_duty_dark += duty_dark;
			out.n_random_dark += random_dark;
		}
	}

	out.rms = sqrt(e2 / (2.0*SCORE_FRAMES));
	out.disturbance_rms = sqrt(phi2 / (2.0*SCORE_FRAMES));
	out.cmd_rms = sqrt(u2 / (2.0*SCORE_FRAMES));
	out.live_rms = sqrt(live2 / (2.0*out.n_live));
	out.dark_rms = out.n_dark ? sqrt(dark2 / (2.0*out.n_dark)) : 0.0;
	out.duty_rms = out.n_duty_dark
		? sqrt(duty2 / (2.0*out.n_duty_dark)) : 0.0;
	out.random_rms = out.n_random_dark
		? sqrt(random2 / (2.0*out.n_random_dark)) : 0.0;
	CHECK(out.n_dark == out.n_duty_dark + out.n_random_dark,
		"random blanks overlapped duty-cycle off frames");
	struct aylp_fsp_data *d = dev.device_data;
	out.n_predicted = d->axis[0].n_dark_predicted;
	double wn = 0.0;
	for (int a = 0; a < 2; a++) for (size_t i = 0; i < d->broad_order; i++)
		wn += d->axis[a].broad_w[i] * d->axis[a].broad_w[i];
	out.wnorm = sqrt(wn);
	CHECK(isfinite(out.wnorm) && out.wnorm < 100.0,
		"unbounded coefficient norm %g", out.wnorm);
	free_fsp(&dev);
	return out;
}

static void print_result(const char *name, struct result r)
{
	printf("%-22s input %.3f px  rms %.3f px  live %.3f px  "
		"duty %.3f px  random %.3f px  peak %.3f px  "
		"u_rms %.6f V  u_peak %.6f V  sat %.1f%%  nrand %zu  "
		"predicted %zu  ||w|| %.4f\n", name,
		r.disturbance_rms*PIXEL_SCALE, r.rms*PIXEL_SCALE,
		r.live_rms*PIXEL_SCALE, r.duty_rms*PIXEL_SCALE,
		r.random_rms*PIXEL_SCALE, r.peak*PIXEL_SCALE, r.cmd_rms,
		r.cmd_peak, 100.0*r.n_clamped/(2.0*SCORE_FRAMES),
		r.n_random_dark, r.n_predicted, r.wnorm);
}

static double run_tone(bool predict, size_t duty_dark_frames, double freq)
{
	struct aylp_device dev;
	CHECK(!make_fsp(&dev, predict), "fsp init failed in tone sweep");
	if (!dev.device_data) return NAN;
	const size_t train = 15*2310, score = 5*2310;
	const double amp = sqrt(2.0)*3.0/PIXEL_SCALE;
	double delayed[2][DELAY] = {{0}}, last_meas[2] = {0};
	double e2 = 0.0;

	for (size_t k = 0; k < train + score; k++) {
		double t = (double)k/FS;
		size_t phase = k % PERIOD;
		bool dark = duty_dark_frames > 0
			&& phase >= PERIOD - duty_dark_frames;
		double truth[2];
		for (int a = 0; a < 2; a++) {
			double phi = amp*sin(2*M_PI*freq*t + a*0.71);
			truth[a] = phi + delayed[a][k % DELAY];
			if (!dark) last_meas[a] = truth[a];
		}
		double in_data[2] = {last_meas[0], last_meas[1]};
		gsl_vector in = {.size = 2, .stride = 1, .data = in_data,
			.block = NULL, .owner = 0};
		struct aylp_state state = {0};
		state.header.type = AYLP_T_VECTOR;
		state.header.units = AYLP_U_MINMAX;
		state.header.status = dark ? AYLP_FRAME_REJECTED : 0;
		state.vector = &in;
		int err = dev.proc(&dev, &state);
		CHECK(!err, "fsp tone proc failed at %.1f Hz", freq);
		if (err) break;
		for (int a = 0; a < 2; a++) {
			double u = state.vector->data[a*state.vector->stride];
			CHECK(isfinite(u) && fabs(u) <= COMMAND_LIMIT + 1e-6,
				"bad %.1f Hz command %g", freq, u);
			delayed[a][k % DELAY] = u;
			if (k >= train) e2 += truth[a]*truth[a];
		}
	}
	free_fsp(&dev);
	return sqrt(e2/(2.0*score))*PIXEL_SCALE;
}

static int run_frequency_sweep(void)
{
	static const double freqs[] = {2, 5, 10, 20, 30, 50, 75, 100,
		150, 200, 250, 260, 275, 300, 320, 330, 340, 350,
		400, 500, 600, 750, 900, 1000};
	puts("freq_hz,no_dark_px,hold20_px,predict20_px,no_dark_db,hold20_db,predict20_db");
	for (size_t i = 0; i < sizeof freqs/sizeof freqs[0]; i++) {
		double f = freqs[i];
		double live = run_tone(true, 0, f);
		double hold = run_tone(false, DUTY_DARK, f);
		double pred = run_tone(true, DUTY_DARK, f);
		printf("%.0f,%.6f,%.6f,%.6f,%.3f,%.3f,%.3f\n", f,
			live, hold, pred, 20*log10(3.0/live),
			20*log10(3.0/hold), 20*log10(3.0/pred));
	}
	return n_fail ? 1 : 0;
}

static void free_replay(struct replay_data *r)
{
	free(r->v[0]);
	free(r->v[1]);
	free(r->dark);
}

static int load_replay(const char *path, struct replay_data *r)
{
	memset(r, 0, sizeof *r);
	FILE *fp = fopen(path, "rb");
	if (!fp) return -1;
	size_t cap = 0;
	struct aylp_header h;
	while (r->n < 10*2310 && fread(&h, sizeof h, 1, fp) == 1) {
		if (h.magic != AYLP_MAGIC) break;
		size_t n = h.log_dim.y*h.log_dim.x;
		size_t elem = h.type == AYLP_T_MATRIX_UCHAR
			|| h.type == AYLP_T_BLOCK_UCHAR ? 1 : sizeof(double);
		if (h.type != AYLP_T_VECTOR || n != 2) {
			if (fseek(fp, (long)(n*elem), SEEK_CUR)) break;
			continue;
		}
		double v[2];
		if (fread(v, sizeof v, 1, fp) != 1) break;
		if (r->n == cap) {
			cap = cap ? 2*cap : 4096;
			for (int a = 0; a < 2; a++) {
				double *p = realloc(r->v[a], cap*sizeof *p);
				if (!p) { fclose(fp); free_replay(r); return -1; }
				r->v[a] = p;
			}
			unsigned char *p = realloc(r->dark, cap);
			if (!p) { fclose(fp); free_replay(r); return -1; }
			r->dark = p;
		}
		for (int a = 0; a < 2; a++) r->v[a][r->n] = v[a];
		r->dark[r->n] = !!(h.status & (AYLP_FRAME_REJECTED
			| AYLP_NO_SIGNAL | AYLP_BEAM_LOST));
		r->n++;
	}
	fclose(fp);
	if (r->n < 8*2310) { free_replay(r); return -1; }

	/* The file contains held coordinates in rejected frames.  Estimate the
	 * unobserved physical path only for scoring, using the two surrounding
	 * accepted samples; the controller still receives a hold and the flag. */
	for (int a = 0; a < 2; a++) {
		for (size_t i = 0; i < r->n;) {
			if (!r->dark[i]) { i++; continue; }
			size_t begin = i;
			while (i < r->n && r->dark[i]) i++;
			size_t end = i;
			double left = begin ? r->v[a][begin - 1]
				: (end < r->n ? r->v[a][end] : 0.0);
			double right = end < r->n ? r->v[a][end] : left;
			for (size_t j = begin; j < end; j++)
				r->v[a][j] = left + (right-left)
					* (double)(j-begin+1)/(double)(end-begin+1);
		}
		double mean = 0.0;
		for (size_t i = 0; i < r->n; i++) mean += r->v[a][i];
		mean /= r->n;
		for (size_t i = 0; i < r->n; i++) r->v[a][i] -= mean;
	}
	return 0;
}

static struct result replay_series(const struct replay_data *r, bool predict)
{
	struct result out = {0};
	struct aylp_device dev;
	CHECK(!make_fsp_config(&dev, predict, true),
		"fsp init failed in AYLP replay");
	if (!dev.device_data) return out;
	const size_t score_at = 7*2310;
	double delayed[2][DELAY] = {{0}}, last[2] = {r->v[0][0], r->v[1][0]};
	double e2[2] = {0}, open2[2] = {0};
	for (size_t k = 0; k < r->n; k++) {
		double truth[2], in_data[2];
		for (int a = 0; a < 2; a++) {
			truth[a] = r->v[a][k] + delayed[a][k % DELAY];
			if (!r->dark[k]) last[a] = truth[a];
			in_data[a] = last[a];
		}
		gsl_vector in = {.size = 2, .stride = 1, .data = in_data,
			.block = NULL, .owner = 0};
		struct aylp_state state = {0};
		state.header.type = AYLP_T_VECTOR;
		state.header.units = AYLP_U_MINMAX;
		state.header.status = r->dark[k] ? AYLP_FRAME_REJECTED : 0;
		state.vector = &in;
		int err = dev.proc(&dev, &state);
		CHECK(!err, "fsp AYLP replay failed at frame %zu", k);
		if (err) break;
		for (int a = 0; a < 2; a++) {
			double u = state.vector->data[a*state.vector->stride];
			delayed[a][k % DELAY] = u;
			if (k >= score_at) {
				e2[a] += truth[a]*truth[a];
				open2[a] += r->v[a][k]*r->v[a][k];
			}
		}
		if (k >= score_at) out.n_dark += r->dark[k];
	}
	size_t nscore = r->n - score_at;
	out.live_rms = sqrt(e2[0]/nscore)*PIXEL_SCALE;
	out.dark_rms = sqrt(e2[1]/nscore)*PIXEL_SCALE;
	out.duty_rms = sqrt((e2[0]+e2[1])/(2.0*nscore))*PIXEL_SCALE;
	out.random_rms = sqrt(open2[0]/nscore)*PIXEL_SCALE;
	out.disturbance_rms = sqrt(open2[1]/nscore)*PIXEL_SCALE;
	free_fsp(&dev);
	return out;
}

static int run_aylp_replay(const char *path)
{
	struct replay_data r;
	if (load_replay(path, &r)) {
		fprintf(stderr, "could not load a 10 s vector replay from %s\n", path);
		return 1;
	}
	struct result hold = replay_series(&r, false);
	struct result pred = replay_series(&r, true);
	size_t nscore = r.n - 7*2310;
	printf("source=%s frames=%zu score_frames=%zu rejected=%zu (%.1f%%)\n",
		path, r.n, nscore, pred.n_dark, 100.0*pred.n_dark/nscore);
	printf("axis,open_px,hold_px,predict_px,improvement_pct\n");
	printf("y,%.6f,%.6f,%.6f,%.2f\n", hold.random_rms,
		hold.live_rms, pred.live_rms,
		100.0*(1.0-pred.live_rms/hold.live_rms));
	printf("x,%.6f,%.6f,%.6f,%.2f\n", hold.disturbance_rms,
		hold.dark_rms, pred.dark_rms,
		100.0*(1.0-pred.dark_rms/hold.dark_rms));
	printf("combined,NA,%.6f,%.6f,%.2f\n", hold.duty_rms,
		pred.duty_rms, 100.0*(1.0-pred.duty_rms/hold.duty_rms));
	free_replay(&r);
	return n_fail ? 1 : 0;
}

/* A ramp is the distinguishing case for the new second slow state: the
 * legacy EWMA must follow it with a finite lag, while the position/rate model
 * has zero steady-state following error.  Exercise the real Smith command
 * ring as well as the observer update. */
static double run_drift_ramp(unsigned order)
{
	struct aylp_device dev = {0};
	char json[1024];
	snprintf(json, sizeof json,
		"{\"type\":\"vector\",\"units\":\"V\",\"fs\":1000,"
		"\"delay\":3,\"delay_frac\":0.25,\"drift_tau\":0.5,"
		"\"drift_order\":%u,\"broad_order\":0,\"start_delay\":0,"
		"\"ramp\":0,\"guard_ratio\":0,\"adapt_period\":0,"
		"\"clamp\":10,\"y\":{\"K\":1,\"r\":1,\"freqs\":[20],"
		"\"zeta\":[0.1],\"q\":[1e-12]},\"x\":{\"K\":1,\"r\":1,"
		"\"freqs\":[20],\"zeta\":[0.1],\"q\":[1e-12]}}", order);
	dev.uri = "anyloop:fsp";
	dev.params = json_tokener_parse(json);
	CHECK(dev.params && !fsp_init(&dev), "drift-order %u init failed", order);
	if (!dev.device_data) return INFINITY;
	double delayed[2][4] = {{0}}, sum = 0.0;
	const size_t total = 30000, score = 10000;
	for (size_t k = 0; k < total; k++) {
		double phi = 0.02 * (double)k / 1000.0;
		double in_data[2] = {phi + delayed[0][k % 4],
			phi + delayed[1][k % 4]};
		gsl_vector in = {.size=2, .stride=1, .data=in_data};
		struct aylp_state state = {0};
		state.header.type = AYLP_T_VECTOR;
		state.header.units = AYLP_U_MINMAX;
		state.vector = &in;
		CHECK(!dev.proc(&dev, &state), "drift-order %u proc failed", order);
		for (int a = 0; a < 2; a++)
			delayed[a][k % 4] = state.vector->data[a];
		if (k >= total-score) sum += fabs(in_data[0]);
	}
	double mae = sum / score;
	free_fsp(&dev);
	return mae;
}

int main(int argc, char **argv)
{
	if (argc == 2 && !strcmp(argv[1], "--sweep"))
		return run_frequency_sweep();
	if (argc == 3 && !strcmp(argv[1], "--replay-aylp"))
		return run_aylp_replay(argv[2]);
	struct result live_hold = run_case(false, 0.0, 0, 3.0, 4101);
	struct result live_pred = run_case(true, 0.0, 0, 3.0, 4101);
	struct result hold = run_case(false, 0.0, DUTY_DARK, 3.0, 4101);
	struct result pred = run_case(true, 0.0, DUTY_DARK, 3.0, 4101);
	struct result hold1 = run_case(false, 0.01, DUTY_DARK, 3.0, 4101);
	struct result rand1 = run_case(true, 0.01, DUTY_DARK, 3.0, 4101);
	struct result hold5 = run_case(false, 0.05, DUTY_DARK, 3.0, 4101);
	struct result rand5 = run_case(true, 0.05, DUTY_DARK, 3.0, 4101);
	struct result hold10 = run_case(false, 0.10, DUTY_DARK, 3.0, 4101);
	struct result rand10 = run_case(true, 0.10, DUTY_DARK, 3.0, 4101);
	struct result stress_hold = run_case(false, 0.10, DUTY_DARK,
		6.0, 4101);
	struct result stress_pred = run_case(true, 0.10, DUTY_DARK,
		6.0, 4101);
	double drift_one = run_drift_ramp(1);
	double drift_two = run_drift_ramp(2);
	print_result("hold, 0% dark", live_hold);
	print_result("predict, 0% dark", live_pred);
	print_result("hold, duty only", hold);
	print_result("predict, duty only", pred);
	print_result("hold + 1% random", hold1);
	print_result("predict + 1% random", rand1);
	print_result("hold + 5% random", hold5);
	print_result("predict + 5% random", rand5);
	print_result("hold + 10% random", hold10);
	print_result("predict + 10% random", rand10);
	print_result("hold 6px + 10%", stress_hold);
	print_result("predict 6px + 10%", stress_pred);
	printf("drift ramp: order-1 MAE %.6g, order-2 MAE %.6g\n",
		drift_one, drift_two);

	CHECK(pred.rms < hold.rms,
		"dark prediction regressed duty-chop RMS: %g vs hold %g",
		pred.rms, hold.rms);
	CHECK(pred.n_predicted > 0, "dark bank never produced a command");
	CHECK(rand1.n_random_dark > 0 && rand5.n_random_dark > rand1.n_random_dark
		&& rand10.n_random_dark > rand5.n_random_dark,
		"random blank injection did not increase monotonically");
	CHECK(rand1.rms < hold1.rms && rand5.rms < hold5.rms
		&& rand10.rms < hold10.rms,
		"prediction failed to improve a matched random-blank case");
	CHECK(stress_pred.rms < stress_hold.rms,
		"prediction failed the 6 px RMS stress case");
	CHECK(fabs(live_pred.rms - live_hold.rms)*PIXEL_SCALE < 1e-3,
		"dark prediction materially changed an all-live run: %.9g px",
		fabs(live_pred.rms - live_hold.rms)*PIXEL_SCALE);
	CHECK(drift_two < 0.25 * drift_one,
		"position/rate state did not remove ramp lag: %g vs %g",
		drift_two, drift_one);

	if (n_fail) {
		fprintf(stderr, "%u fsp dark integration check(s) failed\n", n_fail);
		return 1;
	}
	puts("all fsp dark integration checks passed");
	return 0;
}
