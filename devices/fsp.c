#include <time.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <sys/mman.h>
#include <gsl/gsl_vector.h>

#include "anyloop.h"
#include "logging.h"
#include "timing.h"
#include "xalloc.h"
#include "fsp.h"


// ---- small dense linear algebra on dim <= AYLP_FSP_MAX_DIM ------------------
// Everything here runs at init and once per adaptation tick, never in the hot
// per-sample path (which uses the O(n_modes) AR recursion directly), so plain
// triple loops on a <=16x16 matrix are fine.

// build the block-diagonal AR(2) state matrix A (dim x dim) from the per-mode
// coefficients: each 2x2 block is [[a1_i, a2_i], [1, 0]]
static void fsp_build_A(const struct aylp_fsp_axis *ax, double *A)
{
	size_t D = ax->dim;
	memset(A, 0, D * D * sizeof(double));
	for (size_t i = 0; i < ax->n_modes; i++) {
		size_t b = 2 * i;
		A[(b + 0) * D + (b + 0)] = ax->a1[i];
		A[(b + 0) * D + (b + 1)] = ax->a2[i];
		A[(b + 1) * D + (b + 0)] = 1.0;
	}
}

// C = A * B (all D x D, row-major)
static void fsp_mm(const double *A, const double *B, double *C, size_t D)
{
	for (size_t i = 0; i < D; i++) for (size_t j = 0; j < D; j++) {
		double s = 0.0;
		for (size_t k = 0; k < D; k++) s += A[i*D+k] * B[k*D+j];
		C[i*D+j] = s;
	}
}

// C = A * B^T (all D x D, row-major)
static void fsp_mmT(const double *A, const double *B, double *C, size_t D)
{
	for (size_t i = 0; i < D; i++) for (size_t j = 0; j < D; j++) {
		double s = 0.0;
		for (size_t k = 0; k < D; k++) s += A[i*D+k] * B[j*D+k];
		C[i*D+j] = s;
	}
}

// Solve the stationary Kalman gain L for one axis by iterating the predicted-
// covariance Riccati recursion to convergence. The measurement is the sum of
// the modal positions (C picks out the even-indexed "x_i(k)" states), plus
// white noise of variance r; the process noise Q is diagonal with q_i on each
// position state. Returns 0 on success.
// Riccati gain solve, refactored into resumable pieces (begin / iterate /
// finalize) so the adaptation can spread it a few iterations per frame while
// init still runs it synchronously. The solve state (A, P, iteration count,
// last trace) lives on the axis so it survives between frames.
#define AYLP_FSP_RICCATI_MAXIT 2000
// Riccati iterations run per frame while a re-identification solve is in
// progress. ~16 D<=16 iterations is a few microseconds -- far below the frame
// budget -- so the solve (typically a few hundred iterations) finishes in tens
// of milliseconds without any single frame stalling.
#define AYLP_FSP_ADAPT_ITERS 16

// snapshot A and seed P = Q, then reset the iteration counter.
static int fsp_solve_gain_begin_scaled(struct aylp_fsp_axis *ax,
	double q_scale)
{
	size_t D = ax->dim;
	if (D == 0 || D > AYLP_FSP_MAX_DIM
			|| !isfinite(q_scale) || q_scale <= 0.0) return -1;
	ax->adapt_q_scale = q_scale;
	fsp_build_A(ax, ax->adapt_A);
	// Q (diagonal, q_i on position states) doubles as the initial P
	memset(ax->adapt_P, 0, D * D * sizeof(double));
	for (size_t i = 0; i < ax->n_modes; i++)
		ax->adapt_P[(2*i) * D + (2*i)] =
			(ax->q[i] > 0.0 ? ax->q[i] : 1e-9) * q_scale;
	ax->adapt_it = 0;
	ax->adapt_prev_trace = 0.0;
	return 0;
}


static int fsp_solve_gain_begin(struct aylp_fsp_axis *ax)
{
	return fsp_solve_gain_begin_scaled(ax, 1.0);
}

// advance the fixed-point iteration by up to `niter` steps.
// returns 1 = converged (or hit the iteration cap), 0 = more work to do,
// -1 = numeric failure. Does not touch ax->L (see finalize).
static int fsp_solve_gain_iterate(struct aylp_fsp_axis *ax, size_t niter)
{
	size_t D = ax->dim;
	double *A = ax->adapt_A;
	double *P = ax->adapt_P;
	double Pf[AYLP_FSP_MAX_DIM * AYLP_FSP_MAX_DIM];
	double AP[AYLP_FSP_MAX_DIM * AYLP_FSP_MAX_DIM];
	double Pn[AYLP_FSP_MAX_DIM * AYLP_FSP_MAX_DIM];
	double PCt[AYLP_FSP_MAX_DIM];		// P C^T
	double CP[AYLP_FSP_MAX_DIM];		// C P
	double Kf[AYLP_FSP_MAX_DIM];		// filter gain wrt predicted state
	for (size_t step = 0; step < niter; step++) {
		if (ax->adapt_it >= AYLP_FSP_RICCATI_MAXIT) return 1;
		// S = C P C^T + r ; PCt[r] = sum_j P[r][2j] ; CP[c] = sum_i P[2i][c]
		double S = ax->r > 0.0 ? ax->r : 1e-9;
		for (size_t rI = 0; rI < D; rI++) {
			double pc = 0.0, cp = 0.0;
			for (size_t i = 0; i < ax->n_modes; i++) {
				pc += P[rI*D + 2*i];
				cp += P[(2*i)*D + rI];
			}
			PCt[rI] = pc;
			CP[rI] = cp;
		}
		for (size_t i = 0; i < ax->n_modes; i++) S += PCt[2*i];
		if (!isfinite(S) || S <= 0.0) return -1;
		for (size_t rI = 0; rI < D; rI++) Kf[rI] = PCt[rI] / S;
		// Pf = P - Kf * CP   (filtered covariance)
		for (size_t i = 0; i < D; i++) for (size_t j = 0; j < D; j++)
			Pf[i*D+j] = P[i*D+j] - Kf[i] * CP[j];
		// Pn = A Pf A^T + Q
		fsp_mm(A, Pf, AP, D);
		fsp_mmT(AP, A, Pn, D);
		for (size_t i = 0; i < ax->n_modes; i++)
			Pn[(2*i)*D + (2*i)] +=
				(ax->q[i] > 0.0 ? ax->q[i] : 1e-9)
				* ax->adapt_q_scale;
		// convergence on the trace
		double tr = 0.0;
		for (size_t i = 0; i < D; i++) tr += Pn[i*D+i];
		memcpy(P, Pn, D * D * sizeof(double));
		size_t it = ax->adapt_it++;
		if (it > 3 && fabs(tr - ax->adapt_prev_trace)
				<= 1e-12 * (1.0 + fabs(tr)))
			return 1;
		ax->adapt_prev_trace = tr;
	}
	return 0;
}

// final gain L = P C^T / (C P C^T + r) from the converged covariance, plus
// each mode's posterior (filtered) position-state variance post_var_i =
// P_ii - L_i * PCt_i (P is symmetric, so C P's column i equals PCt[i]).
// post_var debiases the adaptation's state-energy estimator: xhat[2i] alone
// systematically understates E[x_i^2] by exactly this posterior variance.
static int fsp_solve_gain_finalize(struct aylp_fsp_axis *ax)
{
	size_t D = ax->dim;
	double *P = ax->adapt_P;
	double PCt[AYLP_FSP_MAX_DIM];
	double S = ax->r > 0.0 ? ax->r : 1e-9;
	for (size_t rI = 0; rI < D; rI++) {
		double pc = 0.0;
		for (size_t i = 0; i < ax->n_modes; i++) pc += P[rI*D + 2*i];
		PCt[rI] = pc;
	}
	for (size_t i = 0; i < ax->n_modes; i++) S += PCt[2*i];
	if (!isfinite(S) || S <= 0.0) return -1;
	for (size_t rI = 0; rI < D; rI++) ax->L[rI] = PCt[rI] / S;
	for (size_t i = 0; i < ax->n_modes; i++) {
		double pv = P[(2*i)*D + (2*i)] - ax->L[2*i] * PCt[2*i];
		ax->post_var[i] = pv > 0.0 ? pv : 0.0;
	}
	return 0;
}

// synchronous solve for the init path: run the whole iteration in one call.
static int fsp_solve_gain(struct aylp_fsp_axis *ax)
{
	if (fsp_solve_gain_begin(ax)) return -1;
	int done;
	while ((done = fsp_solve_gain_iterate(ax,
			AYLP_FSP_RICCATI_MAXIT)) == 0)
		;
	if (done < 0) return -1;
	return fsp_solve_gain_finalize(ax);
}

// Build the event-only modal observer without disturbing the normal gain or
// its adaptation workspace.  A physical push is an abrupt increase in modal
// drive, so multiplying Q makes the observer reacquire the ring-down's
// amplitude and phase quickly.  This synchronous solve is init-only; it must
// never migrate into the real-time path.
static int fsp_solve_transient_gain(struct aylp_fsp_axis *ax, double q_scale)
{
	struct aylp_fsp_axis tmp = *ax;
	for (size_t i = 0; i < tmp.n_modes; i++) tmp.q[i] *= q_scale;
	if (fsp_solve_gain(&tmp)) return -1;
	memcpy(ax->transient_L, tmp.L, ax->dim * sizeof(double));
	return 0;
}

// (re)compute the AR(2) coefficients of every mode from its f/zeta at rate fs
static void fsp_build_modes(struct aylp_fsp_axis *ax, double fs)
{
	double Ts = 1.0 / fs;
	for (size_t i = 0; i < ax->n_modes; i++) {
		double z = ax->zeta[i];
		if (z < 0.0) z = 0.0;
		if (z > 0.999) z = 0.999;
		double w = 2.0 * M_PI * ax->f[i] * Ts;
		double rpole = exp(-z * w);
		ax->a1[i] = 2.0 * rpole * cos(w * sqrt(1.0 - z*z));
		ax->a2[i] = -rpole * rpole;
		// stationary variance gain (see fsp.h); guard the (1 -/+ a1
		// - a2) factors, which approach 0 as the pole nears z = 1
		double d1 = 1.0 - ax->a1[i] - ax->a2[i];
		double d2 = 1.0 + ax->a1[i] - ax->a2[i];
		double d3 = 1.0 + ax->a2[i];
		if (d1 < 1e-12) d1 = 1e-12;
		if (d2 < 1e-12) d2 = 1e-12;
		if (d3 < 1e-12) d3 = 1e-12;
		ax->Gv[i] = (1.0 - ax->a2[i]) / (d3 * d1 * d2);
	}
}

// design one axis's burst-guard band-pass (RBJ constant-peak band-pass) at
// the parasitic-loop regeneration frequency fs/(2*(delay+delay_frac)), the
// frequency where any plant-model mismatch rings (~310-340 Hz on this bench).
// Q = 1.5 spans the measured 250-450 Hz hump. Also derives the sustained-
// energy debounce (guard_min_cycles periods of this axis's f0) below which a
// momentary envelope spike -- e.g. low-frequency pointing content leaking
// through the band-pass skirts, a false trip observed in practice -- cannot
// latch a trip the way a genuinely sustained regeneration ring does.
static void fsp_build_guard(struct aylp_fsp_axis *ax,
	const struct aylp_fsp_data *data)
{
	double fs = data->fs;
	double f0 = fs / (2.0 * ((double)ax->delay + ax->delay_frac));
	if (f0 > 0.45 * fs) f0 = 0.45 * fs;
	ax->gd_f0 = f0;
	double w0 = 2.0 * M_PI * f0 / fs;
	double alpha = sin(w0) / (2.0 * 1.5);
	double a0 = 1.0 + alpha;
	ax->gd_b0 = alpha / a0;		// b1 = 0
	ax->gd_b2 = -alpha / a0;
	ax->gd_a1 = -2.0 * cos(w0) / a0;
	ax->gd_a2 = (1.0 - alpha) / a0;
	size_t min_samples = 0;
	if (data->guard_min_cycles > 0.0 && f0 > 0.0)
		min_samples = (size_t)(data->guard_min_cycles * fs / f0 + 0.5);
	ax->gd_min_samples = min_samples;
}

// design the shared command-filter biquad (RBJ cookbook low-pass, Q=1/sqrt2)
static void fsp_build_cmdlp(struct aylp_fsp_data *data)
{
	double w0 = 2.0 * M_PI * data->cmd_fc / data->fs;
	// Q = 1/sqrt2 (Butterworth): alpha = sin(w0)/(2Q) = sin(w0)/sqrt2
	double alpha = sin(w0) * M_SQRT1_2;
	double cw = cos(w0);
	double a0 = 1.0 + alpha;
	data->lp_b0 = (1.0 - cw) / 2.0 / a0;
	data->lp_b1 = (1.0 - cw) / a0;
	data->lp_b2 = (1.0 - cw) / 2.0 / a0;
	data->lp_a1 = -2.0 * cw / a0;
	data->lp_a2 = (1.0 - alpha) / a0;
}

// (re)compute one axis's per-mode pre-compensation for the command filter:
// evaluate H(e^{j2pi f_i/fs}) at each mode center to get the real-valued
// extra prediction horizon n_i = -arg H / omega_i that cancels the filter's
// phase lag at that line, plus the capped gain boost min(1/|H|, 3) that
// cancels its droop. Called at init and on every adaptation tick (f_i wander
// moves the compensation point). See fsp.h for why this is a roll-forward
// and NOT a quadrature rotation.
//
// The mode's full real-valued forward horizon is delay + delay_frac + n_i;
// comp_n[i]/comp_frac[i] store its floor/fractional part so the per-sample
// prediction loop can blend the same way the plain delay_frac case already
// did, rather than rounding n_i to the nearest integer step and discarding
// a fractional phase advance whenever cmd_fc is enabled. With cmd_fc <= 0,
// n_i = 0 for every mode and this reduces exactly to the old delay_frac-only
// blend (comp_n[i] = delay, comp_frac[i] = delay_frac).
static void fsp_build_comp(struct aylp_fsp_axis *ax,
	const struct aylp_fsp_data *data)
{
	ax->max_steps = ax->delay + (ax->delay_frac > 0.0);
	for (size_t i = 0; i < ax->n_modes; i++) {
		double n = 0.0, g = 1.0;
		if (data->cmd_fc > 0.0) {
			double w = 2.0 * M_PI * ax->f[i] / data->fs;
			// H(z) at z = e^{jw}: (b0 + b1 z^-1 + b2 z^-2)/(1 + a1
			// z^-1 + a2 z^-2), evaluated with real/imag parts
			double c1 = cos(w), s1 = -sin(w);	// z^-1
			double c2 = cos(2*w), s2 = -sin(2*w);	// z^-2
			double nr = data->lp_b0 + data->lp_b1*c1 + data->lp_b2*c2;
			double ni = data->lp_b1*s1 + data->lp_b2*s2;
			double dr = 1.0 + data->lp_a1*c1 + data->lp_a2*c2;
			double di = data->lp_a1*s1 + data->lp_a2*s2;
			double dd = dr*dr + di*di;
			double hr = (nr*dr + ni*di) / dd;
			double hi = (ni*dr - nr*di) / dd;
			double mag = hypot(hr, hi);
			g = mag > 1e-9 ? 1.0/mag : 3.0;
			if (g > 3.0) g = 3.0;
			double th = -atan2(hi, hr);	// phase advance (rad, >= 0)
			n = th / w;
		}
		double horizon = (double)ax->delay + ax->delay_frac + n;
		size_t ti = (size_t)horizon;
		double tf = horizon - (double)ti;
		ax->comp_n[i] = ti;
		ax->comp_frac[i] = tf;
		ax->comp_g[i] = g;
		size_t steps = ti + (tf > 0.0);
		if (steps > ax->max_steps) ax->max_steps = steps;
	}
}


// ---- parameter parsing helpers ---------------------------------------------

// pull a JSON array of doubles into dst (up to max); returns count
static size_t fsp_get_darray(struct json_object *val, double *dst, size_t max)
{
	if (!json_object_is_type(val, json_type_array)) return 0;
	size_t n = json_object_array_length(val);
	if (n > max) n = max;
	for (size_t i = 0; i < n; i++)
		dst[i] = json_object_get_double(
			json_object_array_get_idx(val, i)
		);
	return n;
}

// Load a deterministic offline Wiener solution. The text format is one row
// per tap: index y_w y_w_next x_w x_w_next. Blank lines and # comments are
// ignored. Exact contiguous indices prevent silent order/file mismatches.
static int fsp_load_wiener(struct aylp_fsp_data *data)
{
	if (!data->wiener_file) return 0;
	if (!data->broad_order) {
		log_error("fsp: wiener_file requires broad_order > 0.");
		return -1;
	}
	FILE *fp = fopen(data->wiener_file, "r");
	if (!fp) {
		log_error("fsp: could not open Wiener weights \"%s\": %s",
			data->wiener_file, strerror(errno));
		return -1;
	}
	char line[4096];
	size_t row = 0, lineno = 0;
	int ret = 0;
	while (fgets(line, sizeof line, fp)) {
		lineno++;
		char *p = line;
		while (*p == ' ' || *p == '\t') p++;
		if (!*p || *p == '\n' || *p == '#') continue;
		size_t idx;
		double yw, ywn, xw, xwn;
		char extra;
		int n = sscanf(p, "%zu %lf %lf %lf %lf %c",
			&idx, &yw, &ywn, &xw, &xwn, &extra);
		if (n != 5 || idx != row || row >= data->broad_order
				|| !isfinite(yw) || !isfinite(ywn)
				|| !isfinite(xw) || !isfinite(xwn)) {
			log_error("fsp: invalid Wiener row at %s:%zu (expected "
				"contiguous index %zu and four finite weights)",
				data->wiener_file, lineno, row);
			ret = -1;
			break;
		}
		data->axis[0].broad_w[row] = yw;
		data->axis[0].broad_w_next[row] = ywn;
		data->axis[1].broad_w[row] = xw;
		data->axis[1].broad_w_next[row] = xwn;
		row++;
	}
	if (!ret && ferror(fp)) {
		log_error("fsp: error reading Wiener weights \"%s\"",
			data->wiener_file);
		ret = -1;
	}
	fclose(fp);
	if (!ret && row != data->broad_order) {
		log_error("fsp: Wiener file \"%s\" has %zu rows; broad_order "
			"requires exactly %zu", data->wiener_file, row,
			data->broad_order);
		ret = -1;
	}
	if (!ret)
		log_info("fsp: initialized y/x %zu-tap predictors from offline "
			"Wiener solution %s", row, data->wiener_file);
	return ret;
}


// Dump the learned full-band predictor in exactly the format
// fsp_load_wiener() reads, so a converged run can be analysed offline and
// replayed as the next run's wiener_file. The commented header records the
// model the weights were learned UNDER -- taps are only meaningful against
// their own order/horizon/broad_lp/K, and a file loaded under a different
// horizon is silently wrong, so the loader's contiguous-index check is not
// enough on its own. Called from fini, which runs on both the AYLP_DONE and
// the SIGINT path, so an aborted run still saves what it had learned.
static int fsp_save_wiener(struct aylp_fsp_data *data)
{
	if (!data->wiener_out) return 0;
	if (!data->broad_order) {
		log_warn("fsp: wiener_out ignored (broad_order = 0, so there "
			"is no full-band predictor to save)");
		return 0;
	}
	FILE *fp = fopen(data->wiener_out, "w");
	if (!fp) {
		log_error("fsp: could not open \"%s\" to save weights: %s",
			data->wiener_out, strerror(errno));
		return -1;
	}
	// Tap-weight norms are the cheap health check: a converged observer
	// settles, while a diverging one grows without bound.
	double n2[2] = {0.0, 0.0}, n2n[2] = {0.0, 0.0};
	for (int a = 0; a < 2; a++) {
		for (size_t i = 0; i < data->broad_order; i++) {
			n2[a] += data->axis[a].broad_w[i]
				* data->axis[a].broad_w[i];
			n2n[a] += data->axis[a].broad_w_next[i]
				* data->axis[a].broad_w_next[i];
		}
		n2[a] = sqrt(n2[a]);
		n2n[a] = sqrt(n2n[a]);
	}
	time_t now = time(NULL);
	struct tm tm_buf;
	char stamp[32] = "unknown";
	if (localtime_r(&now, &tm_buf))
		strftime(stamp, sizeof stamp, "%Y-%m-%d %H:%M:%S", &tm_buf);
	fprintf(fp, "# fsp full-band predictor weights saved %s\n", stamp);
	fprintf(fp, "# broad_order=%zu broad_mu=%G broad_lp=%zu broad_gd=%zu "
		"fs=%G\n", data->broad_order, data->broad_mu, data->broad_lp,
		data->broad_gd, data->fs);
	fprintf(fp, "# horizon y %zu+%.3G / x %zu+%.3G frames; K y %G / x %G\n",
		data->axis[0].delay, data->axis[0].delay_frac,
		data->axis[1].delay, data->axis[1].delay_frac,
		data->axis[0].K, data->axis[1].K);
	fprintf(fp, "# frames closed=%zu; ||w|| y %.6G / x %.6G; "
		"||w_next|| y %.6G / x %.6G\n", data->n_closed,
		n2[0], n2[1], n2n[0], n2n[1]);
	fprintf(fp, "# freeze_closed=%s init=%s\n",
		data->broad_freeze_closed ? "true" : "false",
		data->wiener_file ? data->wiener_file : "zeros (cold start)");
	fprintf(fp, "# reload with \"wiener_file\", but ONLY into a run with "
		"the same broad_order, horizon, broad_lp and drift_tau\n");
	fprintf(fp, "# index y_w y_w_next x_w x_w_next\n");
	for (size_t i = 0; i < data->broad_order; i++) {
		fprintf(fp, "%zu %.17G %.17G %.17G %.17G\n", i,
			data->axis[0].broad_w[i], data->axis[0].broad_w_next[i],
			data->axis[1].broad_w[i],
			data->axis[1].broad_w_next[i]);
	}
	int ret = 0;
	if (ferror(fp)) {
		log_error("fsp: error writing weights to \"%s\"",
			data->wiener_out);
		ret = -1;
	}
	if (fclose(fp)) {
		log_error("fsp: error closing \"%s\": %s", data->wiener_out,
			strerror(errno));
		ret = -1;
	}
	if (!ret)
		log_info("fsp: saved y/x %zu-tap predictors to %s (||w|| y "
			"%.4G / x %.4G after %zu closed frames)",
			data->broad_order, data->wiener_out, n2[0], n2[1],
			data->n_closed);
	return ret;
}


// Append one convergence sample: where the taps are, and how far they moved
// since the last sample. ||dw|| -> 0 is the observer having converged; a
// settle_time is only defensible if this has flattened before the scored
// window opens. Deliberately a short line rather than a full tap snapshot --
// this runs in the RT path, so the write must stay small.
static void fsp_trace_sample(struct aylp_fsp_data *data, double now)
{
	double wn[2], wnn[2], dwn[2];
	for (int a = 0; a < 2; a++) {
		double s = 0.0, sn = 0.0, sd = 0.0;
		for (size_t i = 0; i < data->broad_order; i++) {
			double w = data->axis[a].broad_w[i];
			double d = w - data->wtrace_prev[a][i];
			s += w * w;
			sn += data->axis[a].broad_w_next[i]
				* data->axis[a].broad_w_next[i];
			sd += d * d;
			data->wtrace_prev[a][i] = w;
		}
		wn[a] = sqrt(s);
		wnn[a] = sqrt(sn);
		dwn[a] = sqrt(sd);
	}
	fprintf(data->wiener_trace_fp,
		"%.3f %zu %.8G %.8G %.8G %.8G %.8G %.8G %.6G %.6G "
		"%zu %zu %zu %zu %.5G\n",
		now - data->t0, data->n_closed,
		wn[0], wnn[0], dwn[0], wn[1], wnn[1], dwn[1],
		data->axis[0].drift_hat, data->axis[1].drift_hat,
		data->axis[0].guard_events, data->axis[1].guard_events,
		data->axis[0].transient_events,
		data->axis[1].transient_events, data->broad_mu_cur);
	// flushed every sample on purpose: the whole point is that the trace
	// survives an abort, and it is one short line every few seconds
	fflush(data->wiener_trace_fp);
}


static const char *fsp_transient_controller_name(
	enum aylp_fsp_transient_controller controller)
{
	switch (controller) {
	case AYLP_FSP_TRANSIENT_PROPORTIONAL: return "proportional";
	case AYLP_FSP_TRANSIENT_INTEGRAL: return "integral";
	case AYLP_FSP_TRANSIENT_PI: return "pi";
	case AYLP_FSP_TRANSIENT_MODAL: return "modal";
	case AYLP_FSP_TRANSIENT_HYBRID: return "hybrid";
	}
	return "unknown";
}


static bool fsp_transient_uses_integral(
	enum aylp_fsp_transient_controller controller)
{
	return controller == AYLP_FSP_TRANSIENT_INTEGRAL
		|| controller == AYLP_FSP_TRANSIENT_PI
		|| controller == AYLP_FSP_TRANSIENT_HYBRID;
}


static bool fsp_transient_uses_modal(
	enum aylp_fsp_transient_controller controller)
{
	return controller == AYLP_FSP_TRANSIENT_MODAL
		|| controller == AYLP_FSP_TRANSIENT_HYBRID;
}


// Learn a candidate stationary predictor without touching the production FIR.
// Errors are measured before this sample's update (prequential/held-out), so a
// shadow cannot qualify merely by fitting the samples it has just seen. A
// promotion means the persistent disturbance has demonstrated a new predictable
// regime; the ordinary event detector still requires the promoted model's
// residual to remain quiet before cross-fading out of recovery.
static void fsp_shadow_update(struct aylp_fsp_axis *ax,
	struct aylp_fsp_data *data, double pe, double pe2,
	double primary_pred1, double primary_pred2,
	double energy1, double energy2, bool suppress, double now, size_t axis)
{
	if (!ax->shadow_active || suppress || ax->transient_saturated) {
		ax->shadow_advantage_start = 0.0;
		return;
	}
	if (ax->transient_recovering) {
		ax->shadow_advantage_start = 0.0;
		return;
	}
	size_t P = data->broad_order;
	double pred1 = 0.0, pred2 = 0.0;
	for (size_t i = 0; i < P; i++) {
		pred1 += ax->shadow_w[i] * ax->broad_xbuf[i];
		pred2 += ax->shadow_w_next[i] * ax->broad_xbuf2[i];
	}
	// Since pe = target - primary_prediction, target - shadow_prediction is
	// pe + primary_prediction - shadow_prediction.
	double spe = pe + primary_pred1 - pred1;
	double spe2 = pe2 + primary_pred2 - pred2;
	if (!isfinite(spe) || !isfinite(spe2)) {
		ax->shadow_active = false;
		log_warn("fsp: %s event shadow predictor became non-finite; "
			"discarding candidate", axis == 0 ? "y" : "x");
		return;
	}
	// Qualify against the exact delay-horizon innovation used by
	// fsp_transient_update() to hold/release the event. The adjacent horizon is
	// still trained and promoted as a pair, but improvement confined to that
	// auxiliary fractional-delay model cannot pass this gate.
	double primary_sample = pe * pe;
	double shadow_sample = spe * spe;
	if (!ax->shadow_frames) {
		ax->shadow_primary_e2 = primary_sample;
		ax->shadow_e2 = shadow_sample;
	} else {
		ax->shadow_primary_e2 += data->transient_shadow_beta
			* (primary_sample - ax->shadow_primary_e2);
		ax->shadow_e2 += data->transient_shadow_beta
			* (shadow_sample - ax->shadow_e2);
	}
	ax->shadow_frames++;

	// Train only the shadow, with the event threshold as the influence bound.
	double learn1 = spe, learn2 = spe2;
	double clip = ax->transient_threshold;
	if (clip > 0.0 && fabs(learn1) > clip)
		learn1 = copysign(clip, learn1);
	if (clip > 0.0 && fabs(learn2) > clip)
		learn2 = copysign(clip, learn2);
	double step1 = data->transient_shadow_mu * learn1 / energy1;
	double step2 = data->transient_shadow_mu * learn2 / energy2;
	double primary_norm2 = 0.0, delta_norm2 = 0.0;
	for (size_t i = 0; i < P; i++) {
		// same exclusion as the production weights: a fabricated tap
		// predicts but never teaches (see broad_syn)
		if (!ax->broad_synbuf[i])
			ax->shadow_w[i] += step1 * ax->broad_xbuf[i];
		if (!ax->broad_synbuf2[i])
			ax->shadow_w_next[i] += step2 * ax->broad_xbuf2[i];
		double d1 = ax->shadow_w[i] - ax->broad_w[i];
		double d2 = ax->shadow_w_next[i] - ax->broad_w_next[i];
		delta_norm2 += d1*d1 + d2*d2;
		primary_norm2 += ax->broad_w[i]*ax->broad_w[i]
			+ ax->broad_w_next[i]*ax->broad_w_next[i];
	}
	if (ax->shadow_promoted) {
		// The candidate has already passed the persistent-regime gate. Keep the
		// active model synchronized with its conservative updates so the normal
		// release detector can actually observe agreement with the new regime.
		// Event recovery still owns the command until the usual quiet hold/ramp.
		memcpy(ax->broad_w, ax->shadow_w, P * sizeof(double));
		memcpy(ax->broad_w_next, ax->shadow_w_next, P * sizeof(double));
		return;
	}

	double ratio2 = data->transient_shadow_ratio
		* data->transient_shadow_ratio;
	double norm = sqrt(primary_norm2);
	double norm_limit = data->transient_shadow_norm_ratio * fmax(norm, 0.1);
	bool eligible = now - ax->transient_start
			>= data->transient_shadow_min_duration
		&& !ax->transient_saturated
		&& ax->shadow_primary_e2 > 1e-12
		&& ax->shadow_e2 <= ratio2 * ax->shadow_primary_e2
		&& sqrt(delta_norm2) <= norm_limit;
	if (!eligible) {
		ax->shadow_advantage_start = 0.0;
		return;
	}
	if (ax->shadow_advantage_start <= 0.0) {
		ax->shadow_advantage_start = now;
		return;
	}
	if (now - ax->shadow_advantage_start < data->transient_shadow_hold)
		return;

	memcpy(ax->broad_w, ax->shadow_w, P * sizeof(double));
	memcpy(ax->broad_w_next, ax->shadow_w_next, P * sizeof(double));
	ax->shadow_error_ratio = sqrt(ax->shadow_e2 / ax->shadow_primary_e2);
	ax->shadow_promoted = true;
	ax->shadow_promotions++;
	data->shadow_promotions++;
	log_warn("fsp: %s event shadow model promoted after %.3G s: held-out "
		"RMS ratio %.3G, coefficient change %.3G (limit %.3G); event "
		"recovery remains active until the new model is quiet",
		axis == 0 ? "y" : "x", now - ax->transient_start,
		ax->shadow_error_ratio, sqrt(delta_norm2), norm_limit);
}


// Update the innovation-triggered recovery state and return the fraction of
// event recovery to use (0 = normal predictor, 1 = selected event controller).
static double fsp_transient_update(struct aylp_fsp_axis *ax,
	struct aylp_fsp_data *data, double resid, double authority,
	double now, size_t axis)
{
	if (data->transient_sigma <= 0.0 || !isfinite(resid)) return 0.0;
	// During post-close warm-up, learn the residual baseline but do not trip.
	// This is deliberately not an open-loop/frozen calibration: the EWMA keeps
	// tracking whenever no event is active, and the broadband weights continue
	// adapting in closed loop. Early cold-predictor errors are forgotten over
	// transient_tau before the detector arms.
	bool armed = now - data->t0 >=
		data->start_delay + data->transient_arm_delay;
	double sigma = sqrt(fmax(ax->transient_var, 0.0));
	double threshold = data->transient_sigma * sigma;
	if (threshold < data->transient_floor)
		threshold = data->transient_floor;
	ax->transient_threshold = threshold;
	bool outside_model = armed && fabs(resid) > threshold;
	// Authority gates entry only. Once recovery is active, a burst-guard hold
	// or another temporary authority reduction must not masquerade as model
	// agreement and release the event while innovation is still high.
	// Online modal re-identification refreshes the event-only gain after the
	// normal gain.  Do not enter a modal event during that short solve window:
	// its AR coefficients already describe the new model while transient_L
	// still describes the old one.
	enum aylp_fsp_transient_controller entry_controller =
		data->transient_controller;
	if (data->transient_modal_ab && ((ax->transient_events + 1) & 1))
		entry_controller = AYLP_FSP_TRANSIENT_PROPORTIONAL;
	bool modal_ready = !fsp_transient_uses_modal(entry_controller)
		|| ax->transient_gain_current;
	bool event_trigger = authority > 0.0 && outside_model && modal_ready;
	if (ax->transient_blocked) {
		if (outside_model) return 0.0;
		ax->transient_blocked = false;
		log_info("fsp: %s transient detector re-armed after innovation "
			"returned below threshold", axis == 0 ? "y" : "x");
	}

	if (!ax->transient_active && event_trigger) {
		ax->transient_active = true;
		ax->transient_recovering = false;
		ax->transient_t_event = now;
		ax->transient_events++;
		// entry_controller also handles backward-compatible same-run push A/B:
		// odd events proportional, even events the configured controller.
		ax->transient_controller = entry_controller;
		ax->transient_modal = fsp_transient_uses_modal(
			ax->transient_controller);
		ax->transient_i = 0.0;
		ax->transient_i_prev = 0.0;
		ax->transient_i_step = 0.0;
		ax->transient_peak_i = 0.0;
		ax->transient_start = now;
		ax->transient_error_t_last = now;
		ax->transient_peak_error = 0.0;
		ax->transient_peak_command = 0.0;
		ax->transient_error2 = 0.0;
		ax->transient_frames = 0;
		ax->shadow_active = data->transient_shadow_mu > 0.0
			&& data->broad_order && ax->shadow_w && ax->shadow_w_next;
		ax->shadow_promoted = false;
		ax->transient_saturated = false;
		ax->shadow_primary_e2 = 0.0;
		ax->shadow_e2 = 0.0;
		ax->shadow_advantage_start = 0.0;
		ax->shadow_error_ratio = NAN;
		ax->shadow_frames = 0;
		if (ax->shadow_active) {
			memcpy(ax->shadow_w, ax->broad_w,
				data->broad_order * sizeof(double));
			memcpy(ax->shadow_w_next, ax->broad_w_next,
				data->broad_order * sizeof(double));
		}
		if (fsp_transient_uses_integral(ax->transient_controller))
			log_warn("fsp: transient recovery on %s axis (event %zu): "
				"|innovation| %.4G > %.4G; using %s recovery, "
				"Kp=%G, Ki=%G/s", axis == 0 ? "y" : "x",
				ax->transient_events, fabs(resid), threshold,
				fsp_transient_controller_name(ax->transient_controller),
				data->transient_kp, data->transient_ki);
		else
			log_warn("fsp: transient recovery on %s axis (event %zu): "
				"|innovation| %.4G > %.4G; using %s recovery, Kp=%G",
				axis == 0 ? "y" : "x", ax->transient_events,
				fabs(resid), threshold,
				fsp_transient_controller_name(ax->transient_controller),
				data->transient_kp);
	} else if (ax->transient_active && outside_model) {
		ax->transient_t_event = now;
		ax->transient_recovering = false;
	} else if (!ax->transient_active && !outside_model) {
		// Once armed, an out-of-model sample is never part of the quiet baseline,
		// even if a burst-guard hold or a pending modal-gain solve temporarily
		// prevents event entry.  Otherwise the guard can teach the detector that
		// the very disturbance it should catch is ordinary noise.
		ax->transient_var += data->transient_beta
			* (resid * resid - ax->transient_var);
	}

	if (!armed && !ax->transient_active) return 0.0;
	if (!ax->transient_active) return 0.0;
	double quiet = now - ax->transient_t_event;
	if (quiet <= data->transient_hold) return 1.0;
	if (!ax->transient_recovering) {
		ax->transient_recovering = true;
		log_info("fsp: %s transient quiet; cross-fading back to the "
			"predictor over %G s", axis == 0 ? "y" : "x",
			data->transient_ramp);
	}
	double ramp = data->transient_ramp > 0.0
		? data->transient_ramp : 1e-9;
	double mix = 1.0 - (quiet - data->transient_hold) / ramp;
	if (mix <= 0.0) {
		double duration = now - ax->transient_start;
		double error_settle =
			ax->transient_error_t_last - ax->transient_start;
		double innovation_quiet =
			ax->transient_t_event - ax->transient_start;
		double rms = ax->transient_frames
			? sqrt(ax->transient_error2 / (double)ax->transient_frames)
			: 0.0;
		log_info("fsp: %s transient %zu (%s) recovered: output settle %.4G s, "
			"innovation quiet %.4G s, handoff %.4G s, peak |e| %.4G, "
			"rms(e) %.4G, "
			"peak |u| %.4G, peak |I| %.4G",
			axis == 0 ? "y" : "x", ax->transient_events,
			fsp_transient_controller_name(ax->transient_controller), error_settle,
			innovation_quiet, duration, ax->transient_peak_error, rms,
			ax->transient_peak_command, ax->transient_peak_i);
		if (data->transient_log_fp) {
			fprintf(data->transient_log_fp,
				"%s,%zu,%s,%.9G,%.9G,%.9G,%.9G,%.9G,%.9G,%.9G,%.9G,%zu,%d,%.9G\n",
				axis == 0 ? "y" : "x", ax->transient_events,
				fsp_transient_controller_name(ax->transient_controller),
				ax->transient_start - data->t0, error_settle,
				innovation_quiet, duration,
				ax->transient_peak_error, rms,
				ax->transient_peak_command, ax->transient_peak_i,
				ax->transient_frames, ax->shadow_promoted,
				ax->shadow_error_ratio);
			fflush(data->transient_log_fp);
		}
		ax->transient_active = false;
		ax->transient_recovering = false;
		ax->shadow_active = false;
		ax->transient_i = 0.0;
		ax->transient_i_step = 0.0;
		return 0.0;
	}
	return mix;
}


// Jury stability test for z^2 + c1*z + c2.  The same test is applied to
// plant_a (forward filter poles) and plant_b/b0 (inverse-filter poles), making
// the configured plant shape both stable and minimum phase.
static bool fsp_biquad_stable(double c1, double c2)
{
	return isfinite(c1) && isfinite(c2) && fabs(c2) < 1.0
		&& 1.0 + c1 + c2 > 0.0
		&& 1.0 - c1 + c2 > 0.0;
}

// Direct-form II transposed biquad used for both the delayed-command plant
// model H(z) and the matched command prefilter H^-1(z).
static inline double fsp_biquad(double x, const double b[3],
	const double a[3], double *z1, double *z2)
{
	double y = b[0]*x + *z1;
	*z1 = b[1]*x - a[1]*y + *z2;
	*z2 = b[2]*x - a[2]*y;
	return y;
}

// parse the per-axis mode/plant params out of a nested JSON object
static int fsp_parse_axis(struct aylp_fsp_axis *ax, struct json_object *obj)
{
	double q_scalar = -1.0;
	size_t nf = 0, nz = 0, nq = 0, nb = 0, na = 0;
	bool have_plant_b = false, have_plant_a = false;
	// sentinels: inherit the global delay/delay_frac unless the axis sets them
	ax->delay = 0;
	ax->delay_auto = false;
	ax->clamp_lo = NAN;		// NAN = inherit the global bound
	ax->clamp_hi = NAN;
	ax->delay_frac = -1.0;
	ax->plant_b[0] = ax->plant_a[0] = 1.0;
	json_object_object_foreach(obj, key, val) {
		if (key[0] == '_') {
		} else if (!strcmp(key, "K") || !strcmp(key, "plant_gain")) {
			ax->K = json_object_get_double(val);
		} else if (!strcmp(key, "delay")) {
			if (json_object_is_type(val, json_type_string)
					&& !strcmp(json_object_get_string(val),
						"auto"))
				ax->delay_auto = true;
			else ax->delay = json_object_get_uint64(val);
		} else if (!strcmp(key, "delay_frac")) {
			ax->delay_frac = json_object_get_double(val);
		} else if (!strcmp(key, "clamp_min")) {
			ax->clamp_lo = json_object_get_double(val);
		} else if (!strcmp(key, "clamp_max")) {
			ax->clamp_hi = json_object_get_double(val);
		} else if (!strcmp(key, "plant_b")) {
			have_plant_b = true;
			if (json_object_is_type(val, json_type_array)) {
				nb = json_object_array_length(val);
				if (nb == 3) fsp_get_darray(val, ax->plant_b, 3);
			}
		} else if (!strcmp(key, "plant_a")) {
			have_plant_a = true;
			if (json_object_is_type(val, json_type_array)) {
				na = json_object_array_length(val);
				if (na == 3) fsp_get_darray(val, ax->plant_a, 3);
			}
		} else if (!strcmp(key, "r")) {
			ax->r = json_object_get_double(val);
		} else if (!strcmp(key, "freqs")) {
			nf = fsp_get_darray(val, ax->f, AYLP_FSP_MAX_MODES);
		} else if (!strcmp(key, "zeta")) {
			nz = fsp_get_darray(val, ax->zeta, AYLP_FSP_MAX_MODES);
		} else if (!strcmp(key, "q")) {
			if (json_object_is_type(val, json_type_array))
				nq = fsp_get_darray(val, ax->q,
					AYLP_FSP_MAX_MODES);
			else
				q_scalar = json_object_get_double(val);
		} else {
			log_warn("fsp: unknown axis parameter \"%s\"", key);
		}
	}
	ax->n_modes = nf;
	ax->dim = 2 * nf;
	if (!nf) {
		log_error("fsp: each axis needs a non-empty \"freqs\" array.");
		return -1;
	}
	if (ax->dim > AYLP_FSP_MAX_DIM) {
		log_error("fsp: at most %d modes per axis.", AYLP_FSP_MAX_MODES);
		return -1;
	}
	// fill defaults: zeta 0.002 (lightly damped line), q from the scalar or
	// 1e-5 -- q is the DRIVE variance in (normalized error units)^2 per
	// sample and must be physically scaled against r; q/r near 1 puts the
	// whole loop in the waterbed-amplifying regime (see doc/devices/fsp.md)
	for (size_t i = 0; i < nf; i++) {
		if (i >= nz) ax->zeta[i] = 0.002;
		if (i >= nq) ax->q[i] = q_scalar > 0.0 ? q_scalar : 1e-5;
	}
	if (ax->r <= 0.0) ax->r = 1.0;
	if (ax->K == 0.0) {
		log_error("fsp: each axis needs a nonzero plant gain \"K\".");
		return -1;
	}
	if (have_plant_b != have_plant_a
			|| (have_plant_b && (nb != 3 || na != 3))) {
		log_error("fsp: plant_b and plant_a must both contain exactly "
			"three coefficients.");
		return -1;
	}
	if (have_plant_b) {
		for (size_t i = 0; i < 3; i++) {
			if (!isfinite(ax->plant_b[i]) || !isfinite(ax->plant_a[i])) {
				log_error("fsp: plant biquad coefficients must be finite.");
				return -1;
			}
		}
		if (fabs(ax->plant_a[0]) < 1e-12) {
			log_error("fsp: plant_a[0] must be nonzero.");
			return -1;
		}
		// Normalize the forward denominator to a0=1.
		double a0 = ax->plant_a[0];
		for (size_t i = 0; i < 3; i++) {
			ax->plant_b[i] /= a0;
			ax->plant_a[i] /= a0;
		}
		if (fabs(ax->plant_b[0]) < 1e-12
				|| !fsp_biquad_stable(ax->plant_a[1],
					ax->plant_a[2])
				|| !fsp_biquad_stable(ax->plant_b[1]/ax->plant_b[0],
					ax->plant_b[2]/ax->plant_b[0])) {
			log_error("fsp: plant biquad and its inverse must both be "
				"stable (all poles and zeros strictly inside the unit "
				"circle).");
			return -1;
		}
		double bdc = ax->plant_b[0] + ax->plant_b[1]
			+ ax->plant_b[2];
		double adc = 1.0 + ax->plant_a[1] + ax->plant_a[2];
		if (fabs(bdc) < 1e-9 || fabs(adc) < 1e-9) {
			log_error("fsp: plant biquad must have finite, nonzero DC "
				"gain.");
			return -1;
		}
		// H^-1(z) = A(z)/B(z), normalized by the leading b0.
		ax->plant_ib[0] = 1.0 / ax->plant_b[0];
		ax->plant_ib[1] = ax->plant_a[1] / ax->plant_b[0];
		ax->plant_ib[2] = ax->plant_a[2] / ax->plant_b[0];
		ax->plant_ia[0] = 1.0;
		ax->plant_ia[1] = ax->plant_b[1] / ax->plant_b[0];
		ax->plant_ia[2] = ax->plant_b[2] / ax->plant_b[0];
		ax->plant_shaped = true;
	}
	return 0;
}


static void fsp_adapt(struct aylp_fsp_data *data);
static int fsp_adapt_advance(struct aylp_fsp_axis *ax,
	struct aylp_fsp_data *data, size_t axis, size_t niter);

int fsp_init(struct aylp_device *self)
{
	self->proc = &fsp_proc;
	self->fini = &fsp_fini;
	self->device_data = xcalloc(1, sizeof(struct aylp_fsp_data));
	struct aylp_fsp_data *data = self->device_data;

	// defaults
	data->delay = 5;		// ~2 ms at ~2310 Hz, matches the loop
	data->delay_frac = 0.0;
	data->delay_auto_max = 16;
	data->delay_auto_settle = 0.5;
	data->delay_auto_window = 1.0;
	data->delay_auto_bias = 1.5;	// 1 frame of sensor + half a frame of ZOH
	data->fs = 2310.0;
	data->clamp = 1.0;
	// NAN = not explicitly configured; resolved from `clamp` after parsing
	data->clamp_lo = NAN;
	data->clamp_hi = NAN;
	data->start_delay = 0.0;
	data->precenter_delay = 0.0;
	data->precenter_ramp = 0.0;
	data->precenter_clamp = 0.0;
	data->ramp = 10.0;
	data->adapt_period = 0.0;	// fixed FSP unless asked
	data->adapt_df_max = 0.5;	// Hz per update
	data->adapt_tau = 5.0;
	data->cmd_fc = 0.0;		// raw minimum-variance command unless asked
	data->broad_order = 0;		// modal-only unless explicitly identified
	data->broad_mu = 0.03;
	data->broad_mu_init = 0.0;	// 0 = schedule disabled
	data->broad_mu_tau = 30.0;
	data->broad_lp = 0;		// raw broadband observer unless asked
	data->broad_freeze_closed = true;
	data->dark_predict = true;	// predict through a dark frame
	data->dark_predict_max = 0;	// 0 = auto, 2*(delay + broad_gd)
	data->dark_train_min_real = 0.5;
	data->dark_fab_frames = 0;	// 0 = auto, broad_hist_len/2
	data->dark_bank_train = 4;	// bank horizons trained per live frame
	data->drift_tau = 0.0;		// compound predictor unless asked
	data->drift_order = 1;		// legacy one-state EWMA unless asked
	data->transient_sigma = 0.0;	// transient fallback is opt-in
	data->transient_floor = 0.02;
	data->transient_settle_error = 0.03;
	data->transient_kp = 0.25;
	data->transient_ki = 5.0;
	data->transient_i_leak = 0.0;
	data->transient_i_limit = 1.0;
	data->transient_mode = AYLP_FSP_TRANSIENT_MANUAL;
	data->transient_controller = AYLP_FSP_TRANSIENT_PROPORTIONAL;
	data->transient_modal_q_scale = 0.0;
	data->transient_tau = 5.0;
	data->transient_hold = 0.10;
	data->transient_ramp = 0.25;
	data->transient_clamp_lo = NAN;
	data->transient_clamp_hi = NAN;
	data->transient_trip_command = NAN;
	data->transient_max_duration = 0.0;
	data->transient_arm_delay = 0.0;	// arm at close unless asked
	data->transient_shadow_mu = 0.0;	// opt-in persistent-regime recovery
	data->transient_shadow_tau = 0.5;
	data->transient_shadow_min_duration = 2.0;
	data->transient_shadow_hold = 1.0;
	data->transient_shadow_ratio = 0.70;
	data->transient_shadow_norm_ratio = 1.0;
	data->trip_frames = 8;
	data->beam_recover_ramp = 0.5;
	// burst guard defaults: on. Quiet-bench 250-450 Hz envelope is ~0.003
	// normalized; the floor keeps the bar at 4 * 0.008 = 0.032 (~0.5 px),
	// well below the multi-px bursts and well above ambient.
	data->guard_ratio = 4.0;
	data->guard_floor = 0.008;
	data->guard_hold = 0.25;
	data->guard_ramp = 1.0;
	data->guard_tick = 10.0;
	// require ~3 cycles of sustained envelope before latching a trip;
	// <= 0 reverts to the old single-sample trigger
	data->guard_min_cycles = 3.0;
	// stall-gap patch default: on. The response is proportional to the
	// frames missed (pad + continue, no authority change), so the trigger
	// can sit just above normal cadence jitter: 2 ms is ~7 frames at
	// 3788 Hz, and it also catches the ~0.3-0.6 s camera-stall hole
	// (see fsp.h).
	data->gap_trip = 0.002;

	if (!self->params) {
		log_error("fsp: no params object found.");
		return -1;
	}
	struct json_object *axy = NULL, *axx = NULL;
	json_object_object_foreach(self->params, key, val) {
		if (key[0] == '_') {
		} else if (!strcmp(key, "type")) {
			const char *s = json_object_get_string(val);
			if (!strcmp(s, "vector")) data->type = AYLP_T_VECTOR;
			else log_error("fsp: unrecognized type: %s", s);
		} else if (!strcmp(key, "units")) {
			data->units = aylp_units_from_string(
				json_object_get_string(val));
		} else if (!strcmp(key, "delay")) {
			if (json_object_is_type(val, json_type_string)
					&& !strcmp(json_object_get_string(val),
						"auto"))
				data->delay_auto = true;
			else data->delay = json_object_get_uint64(val);
		} else if (!strcmp(key, "delay_frac")) {
			data->delay_frac = json_object_get_double(val);
		} else if (!strcmp(key, "delay_auto_max")) {
			data->delay_auto_max = json_object_get_uint64(val);
		} else if (!strcmp(key, "delay_auto_settle")) {
			data->delay_auto_settle = json_object_get_double(val);
		} else if (!strcmp(key, "delay_auto_window")) {
			data->delay_auto_window = json_object_get_double(val);
		} else if (!strcmp(key, "delay_auto_bias")) {
			if (json_object_is_type(val, json_type_string)
					&& !strcmp(json_object_get_string(val),
						"auto"))
				data->delay_auto_bias_auto = true;
			else data->delay_auto_bias =
				json_object_get_double(val);
		} else if (!strcmp(key, "delay_ident_ms")) {
			if (json_object_is_type(val, json_type_string)
					&& !strcmp(json_object_get_string(val),
						"auto"))
				data->delay_ident_auto = true;
			else data->delay_ident_ms = json_object_get_double(val);
		} else if (!strcmp(key, "fs")) {
			if (json_object_is_type(val, json_type_string)
					&& !strcmp(json_object_get_string(val),
						"auto"))
				data->fs_auto = true;
			else data->fs = json_object_get_double(val);
		} else if (!strcmp(key, "clamp")) {
			data->clamp = fabs(json_object_get_double(val));
		} else if (!strcmp(key, "clamp_min")) {
			data->clamp_lo = json_object_get_double(val);
		} else if (!strcmp(key, "clamp_max")) {
			data->clamp_hi = json_object_get_double(val);
		} else if (!strcmp(key, "start_delay")) {
			data->start_delay = json_object_get_double(val);
			if (data->start_delay < 0) data->start_delay = 0;
		} else if (!strcmp(key, "precenter_delay")) {
			data->precenter_delay = json_object_get_double(val);
		} else if (!strcmp(key, "precenter_ramp")) {
			data->precenter_ramp = json_object_get_double(val);
		} else if (!strcmp(key, "precenter_clamp")) {
			data->precenter_clamp = fabs(json_object_get_double(val));
		} else if (!strcmp(key, "ramp")) {
			data->ramp = json_object_get_double(val);
			if (data->ramp < 0) data->ramp = 0;
		} else if (!strcmp(key, "adapt_period")) {
			data->adapt_period = json_object_get_double(val);
		} else if (!strcmp(key, "adapt_df_max")) {
			data->adapt_df_max = fabs(json_object_get_double(val));
		} else if (!strcmp(key, "adapt_tau")) {
			data->adapt_tau = json_object_get_double(val);
		} else if (!strcmp(key, "cmd_fc")) {
			data->cmd_fc = json_object_get_double(val);
		} else if (!strcmp(key, "broad_order")) {
			data->broad_order = json_object_get_uint64(val);
		} else if (!strcmp(key, "broad_mu")) {
			data->broad_mu = json_object_get_double(val);
		} else if (!strcmp(key, "broad_mu_init")) {
			data->broad_mu_init = json_object_get_double(val);
		} else if (!strcmp(key, "broad_mu_tau")) {
			data->broad_mu_tau = json_object_get_double(val);
		} else if (!strcmp(key, "wiener_file")) {
			xfree(data->wiener_file);
			data->wiener_file =
				xstrdup(json_object_get_string(val));
		} else if (!strcmp(key, "wiener_out")) {
			xfree(data->wiener_out);
			data->wiener_out =
				xstrdup(json_object_get_string(val));
		} else if (!strcmp(key, "wiener_trace")) {
			xfree(data->wiener_trace);
			data->wiener_trace =
				xstrdup(json_object_get_string(val));
		} else if (!strcmp(key, "wiener_trace_period")) {
			data->wiener_trace_period =
				json_object_get_double(val);
		} else if (!strcmp(key, "broad_lp")) {
			data->broad_lp = json_object_get_uint64(val);
		} else if (!strcmp(key, "dark_predict")) {
			data->dark_predict = json_object_get_boolean(val);
		} else if (!strcmp(key, "dark_predict_max")) {
			data->dark_predict_max = json_object_get_uint64(val);
		} else if (!strcmp(key, "dark_train_min_real")) {
			data->dark_train_min_real = json_object_get_double(val);
		} else if (!strcmp(key, "dark_fab_frames")) {
			data->dark_fab_frames = json_object_get_uint64(val);
		} else if (!strcmp(key, "dark_bank_train")) {
			data->dark_bank_train = json_object_get_uint64(val);
		} else if (!strcmp(key, "drift_tau")) {
			data->drift_tau = json_object_get_double(val);
		} else if (!strcmp(key, "drift_order")) {
			data->drift_order = json_object_get_uint64(val);
		} else if (!strcmp(key, "transient_mode")) {
			const char *s = json_object_get_string(val);
			if (!strcmp(s, "manual"))
				data->transient_mode = AYLP_FSP_TRANSIENT_MANUAL;
			else if (!strcmp(s, "off"))
				data->transient_mode = AYLP_FSP_TRANSIENT_OFF;
			else if (!strcmp(s, "auto"))
				data->transient_mode = AYLP_FSP_TRANSIENT_AUTO;
			else {
				log_error("fsp: transient_mode must be manual, off, or auto "
					"(got \"%s\").", s);
				return -1;
			}
		} else if (!strcmp(key, "transient_sigma")) {
			data->transient_sigma = json_object_get_double(val);
		} else if (!strcmp(key, "transient_floor")) {
			data->transient_floor =
				fabs(json_object_get_double(val));
		} else if (!strcmp(key, "transient_settle_error")) {
			data->transient_settle_error =
				fabs(json_object_get_double(val));
		} else if (!strcmp(key, "transient_kp")) {
			data->transient_kp = json_object_get_double(val);
		} else if (!strcmp(key, "transient_ki")) {
			data->transient_ki = json_object_get_double(val);
		} else if (!strcmp(key, "transient_i_leak")) {
			data->transient_i_leak = json_object_get_double(val);
		} else if (!strcmp(key, "transient_i_limit")) {
			data->transient_i_limit =
				fabs(json_object_get_double(val));
		} else if (!strcmp(key, "transient_controller")) {
			const char *s = json_object_get_string(val);
			data->transient_controller_set = true;
			if (!strcmp(s, "proportional"))
				data->transient_controller =
					AYLP_FSP_TRANSIENT_PROPORTIONAL;
			else if (!strcmp(s, "integral"))
				data->transient_controller = AYLP_FSP_TRANSIENT_INTEGRAL;
			else if (!strcmp(s, "pi"))
				data->transient_controller = AYLP_FSP_TRANSIENT_PI;
			else if (!strcmp(s, "modal"))
				data->transient_controller = AYLP_FSP_TRANSIENT_MODAL;
			else if (!strcmp(s, "hybrid"))
				data->transient_controller = AYLP_FSP_TRANSIENT_HYBRID;
			else {
				log_error("fsp: unknown transient_controller \"%s\".", s);
				return -1;
			}
		} else if (!strcmp(key, "transient_modal_q_scale")) {
			data->transient_modal_q_scale =
				json_object_get_double(val);
		} else if (!strcmp(key, "transient_modal_ab")) {
			data->transient_modal_ab =
				json_object_get_boolean(val);
		} else if (!strcmp(key, "transient_log")) {
			xfree(data->transient_log);
			data->transient_log =
				xstrdup(json_object_get_string(val));
		} else if (!strcmp(key, "transient_tau")) {
			data->transient_tau = json_object_get_double(val);
		} else if (!strcmp(key, "transient_hold")) {
			data->transient_hold =
				fabs(json_object_get_double(val));
		} else if (!strcmp(key, "transient_ramp")) {
			data->transient_ramp =
				fabs(json_object_get_double(val));
		} else if (!strcmp(key, "transient_clamp")) {
			double limit = fabs(json_object_get_double(val));
			data->transient_clamp_lo = -limit;
			data->transient_clamp_hi = limit;
		} else if (!strcmp(key, "transient_clamp_min")) {
			data->transient_clamp_lo = json_object_get_double(val);
		} else if (!strcmp(key, "transient_clamp_max")) {
			data->transient_clamp_hi = json_object_get_double(val);
		} else if (!strcmp(key, "transient_trip_command")) {
			data->transient_trip_command =
				fabs(json_object_get_double(val));
		} else if (!strcmp(key, "transient_max_duration")) {
			data->transient_max_duration =
				fabs(json_object_get_double(val));
		} else if (!strcmp(key, "transient_arm_delay")) {
			data->transient_arm_delay =
				fabs(json_object_get_double(val));
		} else if (!strcmp(key, "transient_shadow_mu")) {
			data->transient_shadow_mu = json_object_get_double(val);
		} else if (!strcmp(key, "transient_shadow_tau")) {
			data->transient_shadow_tau = json_object_get_double(val);
		} else if (!strcmp(key, "transient_shadow_min_duration")) {
			data->transient_shadow_min_duration = json_object_get_double(val);
		} else if (!strcmp(key, "transient_shadow_hold")) {
			data->transient_shadow_hold = json_object_get_double(val);
		} else if (!strcmp(key, "transient_shadow_ratio")) {
			data->transient_shadow_ratio = json_object_get_double(val);
		} else if (!strcmp(key, "transient_shadow_norm_ratio")) {
			data->transient_shadow_norm_ratio = json_object_get_double(val);
		} else if (!strcmp(key, "broad_freeze_closed")) {
			data->broad_freeze_closed = json_object_get_boolean(val);
		} else if (!strcmp(key, "trip_error")) {
			data->trip_error = fabs(json_object_get_double(val));
		} else if (!strcmp(key, "trip_command")) {
			data->trip_command = fabs(json_object_get_double(val));
		} else if (!strcmp(key, "trip_frames")) {
			data->trip_frames = json_object_get_uint64(val);
		} else if (!strcmp(key, "beam_recover_ramp")) {
			data->beam_recover_ramp =
				fabs(json_object_get_double(val));
		} else if (!strcmp(key, "guard_ratio")) {
			data->guard_ratio = json_object_get_double(val);
		} else if (!strcmp(key, "guard_floor")) {
			data->guard_floor = fabs(json_object_get_double(val));
		} else if (!strcmp(key, "guard_hold")) {
			data->guard_hold = fabs(json_object_get_double(val));
		} else if (!strcmp(key, "guard_ramp")) {
			data->guard_ramp = fabs(json_object_get_double(val));
		} else if (!strcmp(key, "guard_tick")) {
			data->guard_tick = fabs(json_object_get_double(val));
		} else if (!strcmp(key, "guard_min_cycles")) {
			data->guard_min_cycles = json_object_get_double(val);
		} else if (!strcmp(key, "gap_trip")) {
			if (json_object_is_type(val, json_type_string)
					&& !strcmp(json_object_get_string(val),
						"auto"))
				data->gap_trip_auto = true;
			else data->gap_trip = json_object_get_double(val);
		} else if (!strcmp(key, "y") || !strcmp(key, "axis_y")) {
			axy = val;
		} else if (!strcmp(key, "x") || !strcmp(key, "axis_x")) {
			axx = val;
		} else {
			log_warn("fsp: unknown parameter \"%s\"", key);
		}
	}

	if (data->type != AYLP_T_VECTOR) {
		log_error("fsp: type must be \"vector\".");
		return -1;
	}
	if (data->transient_mode == AYLP_FSP_TRANSIENT_OFF) {
		// Keep the tuning in the configuration, but make off unambiguous even
		// when an old profile still contains nonzero recovery settings.
		data->transient_sigma = 0.0;
		data->transient_shadow_mu = 0.0;
	} else if (data->transient_mode == AYLP_FSP_TRANSIENT_AUTO
			&& data->transient_sigma <= 0.0) {
		// Conservative, documented default.  The remaining parameters retain
		// their explicit values or the normal FSP defaults.
		data->transient_sigma = 6.0;
	}
	// Preserve the pre-selector behavior: configuring an event-only modal
	// gain used modal recovery unless a controller was named explicitly.
	if (!data->transient_controller_set
			&& data->transient_modal_q_scale > 0.0)
		data->transient_controller = AYLP_FSP_TRANSIENT_MODAL;
	// Resolve the output bounds: clamp_min/clamp_max win over the symmetric
	// `clamp` shorthand, whichever order they appeared in the config.
	if (isnan(data->clamp_lo)) data->clamp_lo = -data->clamp;
	if (isnan(data->clamp_hi)) data->clamp_hi = data->clamp;
	if (!isfinite(data->clamp_lo) || !isfinite(data->clamp_hi)
	|| data->clamp_lo >= data->clamp_hi) {
		log_error("fsp: clamp bounds must be finite with clamp_min < "
			"clamp_max (got %G, %G).",
			data->clamp_lo, data->clamp_hi);
		return -1;
	}
	// Zero must be reachable: the start_delay hold parks the command at 0,
	// `ramp` blends 0 -> full authority, and the burst guard and the
	// non-finite fallback both force u = 0. A window excluding 0 would make
	// all four of those unrepresentable.
	if (data->clamp_lo > 0.0 || data->clamp_hi < 0.0) {
		log_error("fsp: clamp window [%G, %G] must contain 0 -- the "
			"startup hold, ramp and burst guard all drive the "
			"command to zero.", data->clamp_lo, data->clamp_hi);
		return -1;
	}
	if (data->clamp_lo != -data->clamp_hi)
		log_info("fsp: asymmetric command limit [%G, %G] -- use this "
			"when the actuator's reachable range is not centred "
			"on the command origin", data->clamp_lo,
			data->clamp_hi);
	if (isnan(data->transient_clamp_lo))
		data->transient_clamp_lo = data->clamp_lo;
	if (isnan(data->transient_clamp_hi))
		data->transient_clamp_hi = data->clamp_hi;
	if (isnan(data->transient_trip_command))
		data->transient_trip_command = data->trip_command;
	if (!isfinite(data->transient_clamp_lo)
	|| !isfinite(data->transient_clamp_hi)
	|| data->transient_clamp_lo >= data->transient_clamp_hi
	|| data->transient_clamp_lo > data->clamp_lo
	|| data->transient_clamp_hi < data->clamp_hi
	|| data->transient_clamp_lo > 0.0
	|| data->transient_clamp_hi < 0.0) {
		log_error("fsp: transient clamp [%G, %G] must be finite, contain "
			"zero, and contain the normal clamp [%G, %G].",
			data->transient_clamp_lo, data->transient_clamp_hi,
			data->clamp_lo, data->clamp_hi);
		return -1;
	}
	if (data->delay < 1) {
		log_error("fsp: delay must be >= 1 sample.");
		return -1;
	}
	if (data->delay_frac < 0.0 || data->delay_frac >= 1.0) {
		log_error("fsp: delay_frac must satisfy 0 <= delay_frac < 1.");
		return -1;
	}
	if (data->delay_auto_max < 1 || data->delay_auto_max > 1024) {
		log_error("fsp: delay_auto_max must be in [1, 1024].");
		return -1;
	}
	if (!isfinite(data->delay_auto_settle) || data->delay_auto_settle < 0.0
	|| !isfinite(data->delay_auto_window) || data->delay_auto_window <= 0.0
	|| !isfinite(data->delay_auto_bias) || data->delay_auto_bias < 0.0) {
		log_error("fsp: auto delay needs delay_auto_settle >= 0, "
			"delay_auto_window > 0 and delay_auto_bias >= 0.");
		return -1;
	}
	if (data->fs <= 0.0) {
		log_error("fsp: fs must be > 0.");
		return -1;
	}
	if (data->delay_auto_bias_auto && !data->delay_ident_auto
			&& !(data->delay_ident_ms > 0.0)) {
		log_error("fsp: delay_auto_bias \"auto\" needs delay_ident_ms, "
			"the last identified transport delay in ms -- it is "
			"what the bias is recovered from.");
		return -1;
	}

	if (data->broad_order > 4096) {
		log_error("fsp: broad_order must be <= 4096.");
		return -1;
	}
	if (data->broad_order && (data->broad_mu <= 0.0
			|| data->broad_mu >= 2.0)) {
		log_error("fsp: broad_mu must satisfy 0 < broad_mu < 2.");
		return -1;
	}
	if (data->broad_order && data->broad_mu_init > 0.0
			&& (data->broad_mu_init < data->broad_mu
				|| data->broad_mu_init >= 2.0)) {
		log_error("fsp: broad_mu_init must satisfy broad_mu <= "
			"broad_mu_init < 2 (it is the FASTER initial step).");
		return -1;
	}
	if (data->broad_mu_init > 0.0 && data->broad_mu_tau <= 0.0) {
		log_error("fsp: broad_mu_tau must be > 0 when broad_mu_init is "
			"set.");
		return -1;
	}
	if (data->broad_order && (!isfinite(data->dark_train_min_real)
			|| data->dark_train_min_real <= 0.0
			|| data->dark_train_min_real > 1.0)) {
		log_error("fsp: dark_train_min_real must satisfy 0 < value <= 1.");
		return -1;
	}
	if (data->dark_predict_max > 4096) {
		log_error("fsp: dark_predict_max must be <= 4096 frames.");
		return -1;
	}
	if (!isfinite(data->transient_sigma)
	|| !isfinite(data->transient_floor)
	|| !isfinite(data->transient_settle_error)
	|| !isfinite(data->transient_kp)
	|| !isfinite(data->transient_ki)
	|| !isfinite(data->transient_i_leak)
	|| !isfinite(data->transient_i_limit)
	|| !isfinite(data->transient_tau)
	|| !isfinite(data->transient_hold)
	|| !isfinite(data->transient_ramp)
	|| !isfinite(data->transient_trip_command)
	|| !isfinite(data->transient_max_duration)
	|| !isfinite(data->transient_arm_delay)
	|| !isfinite(data->transient_shadow_mu)
	|| !isfinite(data->transient_shadow_tau)
	|| !isfinite(data->transient_shadow_min_duration)
	|| !isfinite(data->transient_shadow_hold)
	|| !isfinite(data->transient_shadow_ratio)
	|| !isfinite(data->transient_shadow_norm_ratio)
	|| data->transient_floor < 0.0
	|| data->transient_settle_error < 0.0
	|| data->transient_ki < 0.0
	|| data->transient_i_leak < 0.0
	|| data->transient_i_limit < 0.0
	|| data->transient_hold < 0.0
	|| data->transient_ramp < 0.0
	|| data->transient_trip_command < 0.0
	|| data->transient_max_duration < 0.0
	|| data->transient_arm_delay < 0.0
	|| data->transient_shadow_mu < 0.0
	|| data->transient_shadow_tau < 0.0
	|| data->transient_shadow_min_duration < 0.0
	|| data->transient_shadow_hold < 0.0
	|| data->transient_shadow_norm_ratio < 0.0) {
		log_error("fsp: transient parameters must be finite and floor, "
			"settle_error, Ki, I leak/limit, hold, ramp, trip, "
			"max_duration and arm_delay "
			"must be >= 0.");
		return -1;
	}
	if (data->transient_shadow_mu > 0.0
			&& (data->transient_sigma <= 0.0
				|| !data->broad_order
				|| data->transient_shadow_mu >= 2.0
				|| data->transient_shadow_tau <= 0.0
				|| data->transient_shadow_ratio <= 0.0
				|| data->transient_shadow_ratio >= 1.0
				|| data->transient_shadow_norm_ratio <= 0.0)) {
		log_error("fsp: transient shadow recovery requires transient_sigma > 0, "
			"broad_order > 0, "
			"0 < shadow_mu < 2, shadow_tau > 0, 0 < shadow_ratio < 1, "
			"and shadow_norm_ratio > 0.");
		return -1;
	}
	if (data->transient_sigma > 0.0 && data->transient_tau <= 0.0) {
		log_error("fsp: enabled transient recovery requires "
			"transient_tau > 0.");
		return -1;
	}
	if (data->transient_sigma > 0.0
			&& data->transient_controller != AYLP_FSP_TRANSIENT_INTEGRAL
			&& (data->transient_kp <= 0.0
				|| data->transient_kp >= 1.0)) {
		log_error("fsp: proportional, PI, modal and hybrid transient "
			"controllers require 0 < transient_kp < 1.");
		return -1;
	}
	if (data->transient_sigma > 0.0
			&& fsp_transient_uses_integral(data->transient_controller)
			&& (data->transient_ki <= 0.0
				|| data->transient_i_limit <= 0.0)) {
		log_error("fsp: integral, PI and hybrid transient controllers "
			"require transient_ki > 0 and transient_i_limit > 0.");
		return -1;
	}
	if (data->transient_modal_q_scale < 0.0
			|| !isfinite(data->transient_modal_q_scale)) {
		log_error("fsp: transient_modal_q_scale must be finite and >= 0.");
		return -1;
	}
	if (fsp_transient_uses_modal(data->transient_controller)
			&& data->transient_modal_q_scale <= 0.0) {
		log_error("fsp: modal and hybrid transient controllers require "
			"transient_modal_q_scale > 0.");
		return -1;
	}
	if (data->transient_modal_ab
			&& !fsp_transient_uses_modal(data->transient_controller)) {
		log_error("fsp: transient_modal_ab requires transient_controller "
			"modal or hybrid.");
		return -1;
	}
	if (data->broad_lp) {
		if (!data->broad_order) {
			log_error("fsp: broad_lp requires broad_order > 0.");
			return -1;
		}
		if (data->broad_lp < 3 || data->broad_lp > 63) {
			log_error("fsp: broad_lp must be 3..63 taps.");
			return -1;
		}
		if (!(data->broad_lp & 1)) {
			data->broad_lp += 1;
			log_warn("fsp: broad_lp must be odd for an integer "
				"group delay; using %zu", data->broad_lp);
		}
	}
	data->broad_gd = data->broad_lp ? (data->broad_lp - 1) / 2 : 0;
	if ((data->trip_error > 0.0 || data->trip_command > 0.0)
			&& !data->trip_frames) {
		log_error("fsp: trip_frames must be >= 1 when a trip is enabled.");
		return -1;
	}
	if (!isfinite(data->beam_recover_ramp)) {
		log_error("fsp: beam_recover_ramp must be finite and >= 0.");
		return -1;
	}
	if (data->drift_order < 1 || data->drift_order > 2) {
		log_error("fsp: drift_order must be 1 or 2.");
		return -1;
	}
	if (!isfinite(data->precenter_delay) || data->precenter_delay < 0.0
			|| !isfinite(data->precenter_ramp) || data->precenter_ramp < 0.0
			|| !isfinite(data->precenter_clamp)
			|| (data->precenter_clamp > 0.0
				&& (data->drift_tau <= 0.0
					|| data->start_delay <= data->precenter_delay))) {
		log_error("fsp: precentering needs finite non-negative settings, "
			"drift_tau > 0, and start_delay > precenter_delay");
		return -1;
	}
	if (!axy || !axx) {
		log_error("fsp: both \"y\" and \"x\" axis objects are required.");
		return -1;
	}
	if (fsp_parse_axis(&data->axis[0], axy)) return -1;
	if (fsp_parse_axis(&data->axis[1], axx)) return -1;

	// resolve per-axis delays and clamp bounds: unset values inherit global
	for (int a = 0; a < 2; a++) {
		struct aylp_fsp_axis *ax = &data->axis[a];
		if (!ax->delay && !ax->delay_auto) {
			// no delay of its own: inherit the global one, "auto"
			// included
			ax->delay = data->delay;
			ax->delay_auto = data->delay_auto;
		} else if (ax->delay_auto) {
			ax->delay = data->delay;	// provisional
		}
		if (ax->delay_auto && ax->delay > data->delay_auto_max)
			ax->delay = data->delay_auto_max;
		if (isnan(ax->clamp_lo)) ax->clamp_lo = data->clamp_lo;
		if (isnan(ax->clamp_hi)) ax->clamp_hi = data->clamp_hi;
		if (!isfinite(ax->clamp_lo) || !isfinite(ax->clamp_hi)
		|| ax->clamp_lo >= ax->clamp_hi) {
			log_error("fsp: %s axis clamp bounds must be finite "
				"with clamp_min < clamp_max (got %G, %G).",
				a == 0 ? "y" : "x", ax->clamp_lo, ax->clamp_hi);
			return -1;
		}
		if (ax->clamp_lo > 0.0 || ax->clamp_hi < 0.0) {
			log_error("fsp: %s axis clamp window [%G, %G] must "
				"contain 0 -- the startup hold, ramp and burst "
				"guard all drive the command to zero.",
				a == 0 ? "y" : "x", ax->clamp_lo, ax->clamp_hi);
			return -1;
		}
		if (data->transient_clamp_lo > ax->clamp_lo
		|| data->transient_clamp_hi < ax->clamp_hi) {
			log_error("fsp: transient clamp [%G, %G] must contain the %s "
				"axis clamp [%G, %G].", data->transient_clamp_lo,
				data->transient_clamp_hi, a == 0 ? "y" : "x",
				ax->clamp_lo, ax->clamp_hi);
			return -1;
		}
		if (ax->clamp_lo != -ax->clamp_hi)
			log_info("fsp: %s axis asymmetric command limit "
				"[%G, %G]", a == 0 ? "y" : "x",
				ax->clamp_lo, ax->clamp_hi);
		if (ax->delay_frac < 0.0) ax->delay_frac = data->delay_frac;
		if (ax->delay_frac >= 1.0) {
			log_error("fsp: %s axis delay_frac must satisfy "
				"0 <= delay_frac < 1.", a == 0 ? "y" : "x");
			return -1;
		}
		// Every delay-dependent ring is sized here, once. An auto delay
		// is not known until the loop has run, so its rings are sized
		// for the ceiling instead and the estimate is clamped to it.
		ax->delay_alloc = ax->delay_auto
			? data->delay_auto_max : ax->delay;
	}

	// delay_auto now means "a startup measurement is still needed", which is
	// also true when only one axis asked for it, or when nothing asked for an
	// auto delay but fs or the bias is resolved from the same measurement
	data->delay_auto = data->axis[0].delay_auto || data->axis[1].delay_auto
		|| data->fs_auto || data->delay_auto_bias_auto
		|| data->gap_trip_auto;
	if (data->delay_auto) {
		double need = data->delay_auto_settle + data->delay_auto_window;
		if (data->start_delay < need) {
			log_warn("fsp: raising start_delay %G -> %G s: the auto "
				"delay is measured during the hold and has to "
				"be installed before the loop closes.",
				data->start_delay, need);
			data->start_delay = need;
		}
		// one period sample per frame, one compute sample per timing
		// window published by the source
		size_t cap = (size_t)(data->delay_auto_window * data->fs) + 16;
		if (cap > 8192) cap = 8192;
		data->auto_cap_dt = cap;
		data->auto_cap_c = 64;
		data->auto_dt = xcalloc(cap, sizeof(double));
		data->auto_compute = xcalloc(data->auto_cap_c, sizeof(double));
		log_info("fsp: delay is auto: measuring for %G s after a %G s "
			"settle, then locking in at the end of the %G s hold "
			"(ceiling %zu frames, assumed sensor/ZOH bias %G "
			"frames)", data->delay_auto_window,
			data->delay_auto_settle, data->start_delay,
			data->delay_auto_max, data->delay_auto_bias);
	} else {
		data->delay_auto_done = true;
	}

	if (data->cmd_fc > 0.0) {
		if (data->cmd_fc >= data->fs / 2.0) {
			log_error("fsp: cmd_fc must be < fs/2.");
			return -1;
		}
		fsp_build_cmdlp(data);
		log_info("fsp: command robustness low-pass at %G Hz, per-mode "
			"phase/gain pre-compensated", data->cmd_fc);
	}

	if (data->guard_ratio > 0.0) {
		// fast envelope ~10 ms (rides the burst onset), baseline ~10 s
		data->guard_beta_fast = 1.0 - exp(-1.0 / (0.010 * data->fs));
		data->guard_beta_slow = 1.0 - exp(-1.0 / (10.0 * data->fs));
	}
	// stall-gap DC pad: ~0.5 s EWMA -- slow enough to average the
	// vibration lines out, fast enough to track intra-run drift
	data->gap_dc_beta = 1.0 - exp(-1.0 / (0.5 * data->fs));
	// frequency-demodulator time constant: sized to adapt_df_max (a
	// first-order PLL's capture range is ~1/(2*pi*tau)), not to
	// adapt_tau, so a large sudden line shift is actually within the
	// range the correction cap was meant to chase, instead of reading
	// out a near-arbitrary phase once adapt_tau makes the capture range
	// far narrower than adapt_df_max advertises.
	data->demod_tau = data->adapt_df_max > 0.0
		? 1.0 / (2.0 * M_PI * data->adapt_df_max)
		: (data->adapt_tau > 0.0 ? data->adapt_tau : 5.0);
	data->demod_beta = 1.0 - exp(-1.0 / (data->demod_tau * data->fs));
	if (data->adapt_period > 0.0 && data->adapt_df_max > 0.0)
		log_info("fsp: frequency-demod time constant %.3G s (capture "
			"range ~%.3G Hz, cap %G Hz/update)", data->demod_tau,
			1.0 / (2.0 * M_PI * data->demod_tau), data->adapt_df_max);
	if (data->drift_tau > 0.0)
		data->drift_beta =
			1.0 - exp(-1.0 / (data->drift_tau * data->fs));
	if (data->transient_sigma > 0.0)
		data->transient_beta =
			1.0 - exp(-1.0 / (data->transient_tau * data->fs));
	if (data->transient_shadow_mu > 0.0)
		data->transient_shadow_beta = 1.0
			- exp(-1.0 / (data->transient_shadow_tau * data->fs));
	if (data->transient_log) {
		data->transient_log_fp = fopen(data->transient_log, "w");
		if (!data->transient_log_fp) {
			log_error("fsp: fopen %s: %m", data->transient_log);
			return -1;
		}
		fprintf(data->transient_log_fp,
			"axis,event,recovery,start_s,error_settle_s,innovation_quiet_s,"
			"recovery_s,peak_abs_error,error_rms,peak_abs_command,"
			"peak_abs_integral,samples,shadow_promoted,shadow_error_ratio\n");
		fflush(data->transient_log_fp);
	}

	for (int a = 0; a < 2; a++) {
		struct aylp_fsp_axis *ax = &data->axis[a];
		fsp_build_modes(ax, data->fs);
		fsp_build_comp(ax, data);
		if (data->guard_ratio > 0.0) fsp_build_guard(ax, data);
		if (fsp_solve_gain(ax)) {
			log_error("fsp: Riccati solve failed on %s axis; check "
				"q/r and mode params.", a == 0 ? "y" : "x");
			return -1;
		}
		if (data->transient_modal_q_scale > 0.0
				&& fsp_solve_transient_gain(ax,
					data->transient_modal_q_scale)) {
			log_error("fsp: transient modal Riccati solve failed on %s "
				"axis; check transient_modal_q_scale and mode params.",
				a == 0 ? "y" : "x");
			return -1;
		}
		ax->transient_gain_current =
			data->transient_modal_q_scale > 0.0;
		// One extra command is retained for fractional-delay interpolation.
		// delay_alloc, not delay: with an auto delay these lengths are
		// fixed before the delay is known (see fsp.h).
		ax->ucmd = xcalloc(ax->delay_alloc + 1, sizeof(double));
		if (data->broad_order) {
			// twice the horizon the main filter runs
			ax->dark_predict_max = data->dark_predict_max
				? data->dark_predict_max
				: 2 * (ax->delay_alloc + data->broad_gd);
			if (ax->dark_predict_max < 4) ax->dark_predict_max = 4;
			ax->dark_bank_max = ax->dark_predict_max
				+ ax->delay_alloc + data->broad_gd + 1;
			// The dark bank's longest training window ends dark_bank_max
			// samples back, which can reach deeper
			// into the ring than the main filter's delay + broad_gd
			// window does; the ring must hold whichever is longer.
			size_t reach = ax->delay_alloc + data->broad_gd;
			if (data->dark_predict && ax->dark_bank_max > reach)
				reach = ax->dark_bank_max;
			ax->broad_hist_len = data->broad_order + reach + 2;
			ax->broad_hist = xcalloc(ax->broad_hist_len,
				sizeof(double));
			ax->broad_hist_bl = xcalloc(ax->broad_hist_len,
				sizeof(double));
			if (data->dark_predict)
				ax->dark_w = xcalloc(ax->dark_bank_max
					* data->broad_order, sizeof(double));
			ax->broad_w = xcalloc(data->broad_order, sizeof(double));
			ax->broad_w_next = xcalloc(data->broad_order,
				sizeof(double));
			if (data->transient_shadow_mu > 0.0) {
				ax->shadow_w = xcalloc(data->broad_order, sizeof(double));
				ax->shadow_w_next = xcalloc(data->broad_order,
					sizeof(double));
			}
			ax->broad_xbuf = xcalloc(data->broad_order, sizeof(double));
			ax->broad_xbuf2 = xcalloc(data->broad_order, sizeof(double));
			ax->broad_syn = xcalloc(ax->broad_hist_len, 1);
			ax->broad_synbuf = xcalloc(data->broad_order, 1);
			ax->broad_synbuf2 = xcalloc(data->broad_order, 1);
			if (data->broad_lp)
				ax->broad_lpbuf = xcalloc(data->broad_lp,
					sizeof(double));
			if (!data->dark_bank_train)
				data->dark_bank_train = 4;
			if (data->dark_bank_train > ax->dark_bank_max)
				data->dark_bank_train = ax->dark_bank_max;
		}
		ax->r_ewma = ax->r;
		ax->transient_var = 0.0;
		for (size_t i = 0; i < ax->n_modes; i++) {
			// q_ewma tracks the mode's STATE energy; seed it at the
			// stationary energy implied by the configured drive q
			ax->q_ewma[i] = ax->q[i] * ax->Gv[i];
			ax->demod_re[i] = ax->demod_im[i] = 0.0;
			ax->demod_ph[i] = 0.0;
		}
		log_info("fsp: %s axis, %zu modes, K=%G, predicting %zu+%.3G "
			"samples ahead%s%s", a == 0 ? "y" : "x", ax->n_modes,
			ax->K, ax->delay, ax->delay_frac,
			ax->delay_auto ? " (provisional; auto delay lands at "
				"the end of the hold)" : "",
			data->adapt_period > 0 ? ", adaptive" : " (fixed)");
		if (ax->plant_shaped)
			log_info("fsp: %s Bode-shaped plant H(z): "
				"b=[%G,%G,%G], a=[1,%G,%G], matched stable inverse",
				a == 0 ? "y" : "x",
				ax->plant_b[0], ax->plant_b[1], ax->plant_b[2],
				ax->plant_a[1], ax->plant_a[2]);
	}
	data->broad_mu_cur = data->broad_mu_init > data->broad_mu
		? data->broad_mu_init : data->broad_mu;
	if (fsp_load_wiener(data)) return -1;
	if (data->wiener_out) {
		// Fail now, not after a 1700 s run: the dump happens in fini,
		// so an unwritable path would otherwise cost the whole run's
		// weights and only announce itself at the very end.
		FILE *probe = fopen(data->wiener_out, "w");
		if (!probe) {
			log_error("fsp: wiener_out \"%s\" is not writable: %s",
				data->wiener_out, strerror(errno));
			// fini still runs on the init-failure path; drop the
			// target so it does not repeat the same error
			xfree(data->wiener_out);
			data->wiener_out = NULL;
			return -1;
		}
		fclose(probe);
		log_info("fsp: will save learned %zu-tap predictors to %s at "
			"exit", data->broad_order, data->wiener_out);
	}
	if (data->wiener_trace) {
		if (!data->broad_order) {
			log_warn("fsp: wiener_trace ignored (broad_order = 0)");
			xfree(data->wiener_trace);
			data->wiener_trace = NULL;
		} else {
			data->wiener_trace_fp = fopen(data->wiener_trace, "w");
			if (!data->wiener_trace_fp) {
				log_error("fsp: wiener_trace \"%s\" is not "
					"writable: %s", data->wiener_trace,
					strerror(errno));
				xfree(data->wiener_trace);
				data->wiener_trace = NULL;
				return -1;
			}
			if (data->wiener_trace_period <= 0.0)
				data->wiener_trace_period = 10.0;
			for (int a = 0; a < 2; a++)
				data->wtrace_prev[a] = xcalloc(
					data->broad_order, sizeof(double));
			fprintf(data->wiener_trace_fp,
				"# fsp observer convergence trace, every %G s\n"
				"# broad_order=%zu broad_mu=%G broad_lp=%zu "
				"fs=%G start_delay=%G ramp=%G\n"
				"# ||dw|| is the tap change since the previous "
				"sample: it flattening is the observer having "
				"converged\n"
				"# t_s n_closed y_wnorm y_wnextnorm y_dwnorm "
				"x_wnorm x_wnextnorm x_dwnorm y_drift x_drift "
				"y_guard x_guard y_trans x_trans mu\n",
				data->wiener_trace_period, data->broad_order,
				data->broad_mu, data->broad_lp, data->fs,
				data->start_delay, data->ramp);
			fflush(data->wiener_trace_fp);
			log_info("fsp: observer convergence trace every %G s "
				"to %s", data->wiener_trace_period,
				data->wiener_trace);
		}
	}
	if (data->broad_order)
		log_info("fsp: full-band %zu-state disturbance predictor, horizon "
			"y %zu+%.3G / x %zu+%.3G frames, NLMS mu %G",
			data->broad_order,
			data->axis[0].delay, data->axis[0].delay_frac,
			data->axis[1].delay, data->axis[1].delay_frac,
			data->broad_mu);
	if (data->broad_order && data->broad_mu_init > data->broad_mu)
		log_info("fsp: NLMS step schedule: mu %G -> %G with %G s time "
			"constant from loop close (nominal constant "
			"order/(mu*fs) = %.2G s at the initial step, %.2G s at "
			"the final)", data->broad_mu_init, data->broad_mu,
			data->broad_mu_tau,
			data->broad_order / (data->broad_mu_init * data->fs),
			data->broad_order / (data->broad_mu * data->fs));
	if (data->broad_lp)
		log_info("fsp: observer band-limit: %zu-tap boxcar prefilter "
			"(first null %.0f Hz), +%zu frames folded into the broad "
			"horizon", data->broad_lp, data->fs / (double)data->broad_lp,
			data->broad_gd);
	else if (data->broad_order)
		log_warn("fsp: observer band-limit DISABLED (broad_lp = 0): the "
			"NLMS can learn/chase HF command echo from K/delay "
			"mismatch (see the 2026-07-22 380 Hz ring)");
	if (data->drift_tau > 0.0)
		log_info("fsp: slow drift separated from vibration prediction "
			"with %G s %s", data->drift_tau,
			data->drift_order == 2 ? "position/rate observer" : "EWMA");
	if (data->transient_sigma > 0.0)
		log_info("fsp: transient path on: trigger %G sigma (floor %G), "
			"controller %s, Kp %G, Ki %G/s, I leak %G/s, I limit %G, "
			"output-settle band %G, quiet hold %G s, "
			"return ramp %G s, armed %G s after close, event clamp "
			"[%G, %G], event command trip %G, timeout %G s%s",
			data->transient_sigma, data->transient_floor,
			fsp_transient_controller_name(data->transient_controller),
			data->transient_kp, data->transient_ki,
			data->transient_i_leak, data->transient_i_limit,
			data->transient_settle_error,
			data->transient_hold,
			data->transient_ramp,
			data->transient_arm_delay,
			data->transient_clamp_lo, data->transient_clamp_hi,
			data->transient_trip_command,
			data->transient_max_duration,
			data->transient_modal_q_scale > 0.0
				? ", event-only modal observer enabled" : "");
	if (data->transient_shadow_mu > 0.0)
		log_info("fsp: event shadow NLMS on: mu %G, validation tau %G s, "
			"minimum event %G s, error ratio <= %G for %G s, coefficient "
			"change <= %Gx primary norm", data->transient_shadow_mu,
			data->transient_shadow_tau,
			data->transient_shadow_min_duration,
			data->transient_shadow_ratio,
			data->transient_shadow_hold,
			data->transient_shadow_norm_ratio);
	if (data->trip_error > 0.0 || data->trip_command > 0.0)
		log_info("fsp: non-latching limit diagnostics: error=%G command=%G "
			"for %zu frames; upstream beam loss holds zero, then ramps "
			"back over %G s",
			data->trip_error, data->trip_command,
			data->trip_frames, data->beam_recover_ramp);
	if (data->guard_ratio > 0.0)
		log_info("fsp: burst guard on: y %.0f Hz / x %.0f Hz detectors, "
			"trigger %Gx over max(baseline, %G) sustained %zu/%zu "
			"samples (%G cycles), hold %G s + ramp %G s, ticker "
			"every %G s", data->axis[0].gd_f0, data->axis[1].gd_f0,
			data->guard_ratio, data->guard_floor,
			data->axis[0].gd_min_samples, data->axis[1].gd_min_samples,
			data->guard_min_cycles, data->guard_hold, data->guard_ramp,
			data->guard_tick);
	else
		log_warn("fsp: burst guard DISABLED (guard_ratio <= 0)");
	if (data->gap_trip > 0.0)
		log_info("fsp: stall-gap patch on: frame gap > %G ms pads the "
			"histories with the held command and continues",
			1e3 * data->gap_trip);
	else
		log_warn("fsp: stall-gap patch DISABLED (gap_trip <= 0)");

	// --- RT hygiene (2026-07-21): lock memory and warm the adaptation path.
	// Under SCHED_FIFO the first re-identification after the loop closes was
	// bursting ~2.8 ms/axis (two back-to-back frame gaps ~8 s in) because its
	// code, stack, and workspace pages were cold / not resident. Lock the
	// process into RAM (no demand-paging faults on the RT loop) and run the
	// exact runtime adapt+solve sequence once here, off the critical path, so
	// every page it touches is already warm when adapt_period first fires.
	if (mlockall(MCL_CURRENT | MCL_FUTURE))
		log_warn("fsp: mlockall failed (%s); page faults may still jitter "
			"the loop -- run with privilege / raised memlock rlimit",
			strerror(errno));
	else
		log_info("fsp: locked process memory (mlockall) against RT page "
			"faults");
	if (data->adapt_period > 0.0) {
		fsp_adapt(data);	// build_modes/comp + solve_begin, both axes
		for (int a = 0; a < 2; a++) {
			struct aylp_fsp_axis *ax = &data->axis[a];
			while (ax->adapt_solving
					&& fsp_adapt_advance(ax, data, (size_t)a,
						AYLP_FSP_RICCATI_MAXIT) == 0)
				;
		}
		log_info("fsp: warmed adaptation path (first re-identification "
			"pre-touched)");
	}

	self->type_in = AYLP_T_VECTOR;
	self->units_in = AYLP_U_ANY;
	self->type_out = AYLP_T_UNCHANGED;
	self->units_out = data->units;
	return 0;
}


// slow adaptation tick: refresh q_i (from projected innovation energy) and r
// (from the broadband innovation floor), optionally nudge each center freq
// toward the locally demodulated line, then rebuild coeffs and the gain.
static void fsp_adapt(struct aylp_fsp_data *data)
{
	for (int a = 0; a < 2; a++) {
		struct aylp_fsp_axis *ax = &data->axis[a];
		for (size_t i = 0; i < ax->n_modes; i++) {
			// process noise tracks the mode's recent energy, with a
			// floor so a quiet line can come back. q_ewma is STATE
			// energy; the Riccati wants the DRIVE variance, smaller
			// by the mode's stationary variance gain Gv (~1e4 for a
			// sharp line). Feeding state energy in directly inflates
			// q/r by that factor and re-creates the ~3x waterbed the
			// fsp_sim.py study measured.
			double qn = ax->q_ewma[i] / ax->Gv[i];
			if (qn < 1e-12) qn = 1e-12;
			ax->q[i] = qn;
			// bounded frequency retune from the demod phasor angle
			if (data->adapt_df_max > 0.0) {
				double mag = hypot(ax->demod_re[i],
					ax->demod_im[i]);
				if (mag > 1e-9) {
					double ang = atan2(ax->demod_im[i],
						ax->demod_re[i]);
					// ang is the residual phase over the
					// window; convert to a small df, capped.
					// Must divide by the SAME time constant
					// the phasor was accumulated with
					// (demod_tau, sized to adapt_df_max),
					// not adapt_tau -- otherwise the readout
					// and the conversion disagree on what
					// window "ang" was measured over.
					double df = ang / (2.0 * M_PI)
						* (1.0 / data->demod_tau);
					if (df > data->adapt_df_max)
						df = data->adapt_df_max;
					if (df < -data->adapt_df_max)
						df = -data->adapt_df_max;
					ax->f[i] += df;
					if (ax->f[i] < 0.1) ax->f[i] = 0.1;
					if (ax->f[i] > 0.45 * data->fs)
						ax->f[i] = 0.45 * data->fs;
				}
			}
			ax->demod_re[i] = ax->demod_im[i] = 0.0;
		}
		if (ax->r_ewma > 1e-12) ax->r = ax->r_ewma;
		fsp_build_modes(ax, data->fs);
		fsp_build_comp(ax, data);
		// The mode coefficients and q/r have changed, so the previous fast
		// event gain is no longer paired with the active physical model.
		if (data->transient_modal_q_scale > 0.0)
			ax->transient_gain_current = false;
		ax->adapt_transient_solving = false;
		// Kick off the resumable Riccati solve instead of running it
		// here: fsp_proc advances it a few iterations per frame so it
		// never bursts. The modes/comp above take effect immediately;
		// the previous gain keeps running until the new one converges
		// (a few tens of ms), a negligible mismatch since the modes
		// moved at most adapt_df_max Hz.
		if (fsp_solve_gain_begin(ax)) {
			log_warn("fsp: adaptive Riccati init failed on %s "
				"axis; keeping previous gain.",
				a == 0 ? "y" : "x");
			ax->adapt_solving = false;
		} else {
			ax->adapt_solving = true;
		}
	}
	log_debug("fsp: re-identification started (gain solve amortized).");
}


// Advance one axis's serial normal-then-transient Riccati refresh.  The event
// solve uses the same persistent workspace, so it remains bounded per frame and
// does not reintroduce the multi-millisecond real-time burst the amortization
// was added to remove.
static int fsp_adapt_advance(struct aylp_fsp_axis *ax,
	struct aylp_fsp_data *data, size_t axis, size_t niter)
{
	int done = fsp_solve_gain_iterate(ax, niter);
	if (done == 0) return 0;
	const char *name = axis == 0 ? "y" : "x";
	if (done < 0) {
		log_warn("fsp: adaptive %s Riccati solve failed on %s axis; "
			"keeping previous gain.", ax->adapt_transient_solving
				? "transient" : "normal", name);
		ax->adapt_solving = false;
		ax->adapt_transient_solving = false;
		return -1;
	}

	if (!ax->adapt_transient_solving) {
		if (fsp_solve_gain_finalize(ax)) {
			log_warn("fsp: adaptive normal Riccati solve failed on %s "
				"axis; keeping previous gain.", name);
			ax->adapt_solving = false;
			return -1;
		}
		if (data->transient_modal_q_scale <= 0.0) {
			ax->adapt_solving = false;
			return 1;
		}
		if (fsp_solve_gain_begin_scaled(ax,
				data->transient_modal_q_scale)) {
			log_warn("fsp: adaptive transient Riccati init failed on %s "
				"axis; modal event entry remains inhibited.", name);
			ax->adapt_solving = false;
			return -1;
		}
		ax->adapt_transient_solving = true;
		return 0;
	}

	// finalize() targets the normal gain and posterior variance. Preserve those
	// already-refreshed values while using the common finalizer for transient_L.
	double normal_L[AYLP_FSP_MAX_DIM];
	double normal_post[AYLP_FSP_MAX_MODES];
	memcpy(normal_L, ax->L, ax->dim * sizeof(double));
	memcpy(normal_post, ax->post_var, ax->n_modes * sizeof(double));
	int err = fsp_solve_gain_finalize(ax);
	if (!err)
		memcpy(ax->transient_L, ax->L, ax->dim * sizeof(double));
	memcpy(ax->L, normal_L, ax->dim * sizeof(double));
	memcpy(ax->post_var, normal_post, ax->n_modes * sizeof(double));
	ax->adapt_solving = false;
	ax->adapt_transient_solving = false;
	if (err) {
		log_warn("fsp: adaptive transient Riccati solve failed on %s axis; "
			"modal event entry remains inhibited.", name);
		return -1;
	}
	ax->transient_gain_current = true;
	log_debug("fsp: adaptive transient gain refreshed on %s axis.", name);
	return 1;
}


// Blended delay / delay+1 broadband prediction over the newest P-sample
// window of broad_hist. Callers check isfinite.
static double fsp_broad_predict(const struct aylp_fsp_axis *ax,
	const struct aylp_fsp_data *data)
{
	size_t P = data->broad_order, H = ax->broad_hist_len;
	size_t idx = ax->broad_head;
	double b = 0.0, bn = 0.0;
	for (size_t i = 0; i < P; i++) {
		b += ax->broad_w[i] * ax->broad_hist[idx];
		bn += ax->broad_w_next[i] * ax->broad_hist[idx];
		idx = idx ? idx - 1 : H - 1;
	}
	return (1.0 - ax->delay_frac) * b + ax->delay_frac * bn;
}


static void fsp_reset_after_beam_loss(struct aylp_fsp_axis *ax,
	struct aylp_fsp_data *data)
{
	memset(ax->xhat, 0, ax->dim * sizeof(double));
	memset(ax->ucmd, 0, (ax->delay + 1) * sizeof(double));
	ax->uhead = 0;
	ax->frac_x1 = ax->frac_y1 = 0.0;
	ax->plant_z1 = ax->plant_z2 = 0.0;
	ax->plant_iz1 = ax->plant_iz2 = 0.0;
	ax->lp_z1 = ax->lp_z2 = 0.0;
	ax->drift_hat = 0.0;
	ax->drift_rate = 0.0;
	ax->phi_dc = 0.0;
	ax->transient_active = false;
	ax->transient_recovering = false;
	ax->transient_blocked = false;
	ax->shadow_active = false;
	ax->transient_i = 0.0;
	ax->transient_i_prev = 0.0;
	ax->transient_i_step = 0.0;
	ax->guard_active = false;
	ax->guard_ramping = false;
	ax->gd_over_count = 0;
	ax->broad_dark_run = 0;
	if (data->broad_order) {
		memset(ax->broad_hist, 0,
			ax->broad_hist_len * sizeof(double));
		memset(ax->broad_hist_bl, 0,
			ax->broad_hist_len * sizeof(double));
		// a zeroed history is fabricated, not measured; broad_seen 0
		// already gates the predictor, but the flags have to agree
		memset(ax->broad_syn, 1, ax->broad_hist_len);
		ax->broad_head = 0;
		ax->broad_seen = 0;
		ax->broad_fab = ax->broad_hist_len;
		ax->broad_lpclean = 0;
		if (data->broad_lp) {
			memset(ax->broad_lpbuf, 0,
				data->broad_lp * sizeof(double));
			ax->broad_lphead = 0;
		}
	}
}


// ---- auto delay ------------------------------------------------------------
// See fsp.h for what is measured and why the compute term has to come from the
// source device rather than from here.

static int fsp_cmp_double(const void *a, const void *b)
{
	double x = *(const double *)a, y = *(const double *)b;
	return x < y ? -1 : (x > y ? 1 : 0);
}


// median, in place (the caller is done with the array either way)
static double fsp_median(double *v, size_t n)
{
	qsort(v, n, sizeof *v, fsp_cmp_double);
	return n & 1 ? v[n/2] : 0.5 * (v[n/2 - 1] + v[n/2]);
}


/** Rebuild everything that was computed from fs, after fs itself is measured.
 * Same pieces init builds, in the same order: the mode AR coefficients (and so
 * their Riccati gains), the per-mode command pre-compensation, the burst
 * guard's regeneration frequency, the command low-pass, and the per-sample EWMA
 * weights that were converted from time constants. Runs once, inside the
 * startup hold, so a Riccati solve here costs a frame nobody is using. */
static void fsp_rebuild_for_fs(struct aylp_fsp_data *data)
{
	for (int a = 0; a < 2; a++) {
		struct aylp_fsp_axis *ax = &data->axis[a];
		fsp_build_modes(ax, data->fs);
		fsp_build_comp(ax, data);
		if (data->guard_ratio > 0.0) fsp_build_guard(ax, data);
		if (fsp_solve_gain(ax))
			log_error("fsp: Riccati solve failed on the %s axis at "
				"the measured fs; the previous gain stays in "
				"place", a == 0 ? "y" : "x");
		if (data->transient_modal_q_scale > 0.0
				&& fsp_solve_transient_gain(ax,
					data->transient_modal_q_scale))
			log_error("fsp: transient modal solve failed on the %s "
				"axis at the measured fs", a == 0 ? "y" : "x");
	}
	if (data->cmd_fc > 0.0) fsp_build_cmdlp(data);
	if (data->drift_tau > 0.0)
		data->drift_beta = 1.0 - exp(-1.0 / (data->drift_tau * data->fs));
	if (data->transient_sigma > 0.0)
		data->transient_beta =
			1.0 - exp(-1.0 / (data->transient_tau * data->fs));
	if (data->transient_shadow_mu > 0.0)
		data->transient_shadow_beta = 1.0
			- exp(-1.0 / (data->transient_shadow_tau * data->fs));
}


/** One measurement sample, taken per frame during the startup hold. */
static void fsp_delay_auto_sample(struct aylp_fsp_data *data, double now)
{
	// the first frames of a run are the source's slow first frame, cold
	// caches and page faults; none of them are the steady state
	if (now - data->t0 < data->delay_auto_settle) return;
	if (data->t_last > 0.0 && data->auto_n_dt < data->auto_cap_dt) {
		double dt = now - data->t_last;
		if (dt > 0.0) data->auto_dt[data->auto_n_dt++] = dt;
	}
	struct aylp_loop_timing t;
	// each publication is already a mean over its own window; record it
	// once, not once per frame, so the median isn't weighted by how long
	// each window happened to stay current
	if (aylp_timing_get(&t) && t.compute_time != data->auto_last_c
			&& data->auto_n_c < data->auto_cap_c) {
		data->auto_last_c = t.compute_time;
		data->auto_compute[data->auto_n_c++] = t.compute_time;
	}
}


/** Turn the samples into a delay and install it. Called once, inside the hold,
 * as soon as the measurement window closes: from here on the delay is in use. */
static void fsp_delay_auto_lock(struct aylp_fsp_data *data)
{
	data->delay_auto_done = true;
	const char *why = 0;
	if (!data->auto_n_c)
		why = "no device in this pipeline published loop timing -- the "
			"source must be one that measures its own idle time "
			"(see libaylp/timing.h)";
	else if (data->auto_n_dt < 32)
		why = "too few frames arrived during the hold to measure a loop "
			"period";
	if (why) {
		log_error("fsp: AUTO DELAY FAILED: %s. Falling back to %zu+%G "
			"frames, which is a default and not a measurement of "
			"this bench. Stop the run rather than close a loop on "
			"it.", why, data->axis[0].delay,
			data->axis[0].delay_frac);
	} else {
		double dt = fsp_median(data->auto_dt, data->auto_n_dt);
		double comp = fsp_median(data->auto_compute, data->auto_n_c);
		double fs_meas = 1.0 / dt;
		if (data->fs_auto) {
			// Everything derived from fs is rebuilt below, so this has
			// to happen before the delay is computed in frames.
			log_info("fsp: auto fs -> %.2f Hz (was %G configured)",
				fs_meas, data->fs);
			data->fs = fs_meas;
			fsp_rebuild_for_fs(data);
		}
		if (data->gap_trip_auto) {
			// 1.5 frames: a normal frame passes, a single missed one
			// (2 periods) trips. Same ratio the 471 Hz configs used
			// when it was written in seconds.
			data->gap_trip = 1.5 * dt;
			log_info("fsp: auto gap_trip -> %.2f ms (1.5 frames at "
				"%.2f Hz)", 1e3 * data->gap_trip, 1.0 / dt);
		}
		if (data->delay_auto_bias_auto && data->delay_ident_auto) {
			// No identification to work from, so MODEL the sensor
			// path from what the source publishes: half an exposure
			// from the integration midpoint to its end, plus a
			// fixed pipeline term.
			//
			// THE 2.5 IS MEASURED, not assumed. Removing the
			// exposure half from this bench's four Bode
			// identifications leaves a remarkably constant
			// remainder across an 8x range of frame rate and a 20x
			// range of exposure:
			//   2026-07-31  3788 Hz  exp 50us   x 2.67  y 2.53 fr
			//   2026-08-07  2310 Hz  exp 400us  x 2.56  y 2.66 fr
			// i.e. the camera pipelines about two frames of readout
			// and transfer (integrate N while N-1 moves over USB,
			// plus the SDK's queue) on top of the half-frame DAC
			// hold. That the remainder is constant in FRAMES rather
			// than in milliseconds is itself the finding: the delay
			// here is dominated by frame-coupled terms.
			//
			// It is still a model. It cannot see actuator lag, and
			// it assumes this camera's pipelining. Run the PRBS
			// delay calibration and the suite replaces
			// delay_ident_ms with a real number.
			struct aylp_loop_timing t;
			double expo = aylp_timing_get(&t) ? t.exposure : 0.0;
			data->delay_auto_bias = 0.5 * expo * data->fs + 2.5;
			data->delay_ident_ms = 1e3 * (data->delay_auto_bias
				+ comp * data->fs) / data->fs;
			log_warn("fsp: delay_ident_ms is auto: MODELLING the "
				"sensor path as %.4G frames (%.3f ms "
				"exposure/2 + 2 frames of camera pipeline + 0.5 "
				"frame ZOH) instead of using an identification. "
				"The 2.5 comes from this bench's own Bode fits "
				"(2.53-2.67 frames of pipeline across 471-3788 "
				"Hz), but it cannot see actuator lag and it "
				"assumes this camera. Run the PRBS delay "
				"calibration for a real number.",
				data->delay_auto_bias, 0.5e3 * expo);
			if (!(expo > 0.0))
				log_warn("fsp: ... and the source published no "
					"exposure, so even the integration term "
					"is missing from that model.");
		} else if (data->delay_auto_bias_auto) {
			// The bias is a physical latency: recover it in SECONDS
			// from the identified delay, then convert at the rate the
			// loop actually achieved. This is what makes the whole
			// arrangement survive a frame-rate change.
			double bias_s = 1e-3 * data->delay_ident_ms - comp;
			if (bias_s < 0.0) {
				log_warn("fsp: auto bias: measured compute %.3f ms "
					"already exceeds the identified delay "
					"%.3f ms; clamping the bias to zero and "
					"trusting the compute term alone. "
					"Re-identify the delay.",
					1e3 * comp, data->delay_ident_ms);
				bias_s = 0.0;
			}
			data->delay_auto_bias = bias_s * data->fs;
			log_info("fsp: auto delay_auto_bias -> %.4G frames "
				"(%.3f ms identified - %.3f ms compute, at "
				"%.2f Hz)", data->delay_auto_bias,
				data->delay_ident_ms, 1e3 * comp, data->fs);
		}
		double comp_frames = comp * data->fs;
		double total = data->delay_auto_bias + comp_frames;
		size_t D = (size_t)total;
		double F = total - (double)D;
		if (D < 1) {
			// a delay below one sample is not representable here and
			// is not physical either: the frame in hand is already
			// one integration old
			D = 1;
			F = 0.0;
		}
		if (D > data->delay_auto_max) {
			log_warn("fsp: auto delay estimate %.3G frames exceeds "
				"delay_auto_max %zu, which is what every ring "
				"was sized for; clamping. Raise "
				"delay_auto_max if the estimate is believable.",
				total, data->delay_auto_max);
			D = data->delay_auto_max;
			F = 0.0;
		}
		if (!data->fs_auto && fabs(fs_meas - data->fs) > 0.1 * data->fs)
			log_warn("fsp: the loop is running at %.1f Hz, %.0f%% "
				"from the configured fs %G. fs sets the mode "
				"AR coefficients, broad_lp, the drift time "
				"constants and the frames this delay is "
				"counted in -- all of them are off by that "
				"factor until fs is corrected.",
				fs_meas, 1e2 * fabs(fs_meas - data->fs)
					/ data->fs, data->fs);
		for (int a = 0; a < 2; a++) {
			struct aylp_fsp_axis *ax = &data->axis[a];
			if (!ax->delay_auto) continue;
			ax->delay = D;
			ax->delay_frac = F;
			// both of these are functions of the delay
			fsp_build_comp(ax, data);
			if (data->guard_ratio > 0.0) fsp_build_guard(ax, data);
			// The command ring just changed length. Everything in
			// it is zero -- the hold parks the output there -- so
			// restarting it clean loses nothing and keeps uhead
			// consistent with the new modulus.
			memset(ax->ucmd, 0,
				(ax->delay_alloc + 1) * sizeof(double));
			ax->uhead = 0;
			ax->frac_x1 = 0.0;
			ax->frac_y1 = 0.0;
		}
		log_info("fsp: auto delay -> %zu+%.3G frames (%.3G total) = "
			"%.3G measured compute (%.3f ms downstream of the "
			"source, %.3f ms frame period, %zu frames / %zu timing "
			"windows) + %G assumed sensor/ZOH bias. THIS IS AN "
			"ESTIMATE: the bias is assumed, and the compute term "
			"is the whole downstream pipeline, so sinks running "
			"after the DAC inflate it. Confirm with a Bode sweep "
			"before trusting an attenuation number.",
			D, F, total, comp_frames, 1e3 * comp, 1e3 * dt,
			data->auto_n_dt, data->auto_n_c,
			data->delay_auto_bias);
	}
	xfree(data->auto_dt);
	xfree(data->auto_compute);
}


int fsp_proc(struct aylp_device *self, struct aylp_state *state)
{
	struct aylp_fsp_data *data = self->device_data;
	gsl_vector *s = state->vector;

	struct timespec tp;
	int err = clock_gettime(CLOCK_MONOTONIC, &tp);
	if (err) {
		log_error("fsp: couldn't get time: %s", strerror(err));
		return -1;
	}
	double now = tp.tv_sec + 1e-9 * tp.tv_nsec;
	if (UNLIKELY(!data->t0)) data->t0 = now;
	bool sensor_lost = state->header.status & AYLP_BEAM_LOST;
	// Both flags mean this iteration carries no trustworthy fresh centroid.
	// center_of_mass uses NO_SIGNAL while holding through a dark pulse;
	// FRAME_REJECTED denotes a distorted-but-present frame. Treat both as a
	// measurement hole so held coordinates cannot train the observers or be
	// integrated as if they were hundreds of identical real samples.
	bool sensor_rejected = state->header.status
		& (AYLP_FRAME_REJECTED | AYLP_NO_SIGNAL);
	if (UNLIKELY(sensor_lost && !data->beam_lost)) {
		data->beam_lost = true;
		data->beam_recover_start = 0.0;
		for (size_t a = 0; a < 2; a++)
			fsp_reset_after_beam_loss(&data->axis[a], data);
		log_error("fsp: BEAM LOST from upstream sensor; holding output at "
			"zero while acquisition searches");
	} else if (UNLIKELY(!sensor_lost && data->beam_lost)) {
		data->beam_lost = false;
		data->beam_recover_start = now;
		log_info("fsp: beam re-acquired; resuming normal operation over "
			"%G s", data->beam_recover_ramp);
	}

	// Convergence sample. Taken from the top of proc so it covers the
	// start_delay hold too -- a nonzero ||dw|| there would mean the
	// observer is learning from data it is supposed to be ignoring.
	if (UNLIKELY(data->wiener_trace_fp
			&& now - data->t_wtrace >= data->wiener_trace_period)) {
		data->t_wtrace = now;
		fsp_trace_sample(data, now);
	}

	if (UNLIKELY(!data->res_v || data->res_v->size != s->size)) {
		if (data->res_v) xfree_type(gsl_vector, data->res_v);
		data->res_v = xmalloc_type(gsl_vector, s->size);
	}
	gsl_vector *r = data->res_v;
	data->n_elem = s->size;

	// closed-loop authority: 0 during the startup hold, then a linear ramp
	// from 0 to 1 over `ramp` seconds so the handover is bumpless
	bool in_hold = (now - data->t0) < data->start_delay;
	// An auto delay is measured during the hold and installed as soon as
	// the measurement window closes -- NOT at the end of the hold. The
	// broadband filter trains against a target delay + broad_gd samples
	// ahead, so everything it learns before the delay is known is learned
	// at the wrong horizon; installing early leaves the rest of the hold
	// (usually all of it) training at the right one. The command ring the
	// install resizes holds nothing but zeros either way: the hold parks
	// the output there.
	if (UNLIKELY(data->delay_auto && !data->delay_auto_done)) {
		if (in_hold) {
			fsp_delay_auto_sample(data, now);
			if (now - data->t0 >= data->delay_auto_settle
					+ data->delay_auto_window)
				fsp_delay_auto_lock(data);
		} else {
			// the hold ended before the window did (start_delay was
			// raised at init to prevent this, so this is a loop
			// running far slower than its configured fs)
			fsp_delay_auto_lock(data);
		}
	}
	double precenter_frac = 0.0;
	if (in_hold && data->precenter_clamp > 0.0
			&& now - data->t0 >= data->precenter_delay) {
		double up = now - data->t0 - data->precenter_delay;
		precenter_frac = data->precenter_ramp > 0.0
			&& up < data->precenter_ramp ? up/data->precenter_ramp : 1.0;
	}
	double frac;
	if (in_hold) {
		frac = 0.0;
	} else {
		if (UNLIKELY(!data->t_close)) {
			data->t_close = now;
			data->t_adapt = now;
		}
		double up = now - data->t_close;
		frac = (data->ramp > 0.0 && up < data->ramp)
			? up / data->ramp : 1.0;
		// NLMS step schedule: decay from broad_mu_init to broad_mu
		// with time constant broad_mu_tau, measured from loop close.
		// Converges the slow eigenmodes during settle_time instead of
		// during the scored window; see fsp.h.
		if (data->broad_mu_init > data->broad_mu)
			data->broad_mu_cur = data->broad_mu
				+ (data->broad_mu_init - data->broad_mu)
				* exp(-up / data->broad_mu_tau);
		if (UNLIKELY(frac > 0.0 && !data->closed)) {
			data->closed = true;
			log_info("fsp: loop closing; blending command in "
				"over %G s", data->ramp);
		}
	}
	if (data->beam_lost) {
		frac = 0.0;
	} else if (data->beam_recover_start > 0.0) {
		double mix = data->beam_recover_ramp > 0.0
			? (now - data->beam_recover_start) / data->beam_recover_ramp
			: 1.0;
		if (mix >= 1.0) {
			mix = 1.0;
			data->beam_recover_start = 0.0;
			log_info("fsp: beam re-acquisition ramp complete");
		}
		frac *= mix;
	}

	// --- stall-gap patch: a proc-to-proc gap beyond gap_trip means the
	// source dropped frames (scheduler hiccup) or stalled outright
	// (camera stream stall + capture restart). The DAC held the last
	// command throughout, so patch each history with what actually
	// happened, sized to the frames missed, and continue -- no hold, no
	// ramp, no authority change (see fsp.h). ---
	if (UNLIKELY(data->gap_trip > 0.0 && data->t_last > 0.0
			&& now - data->t_last >= data->gap_trip)) {
		double gap_s = now - data->t_last;
		size_t miss = (size_t)(gap_s * data->fs + 0.5);
		miss = miss > 1 ? miss - 1 : 1;	// this frame isn't missed
		data->gap_events++;
		log_warn("fsp: %.1f ms frame gap (~%zu frames; source %s); "
			"padding histories with the held command and continuing "
			"(gap event %zu)", 1e3 * gap_s, miss,
			gap_s > 0.1 ? "stall" : "hiccup", data->gap_events);
		for (size_t a = 0; a < 2; a++) {
			struct aylp_fsp_axis *ax = &data->axis[a];
			size_t ring_len = ax->delay + 1;
			// The DAC held the last commanded u (r gets the
			// pre-frac-filter u): run the fractional all-pass and
			// the Smith ring forward on that constant input for
			// the missed frames (capped: past ring_len + settling
			// more iterations change nothing)
			double u_held = ax->frac_x1;
			double a_frac = (1.0 - ax->delay_frac)
				/ (1.0 + ax->delay_frac);
			size_t n = miss < ring_len + 16 ? miss : ring_len + 16;
			for (size_t k = 0; k < n; k++) {
				double u_frac = a_frac * u_held + ax->frac_x1
					- a_frac * ax->frac_y1;
				ax->frac_x1 = u_held;
				ax->frac_y1 = u_frac;
				ax->ucmd[ax->uhead] = u_frac;
				ax->uhead = (ax->uhead + 1) % ring_len;
			}
			// Propagate the modal state through the gap: the AR
			// recursion advances each mode's phase exactly as the
			// model says it should. It is a contraction, so this
			// is safe at any gap length; cap the work at 2 s of
			// frames, beyond which sharp lines have decayed.
			size_t np = miss < (size_t)(2.0 * data->fs)
				? miss : (size_t)(2.0 * data->fs);
			for (size_t k = 0; k < np; k++) {
				for (size_t i = 0; i < ax->n_modes; i++) {
					double p0 = ax->a1[i] * ax->xhat[2*i]
						+ ax->a2[i] * ax->xhat[2*i+1];
					ax->xhat[2*i+1] = ax->xhat[2*i];
					ax->xhat[2*i] = p0;
				}
			}
			// Pad the vibration history with the estimated AC residual.
			// The separately tracked drift state carries the DC correction
			// through the gap; phases across the hole are unknowable.
			// the first live frame after the gap is a step from
			// wherever the disturbance walked to; re-seed the burst
			// detector so the step itself cannot ring it
			ax->gd_precharge = true;
			if (data->broad_order) {
				size_t H = ax->broad_hist_len;
				size_t nh = miss < H ? miss : H;
				double phi_pad = ax->phi_dc - ax->drift_hat;
				for (size_t k = 0; k < nh; k++) {
					ax->broad_head = (ax->broad_head + 1)
						% H;
					ax->broad_hist[ax->broad_head] =
						phi_pad;
					ax->broad_hist_bl[ax->broad_head] =
						phi_pad;
					ax->broad_syn[ax->broad_head] = 1;
				}
				ax->broad_lpclean = 0;
				// pad the observer prefilter ring the same way
				size_t nl = miss < data->broad_lp
					? miss : data->broad_lp;
				for (size_t k = 0; k < nl; k++) {
					ax->broad_lpbuf[ax->broad_lphead] =
						phi_pad;
					ax->broad_lphead = (ax->broad_lphead
						+ 1) % data->broad_lp;
				}
				ax->broad_seen += nh;
				// pause training until the fabricated samples
				// leave the tap window -- unless they are a
				// negligible fraction of it (a few-frame drop
				// trains straight through)
				if (nh > 8) ax->broad_fab = H;
			}
		}
	}
	data->t_last = now;

	double beta = 0.0;	// EWMA weight for adaptation stats
	if (data->adapt_period > 0.0 && data->adapt_tau > 0.0)
		beta = 1.0 - exp(-1.0 / (data->adapt_tau * data->fs));

	for (size_t j = 0; j < s->size; j++) {
		// only elements 0 (y) and 1 (x) are controlled; pass any extras
		if (j > 1) { r->data[j * r->stride] = 0.0; continue; }
		struct aylp_fsp_axis *ax = &data->axis[j];
		double e = s->data[j * s->stride];
		ax->n_total++;

		// per-axis delay bookkeeping (each axis is visited exactly once
		// per frame, so advancing here keeps per-frame cadence)
		size_t delay = ax->delay;
		size_t ring_len = delay + 1;
		size_t slot_older = ax->uhead;	// u(k-delay-1)
		size_t slot = (slot_older + 1) % ring_len;	// u(k-delay)
		if (UNLIKELY(data->beam_lost)) {
			// center_of_mass intentionally holds its last coordinate while
			// searching. Do not treat that stale value as a measurement or
			// let it train either observer; just advance the real zero command.
			ax->ucmd[slot_older] = 0.0;
			ax->uhead = (ax->uhead + 1) % ring_len;
			r->data[j * r->stride] = 0.0;
			continue;
		}
		if (UNLIKELY(sensor_rejected)) {
			// Upstream deliberately held its last coordinate because this
			// frame was distorted, not because the beam disappeared. Do not
			// train observers or correct modal state with the duplicated
			// measurement. Time still advanced, so propagate the autonomous
			// modal state, pad broadband history with its slow baseline, and
			// use the direct bank command below when it is ready;
			// freezing those states would leave every oscillator one frame late.
			// Advance the Bode-shaped plant model on the delayed command even
			// though there is no measurement from which to form phi_meas.
			if (ax->plant_shaped)
				(void)fsp_biquad(ax->ucmd[slot], ax->plant_b, ax->plant_a,
					&ax->plant_z1, &ax->plant_z2);
			for (size_t i = 0; i < ax->n_modes; i++) {
				double p0 = ax->a1[i] * ax->xhat[2*i]
					+ ax->a2[i] * ax->xhat[2*i+1];
				ax->xhat[2*i+1] = ax->xhat[2*i];
				ax->xhat[2*i] = p0;
			}
			ax->n_dark++;
			ax->broad_dark_run++;
			bool dark_commandable = false;
			double dark_broad_hat = 0.0;
			if (data->broad_order) {
				size_t H = ax->broad_hist_len;
				// Keep the ordinary observer baseline-masked through the
				// hole. Dark prediction is used only for the command below,
				// so an estimate can never feed this or a later prediction.
				double phi_pad = ax->phi_dc - ax->drift_hat;
				// the masked history gets what masking-only
				// would have stored: the baseline, boxcar'd
				// like any other raw sample
				double bl_pad = phi_pad;
				if (data->broad_lp) {
					ax->broad_lpbuf[ax->broad_lphead] = phi_pad;
					ax->broad_lphead = (ax->broad_lphead + 1)
						% data->broad_lp;
					double sum = 0.0;
					for (size_t i = 0; i < data->broad_lp; i++)
						sum += ax->broad_lpbuf[i];
					bl_pad = sum / (double)data->broad_lp;
				}
				if (data->dark_predict && ax->dark_w
						&& !ax->broad_fab
						&& ax->broad_seen >= H
						&& ax->broad_dark_run
							<= ax->dark_predict_max) {
					size_t P = data->broad_order;
					size_t k = ax->broad_dark_run;
					// the slot about to be written is
					// now-time, so the window ending at
					// the last real sample starts k slots
					// behind it
					size_t idx = (ax->broad_head + 1 + H
						- k) % H;
					// Command from a single direct predictor: asking the
					// main FIR to predict from a synthetic current sample would
					// cascade two imperfect predictors.  Bank horizon
					// k+bd predicts the same future boxcar target as the
					// live main FIR, directly from the last real sample.
					size_t bd = ax->delay + data->broad_gd;
					size_t hc = k + bd;
					if (hc + 1 <= ax->dark_bank_max) {
						double pc = 0.0, pcn = 0.0;
						const double *wc = ax->dark_w
							+ (hc - 1) * P;
						const double *wcn = wc + P;
						idx = (ax->broad_head + 1 + H - k) % H;
						for (size_t i = 0; i < P; i++) {
							double v = ax->broad_hist_bl[idx];
							pc += wc[i] * v;
							pcn += wcn[i] * v;
							idx = idx ? idx - 1 : H - 1;
						}
						dark_broad_hat = (1.0 - ax->delay_frac) * pc
							+ ax->delay_frac * pcn;
						dark_commandable = isfinite(dark_broad_hat);
						if (dark_commandable) ax->n_dark_predicted++;
					}
				}
				// however good the stand-in, it is not a
				// measurement: mark it so no update ever trains
				// on it, and restart the clean-target count,
				// since the prefilter output is now a boxcar
				// over partly fabricated raw samples
				ax->broad_lpclean = 0;
				ax->broad_head = (ax->broad_head + 1) % H;
				ax->broad_hist[ax->broad_head] = bl_pad;
				ax->broad_hist_bl[ax->broad_head] = bl_pad;
				ax->broad_syn[ax->broad_head] = 1;
				ax->broad_seen++;
				// A hole long enough to fill a serious part of
				// the tap window leaves no regressor worth
				// training against, per-tap exclusion or not:
				// pause learning until it has walked out. This
				// is for an OUTAGE, and the threshold has to be
				// measured against the window rather than set
				// to some small number of frames -- a 13-frame
				// chop trips any small threshold on every
				// streak, and since broad_fab counts down one
				// frame at a time it would then never reach
				// zero and the filter would never train again.
				// Ordinary chops are handled tap by tap.
				size_t fab_at = data->dark_fab_frames
					? data->dark_fab_frames : H / 2;
				if (ax->broad_dark_run == fab_at)
					ax->broad_fab = H;
			}
			// Command through the hole from the direct long-horizon bank.
			// Anything non-nominal falls back to the old held command.
			double u_out = ax->frac_x1;
			if (dark_commandable && !ax->guard_active
					&& !ax->transient_active) {
				double broad_hat = dark_broad_hat;
				if (isfinite(broad_hat)) {
					double v = frac * (-(ax->drift_hat
						+ broad_hat) / ax->K);
					double u = v;
					if (ax->plant_shaped)
						u = fsp_biquad(v, ax->plant_ib,
							ax->plant_ia,
							&ax->plant_iz1,
							&ax->plant_iz2);
					if (data->cmd_fc > 0.0) {
						double uf = data->lp_b0*u
							+ ax->lp_z1;
						ax->lp_z1 = data->lp_b1*u
							- data->lp_a1*uf
							+ ax->lp_z2;
						ax->lp_z2 = data->lp_b2*u
							- data->lp_a2*uf;
						u = uf;
					}
					if (!isfinite(u)) {
						u = 0.0;
						ax->lp_z1 = ax->lp_z2 = 0.0;
						ax->plant_iz1 =
							ax->plant_iz2 = 0.0;
					}
					if (u > ax->clamp_hi) u = ax->clamp_hi;
					if (u < ax->clamp_lo) u = ax->clamp_lo;
					u_out = u;
					ax->n_dark_cmd++;
				}
			}
			double a_frac = (1.0 - ax->delay_frac)
				/ (1.0 + ax->delay_frac);
			double u_frac = a_frac*u_out + ax->frac_x1
				- a_frac*ax->frac_y1;
			ax->frac_x1 = u_out;
			ax->frac_y1 = u_frac;
			ax->ucmd[slot_older] = u_frac;
			ax->uhead = (ax->uhead + 1) % ring_len;
			r->data[j * r->stride] = u_out;
			continue;
		}
		bool dark_exit = ax->broad_dark_run != 0;
		ax->broad_dark_run = 0;
		if (data->broad_order) {
			ax->broad_head = (ax->broad_head + 1)
				% ax->broad_hist_len;
			ax->broad_seen++;
			// a real measurement lands here; phi_bl is written
			// further down, once the Smith reconstruction exists
			ax->broad_syn[ax->broad_head] = 0;
		}

		// Learn the ordinary open-loop operating point.  Use the same slow
		// timescale as the other statistics; attenuation_test provides 500 s
		// of genuine open-loop data, so the preceding zeroed START phase has
		// completely washed out before closing.
		if (in_hold) {
			double tau = data->adapt_tau > 0.0 ? data->adapt_tau : 5.0;
			double b = 1.0 - exp(-1.0 / (tau * data->fs));
			ax->trip_center += b * (e - ax->trip_center);
		}

		// --- burst guard: band-pass the raw error at this axis's
		// regeneration frequency, envelope-detect, and shed authority
		// while a burst is alive (see fsp.h). g_gain multiplies this
		// axis's closed-loop authority; suppress freezes learning so
		// the NLMS/adaptation never trains on the ringing. ---
		double g_gain = 1.0;
		bool suppress = false;
		bool over = false;
		if (data->guard_ratio > 0.0) {
			// First live frame after a run of dark frames or a
			// stall gap: the band-pass state still holds the
			// pre-hole signal, and hitting it with the post-hole
			// value is a step. A band-pass rung by a step at every
			// chop period looks exactly like sustained
			// regeneration (measured on a 25% chop: authority shed
			// to u_rms 0.0056, attenuation 1.91x -> 1.01x), so
			// pre-charge the filter to steady state at this
			// frame's value -- for a constant-peak band-pass
			// (b1 = 0, b2 = -b0) that is z1 = z2 = -b0*e, output
			// exactly 0 -- and let the detector respond only to
			// motion WITHIN contiguous live data.
			if (UNLIKELY(dark_exit || ax->gd_precharge)) {
				ax->gd_precharge = false;
				ax->gd_z1 = ax->gd_z2 = -ax->gd_b0 * e;
			}
			double bp = ax->gd_b0 * e + ax->gd_z1;
			ax->gd_z1 = -ax->gd_a1 * bp + ax->gd_z2;
			ax->gd_z2 = ax->gd_b2 * e - ax->gd_a2 * bp;
			if (!isfinite(bp)) {
				bp = 0.0;
				ax->gd_z1 = ax->gd_z2 = 0.0;
			}
			ax->gd_env += data->guard_beta_fast
				* (bp * bp - ax->gd_env);
			double base = ax->gd_base;
			double fl = data->guard_floor * data->guard_floor;
			if (base < fl) base = fl;
			double thr = data->guard_ratio * data->guard_ratio
				* base;
			over = ax->gd_env > thr;
			if (!ax->guard_active && frac > 0.0 && over) {
				// sustained-energy debounce: require the
				// envelope to stay over threshold for
				// gd_min_samples consecutive samples before
				// latching. A genuine regeneration ring stays
				// over for many cycles; a single spike from
				// low-frequency content leaking through the
				// band-pass skirts (the observed false-trip
				// mode) does not.
				ax->gd_over_count++;
				size_t need = ax->gd_min_samples > 0
					? ax->gd_min_samples : 1;
				if (ax->gd_over_count >= need) {
					size_t sustained = ax->gd_over_count;
					ax->guard_active = true;
					ax->guard_ramping = false;
					ax->gd_over_count = 0;
					ax->guard_t_trip = now;
					ax->guard_events++;
					data->guard_events++;
					if (now - ax->guard_t_log > 2.0) {
						ax->guard_t_log = now;
						log_warn("fsp: burst guard tripped on "
							"%s axis (event %zu): %.0f Hz "
							"envelope %.3G > %.3G sustained "
							"%zu samples; holding authority "
							"at 0 for %G s, then ramping "
							"back over %G s",
							j == 0 ? "y" : "x",
							ax->guard_events, ax->gd_f0,
							sqrt(ax->gd_env), sqrt(thr),
							sustained,
							data->guard_hold,
							data->guard_ramp);
					}
				}
			} else if (!ax->guard_active && frac > 0.0 && !over) {
				ax->gd_over_count = 0;
				// learn the quiet baseline only outside bursts
				ax->gd_base += data->guard_beta_slow
					* (ax->gd_env - ax->gd_base);
			}
		}
		// hold/ramp re-entry, applied outside the detector gate so a
		// stall-gap trip re-enters the same way even when the band-pass
		// detector is disabled (guard_ratio <= 0)
		if (ax->guard_active) {
			double ts = now - ax->guard_t_trip;
			if (ts <= data->guard_hold) {
				g_gain = 0.0;
			} else if (over) {
				// still ringing after the hold: back to
				// zero authority. If recovery had begun
				// this is a fresh regeneration (a new
				// event); a ring that merely outlives
				// the hold extends the same event.
				if (ax->guard_ramping) {
					ax->guard_events++;
					data->guard_events++;
					ax->guard_ramping = false;
				}
				ax->guard_t_trip = now;
				g_gain = 0.0;
			} else {
				double rmp = data->guard_ramp > 0.0
					? data->guard_ramp : 1e-9;
				ax->guard_ramping = true;
				g_gain = (ts - data->guard_hold) / rmp;
				if (g_gain >= 1.0) {
					g_gain = 1.0;
					ax->guard_active = false;
				}
			}
			if (g_gain < 1.0) {
				suppress = true;
				ax->guard_frames++;
			}
		}
		// stall-gap bookkeeping: count down the frames until the
		// fabricated (padded) samples have left the NLMS tap window
		if (UNLIKELY(ax->broad_fab)) {
			ax->broad_fab--;
			suppress = true;
		}

		// --- Smith-predictor core: reconstruct the disturbance by
		// removing our own delayed plant contribution ---
		// Integer delay after a first-order Thiran all-pass fractional
		// delay. Unlike linear interpolation, the Thiran model has unity
		// magnitude through Nyquist and therefore does not leave a false
		// high-frequency command residue in the Smith reconstruction.
		double u_old = ax->ucmd[slot];
		// The Bode-shaped model is evaluated on the command that actually
		// entered the delay line (after clamp and fractional-delay shaping).
		// Since H and the transport delay are both LTI they commute, so this
		// is K H(z) D(z) u with only two extra states per axis.
		double plant_u = u_old;
		if (ax->plant_shaped)
			plant_u = fsp_biquad(u_old, ax->plant_b, ax->plant_a,
				&ax->plant_z1, &ax->plant_z2);
		double phi_meas = e - ax->K * plant_u;
		// slow DC estimate of the reconstructed disturbance, used to
		// pad the NLMS history across a stall gap
		ax->phi_dc += data->gap_dc_beta * (phi_meas - ax->phi_dc);
		// Optional two-timescale model: cancel drift as a separate state and
		// train the vibration predictor only on the zero-mean residual.
		if (data->drift_tau > 0.0) {
			if (data->drift_order == 2) {
				// Alpha-beta realization of a constant-velocity disturbance
				// model.  The position pole supplies native integral/DC
				// rejection; the rate pole removes the following error on a
				// ramp.  beta=alpha^2/2 gives the rate state the same slow
				// scale without making frame noise look like velocity.
				double predicted = ax->drift_hat
					+ ax->drift_rate / data->fs;
				double innovation = phi_meas - predicted;
				ax->drift_hat = predicted
					+ data->drift_beta * innovation;
				ax->drift_rate += 0.5 * data->drift_beta
					* data->drift_beta * data->fs * innovation;
			} else {
				ax->drift_hat += data->drift_beta
					* (phi_meas - ax->drift_hat);
			}
		}
		double phi_vib = phi_meas - ax->drift_hat;

		// Full compound-disturbance observer (Kulcsar/Petit/Meimon):
		// identify the delay-step conditional mean of the reconstructed
		// disturbance. In scalar form the stationary Kalman/Wiener observer
		// is an FIR state-space realization; NLMS tracks slow spectrum drift
		// without a large online covariance matrix.
		double broad_hat = 0.0;
		double transient_mix = 0.0;
		if (data->broad_order) {
			size_t P = data->broad_order, H = ax->broad_hist_len;
			// Observer band-limit (see fsp.h): boxcar the raw phi so
			// the NLMS never sees (nor learns to chase) the HF command
			// echo; its exact integer group delay broad_gd is added to
			// the prediction horizon below, so in-band timing holds.
			double phi_bl = phi_vib;
			if (data->broad_lp) {
				ax->broad_lpbuf[ax->broad_lphead] = phi_vib;
				ax->broad_lphead = (ax->broad_lphead + 1)
					% data->broad_lp;
				double sum = 0.0;
				for (size_t i = 0; i < data->broad_lp; i++)
					sum += ax->broad_lpbuf[i];
				phi_bl = sum / (double)data->broad_lp;
			}
			// phi_bl is the boxcar of the last broad_lp raw
			// samples, so it is only a clean training target once
			// that many real ones have been pushed in a row
			if (ax->broad_lpclean < data->broad_lp + 1)
				ax->broad_lpclean++;
			bool target_clean = !data->broad_lp
				|| ax->broad_lpclean >= data->broad_lp;
			size_t bd = delay + data->broad_gd;
			ax->broad_hist[ax->broad_head] = phi_bl;
			ax->broad_hist_bl[ax->broad_head] = phi_bl;
			if (ax->broad_seen >= H) {
				// The delay-tap and (delay+1)-tap regressors are the
				// same P-sample window shifted by one, sharing P-1
				// taps: read the union (P+1 samples) once instead of
				// walking the ring twice. xbuf/pred1/energy1 belong to
				// the delay predictor (broad_w); xbuf2/pred2/energy2 to
				// the delay+1 predictor (broad_w_next).
				// Provenance rides along with the samples. Every
				// tap feeds the PREDICTION -- a synthesized
				// sample is the best estimate of an instant
				// nobody measured, and dropping it would leave
				// a hole in the regressor instead. Only the
				// UPDATE excludes them, so the filter is never
				// trained toward its own output; energy*_upd is
				// the norm over just the taps that will move,
				// which is what an NLMS step must normalize by.
				size_t idx = (ax->broad_head + H - bd) % H;
				double pred1 = 0.0, energy1 = 1e-12;
				double pred2 = 0.0, energy2 = 1e-12;
				double pred1_upd = 0.0, pred2_upd = 0.0;
				double energy1_upd = 1e-12, energy2_upd = 1e-12;
				size_t n_syn1 = 0, n_syn2 = 0;
				unsigned char sy;
				double v = ax->broad_hist[idx];	// walk1-only tap
				sy = ax->broad_syn[idx];
				ax->broad_xbuf[0] = v;
				ax->broad_synbuf[0] = sy;
				pred1 += ax->broad_w[0] * v;
				energy1 += v * v;
				if (sy) n_syn1++; else {
					energy1_upd += v * v;
					pred1_upd += ax->broad_w[0] * v;
				}
				idx = idx ? idx - 1 : H - 1;
				for (size_t i = 1; i < P; i++) {
					v = ax->broad_hist[idx];
					sy = ax->broad_syn[idx];
					ax->broad_xbuf[i] = v;
					ax->broad_synbuf[i] = sy;
					pred1 += ax->broad_w[i] * v;
					energy1 += v * v;
					ax->broad_xbuf2[i - 1] = v;
					ax->broad_synbuf2[i - 1] = sy;
					pred2 += ax->broad_w_next[i - 1] * v;
					energy2 += v * v;
					if (sy) {
						n_syn1++;
						n_syn2++;
					} else {
						energy1_upd += v * v;
						energy2_upd += v * v;
						pred1_upd += ax->broad_w[i] * v;
						pred2_upd += ax->broad_w_next[i - 1] * v;
					}
					idx = idx ? idx - 1 : H - 1;
				}
				v = ax->broad_hist[idx];		// walk2-only tap
				sy = ax->broad_syn[idx];
				ax->broad_xbuf2[P - 1] = v;
				ax->broad_synbuf2[P - 1] = sy;
				pred2 += ax->broad_w_next[P - 1] * v;
				energy2 += v * v;
				if (sy) n_syn2++; else {
					energy2_upd += v * v;
					pred2_upd += ax->broad_w_next[P - 1] * v;
				}
				ax->n_syn_taps += n_syn1;
				// A partial-update NLMS step normalizes by the
				// norm of the taps it actually moves, so if only
				// a few real taps survive, that norm is tiny and
				// the step is enormous: one small sample would
				// take the whole correction. Refuse to train
				// until enough of the regressor is real. This is
				// the fine-grained version of broad_fab, and the
				// criterion has to be the COUNT rather than the
				// energy -- fabricated samples sit near the DC
				// baseline, so they carry almost no energy and
				// an energy test would call a window of them
				// clean.
				size_t need_real = (size_t)(data->dark_train_min_real
					* (double)P);
				bool enough1 = P - n_syn1 >= need_real;
				bool enough2 = P - n_syn2 >= need_real;

				double pe = phi_bl - pred1;
				transient_mix = fsp_transient_update(ax, data, pe,
					frac * g_gain, now, j);
				// A push is deliberately outside the stationary disturbance
				// distribution. Freeze both FIR horizons for the entire event,
				// including the return cross-fade, then resume on the first fully
				// normal frame. The frozen predictor remains the independent model
				// used by the release detector, so the event cannot teach itself
				// into agreement.
				bool train = (!data->broad_freeze_closed || in_hold)
					&& !suppress && !ax->transient_active
					// the target is a boxcar over the raw
					// samples; if any of those were
					// fabricated, this frame has no honest
					// target to train against
					&& target_clean;
				if (!target_clean) ax->n_train_skip_target++;
				double pe_upd = phi_bl - pred1_upd;
				if (train && enough1 && isfinite(pe_upd)) {
					double step = data->broad_mu_cur * pe_upd / energy1_upd;
					for (size_t i = 0; i < P; i++) {
						if (ax->broad_synbuf[i]) continue;
						ax->broad_w[i] += step * ax->broad_xbuf[i];
					}
				} else if (!isfinite(pe)) {
					memset(ax->broad_w, 0, P * sizeof(double));
				}
				// Identify the adjacent delay+1 predictor. Their weighted
				// combination is the conditional mean at the fractional
				// Bode-fit horizon.
				double pe2 = phi_bl - pred2;
				double pe2_upd = phi_bl - pred2_upd;
				fsp_shadow_update(ax, data, pe_upd, pe2_upd,
					pred1_upd, pred2_upd,
					energy1_upd, energy2_upd,
					suppress || !target_clean
						|| !enough1 || !enough2,
					now, j);
				if (train && enough2 && isfinite(pe2_upd)) {
					double step = data->broad_mu_cur * pe2_upd / energy2_upd;
					for (size_t i = 0; i < P; i++) {
						if (ax->broad_synbuf2[i]) continue;
						ax->broad_w_next[i] += step
							* ax->broad_xbuf2[i];
					}
				} else if (!isfinite(pe2)) {
					memset(ax->broad_w_next, 0, P * sizeof(double));
				}
				broad_hat = fsp_broad_predict(ax, data);
				if (!isfinite(broad_hat)) broad_hat = 0.0;
				// Train the dark bank: horizon-h filter learns
				// to predict this frame's clean target from
				// the baseline-masked history h frames back.
				// Full-update NLMS -- the baseline fills are
				// part of the bank's input signal by design,
				// the same signal its predictions read, so
				// there is nothing to exclude and no
				// partial-update normalization trap.
				// Round-robin dark_bank_train horizons per
				// frame bounds the added cost; broad_xbuf is
				// dead as scratch by this point.
				if (ax->dark_w && train) {
					for (size_t n = 0;
							n < data->dark_bank_train;
							n++) {
						size_t h = ax->dark_bank_h + 1;
						ax->dark_bank_h =
							(ax->dark_bank_h + 1)
							% ax->dark_bank_max;
						double *w = ax->dark_w
							+ (h - 1) * P;
						size_t di = (ax->broad_head
							+ H - h) % H;
						double dpred = 0.0;
						double denergy = 1e-12;
						for (size_t i = 0; i < P; i++) {
							double dv =
								ax->broad_hist_bl[di];
							ax->broad_xbuf[i] = dv;
							dpred += w[i] * dv;
							denergy += dv * dv;
							di = di ? di - 1 : H - 1;
						}
						double dpe = phi_bl - dpred;
						if (isfinite(dpe)) {
							double step =
								data->broad_mu_cur
								* dpe / denergy;
							for (size_t i = 0; i < P; i++)
								w[i] += step
									* ax->broad_xbuf[i];
						} else {
							memset(w, 0, P
								* sizeof(double));
						}
					}
				}
			}
		}

		// --- Kalman filter: predict one step, then correct ---
		// xpred = A xhat (per-mode AR recursion); Cxpred = sum positions
		double xpred[AYLP_FSP_MAX_DIM];
		double Cxpred = 0.0;
		for (size_t i = 0; i < ax->n_modes; i++) {
			double p0 = ax->a1[i]*ax->xhat[2*i]
				+ ax->a2[i]*ax->xhat[2*i+1];
			double p1 = ax->xhat[2*i];
			xpred[2*i] = p0;
			xpred[2*i+1] = p1;
			Cxpred += p0;
		}
		double innov = phi_vib - Cxpred;
		if (!isfinite(innov)) {
			// bad sample: reset the estimate, output 0; still advance
			// this axis's command ring to keep the delay line in step
			memset(ax->xhat, 0, ax->dim * sizeof(double));
			ax->plant_iz1 = ax->plant_iz2 = 0.0;
			ax->ucmd[slot_older] = 0.0;
			ax->uhead = (ax->uhead + 1) % ring_len;
			r->data[j * r->stride] = 0.0;
			continue;
		}
		if (!data->broad_order) {
			transient_mix = fsp_transient_update(ax, data, innov,
				frac * g_gain, now, j);
		}
		// During a push, use the high-Q physical-model gain to identify the
		// newly excited modes from successive position measurements.  Their
		// two-state AR realizations estimate both displacement and velocity;
		// normal operation keeps the quieter stationary gain.
		const double *observer_L = ax->L;
		if (ax->transient_active && ax->transient_modal)
			observer_L = ax->transient_L;
		for (size_t d = 0; d < ax->dim; d++)
			ax->xhat[d] = xpred[d] + observer_L[d] * innov;

		// --- adaptation statistics (cheap, per sample); frozen while
		// the burst guard is active so the ringing never contaminates
		// the identified model ---
		if (beta > 0.0 && !suppress && !ax->transient_active) {
			ax->r_ewma += beta * (innov*innov - ax->r_ewma);
			for (size_t i = 0; i < ax->n_modes; i++) {
				// E[x_i^2] = Var(xhat_i) + post_var_i, not
				// xhat_i^2 alone -- omitting post_var
				// understates true modal energy whenever the
				// posterior variance isn't already negligible,
				// biasing q low and, via the Riccati solve, the
				// gain low right along with it (see post_var's
				// definition in fsp.h).
				double amp = ax->xhat[2*i]*ax->xhat[2*i]
					+ ax->post_var[i];
				ax->q_ewma[i] += beta * (amp - ax->q_ewma[i]);
				// quadrature demod of the reconstructed
				// disturbance at the nominal line, to sense
				// frequency drift. Uses demod_beta, not beta:
				// the demodulator is a first-order PLL whose
				// capture range is ~1/(2*pi*time constant), so
				// it needs a time constant sized to
				// adapt_df_max (see demod_beta in fsp.h), not
				// to the (typically much slower) adapt_tau
				// used for r_ewma/q_ewma above.
				ax->demod_ph[i] += 2.0*M_PI*ax->f[i]/data->fs;
				if (ax->demod_ph[i] > M_PI)
					ax->demod_ph[i] -= 2.0*M_PI;
				ax->demod_re[i] += data->demod_beta * (phi_vib
					* cos(ax->demod_ph[i]) - ax->demod_re[i]);
				ax->demod_im[i] += data->demod_beta * (phi_vib
					* sin(ax->demod_ph[i]) - ax->demod_im[i]);
			}
		}

		// --- delay-step prediction: run the AR mean forward, harvesting
		// mode i's (boosted) contribution at its own real-valued
		// horizon delay + delay_frac + n_i (comp_n[i]/comp_frac[i],
		// its floor/fractional part -- see fsp_build_comp), so the
		// command filter's phase lag at that line nets out to ~0; the
		// minimum-variance command cancels the predicted sum ---
		double p_now[AYLP_FSP_MAX_DIM];
		memcpy(p_now, ax->xhat, ax->dim * sizeof(double));
		double phi_hat = 0.0;
		for (size_t step = 1; step <= ax->max_steps; step++) {
			for (size_t i = 0; i < ax->n_modes; i++) {
				double p0 = ax->a1[i]*p_now[2*i]
					+ ax->a2[i]*p_now[2*i+1];
				p_now[2*i+1] = p_now[2*i];
				p_now[2*i] = p0;
				size_t target = ax->comp_n[i];
				double tf = ax->comp_frac[i];
				if (step == target)
					phi_hat += ax->comp_g[i]
						* (1.0 - tf) * p0;
				else if (tf > 0.0 && step == target + 1)
					phi_hat += ax->comp_g[i] * tf * p0;
			}
		}

		// The full-band predictor and modal predictor estimate the same
		// disturbance. Select, do not sum, to avoid double cancellation.
		double vibration_hat = data->broad_order ? broad_hat : phi_hat;
		double drift_predictive = ax->drift_hat;
		if (data->drift_tau > 0.0 && data->drift_order == 2)
			drift_predictive += (ax->delay + ax->delay_frac)
				* ax->drift_rate / data->fs;
		double predictive_hat = drift_predictive + vibration_hat;
		// Normal operation cancels the separately estimated drift plus the
		// learned full-band prediction. Event recovery is selected independently.
		// It deliberately uses this same axis's Bode-derived K and the modal
		// prediction above uses this same axis's delay/delay_frac horizon; there
		// is no second transient plant gain or delay to keep synchronized.
		// P, I, PI, fast modal+P, or fast modal+PI (hybrid).
		double v_predictive = -predictive_hat / ax->K;
		double transient_vibration = ax->transient_modal
			? phi_hat : 0.0;
		bool transient_integral = ax->transient_active
			&& fsp_transient_uses_integral(ax->transient_controller);
		bool transient_proportional = ax->transient_controller
			!= AYLP_FSP_TRANSIENT_INTEGRAL;
		ax->transient_i_prev = ax->transient_i;
		ax->transient_i_step = 0.0;
		if (transient_integral && !ax->transient_recovering) {
			double next = ax->transient_i + (-data->transient_ki * e / ax->K
				- data->transient_i_leak * ax->transient_i) / data->fs;
			if (next > data->transient_i_limit)
				next = data->transient_i_limit;
			if (next < -data->transient_i_limit)
				next = -data->transient_i_limit;
			ax->transient_i_step = next - ax->transient_i;
			ax->transient_i = next;
		}
		double v_transient =
			-(ax->drift_hat + transient_vibration
				+ (transient_proportional
					? data->transient_kp * e : 0.0)) / ax->K
			+ (transient_integral ? ax->transient_i : 0.0);
		// During the startup hold, optionally move only the learned DC bias.
		// The modal/broadband prediction remains at zero authority until close,
		// so vibration control is trained around the centered operating point.
		double v;
		if (precenter_frac > 0.0) {
			v = -precenter_frac * ax->drift_hat / ax->K;
		} else if (data->precenter_clamp > 0.0 && !in_hold
				&& !ax->transient_active) {
			// Keep the learned centering bias continuously applied while
			// only the vibration term ramps in at loop close.
			v = -drift_predictive / ax->K
				+ frac * g_gain * (-vibration_hat / ax->K);
		} else {
			v = frac * g_gain
				* ((1.0 - transient_mix) * v_predictive
					+ transient_mix * v_transient);
		}
		double u = v;
		if (ax->plant_shaped)
			u = fsp_biquad(v, ax->plant_ib, ax->plant_ia,
				&ax->plant_iz1, &ax->plant_iz2);
		// A static centroid offset is normal and may require substantial DC
		// command to remove. These thresholds diagnose an implausible request;
		// the clamps contain it, while only explicit upstream beam loss latches
		// the output to zero.
		// A detected push is expected to violate the normal pointing-error
		// threshold. During that explicitly bounded window, use its separate
		// command threshold and hard timeout instead.
		bool over_error = !ax->transient_active && data->trip_error > 0.0
			&& fabs(e) > fabs(ax->trip_center) + data->trip_error;
		double command_trip = ax->transient_active
			? data->transient_trip_command : data->trip_command;
		bool over_command = command_trip > 0.0 && fabs(u) > command_trip;
		bool event_timeout = ax->transient_active
			&& data->transient_max_duration > 0.0
			&& now - ax->transient_start > data->transient_max_duration;
		if (frac > 0.0 && (over_error || over_command)) {
			ax->trip_count++;
		} else {
			ax->trip_count = 0;
			ax->trip_warned = false;
		}
		if (!ax->trip_warned && ax->trip_count >= data->trip_frames) {
			ax->trip_warned = true;
			log_warn("fsp: command/error limit on %s axis: e=%G "
				"(open baseline=%G), requested u=%G; output remains "
				"clamped", j == 0 ? "y" : "x", e,
				ax->trip_center, u);
		}
		if (event_timeout) {
			ax->transient_active = false;
			ax->transient_recovering = false;
			ax->shadow_active = false;
			ax->transient_blocked = true;
			ax->transient_i = 0.0;
			ax->transient_i_step = 0.0;
			log_warn("fsp: %s transient exceeded %G s; returning to "
				"normal authority until innovation is quiet",
				j == 0 ? "y" : "x", data->transient_max_duration);
		}
		if (data->beam_lost) u = 0.0;
		// command robustness low-pass (DF2T biquad); see fsp.h
		if (data->cmd_fc > 0.0) {
			double uf = data->lp_b0*u + ax->lp_z1;
			ax->lp_z1 = data->lp_b1*u - data->lp_a1*uf + ax->lp_z2;
			ax->lp_z2 = data->lp_b2*u - data->lp_a2*uf;
			u = uf;
		}
		if (data->beam_lost) {
			u = 0.0;
			ax->lp_z1 = ax->lp_z2 = 0.0;
			ax->plant_iz1 = ax->plant_iz2 = 0.0;
		}
		if (!isfinite(u)) {
			u = 0.0;
			ax->lp_z1 = ax->lp_z2 = 0.0;
			ax->plant_iz1 = ax->plant_iz2 = 0.0;
		}
		// Per-axis, asymmetric-capable output limit. This clamped value
		// is BOTH what leaves the device and what enters the delay ring
		// below, so the plant model and the actuator see the same
		// command -- which is why the limit belongs here rather than in
		// a downstream clamp stage. Per-axis because the two axes can
		// have opposite command->voltage signs, which mirrors an
		// asymmetric window (see aylp_fsp_axis).
		double clamp_hi = ax->transient_active
			? data->transient_clamp_hi : ax->clamp_hi;
		double clamp_lo = ax->transient_active
			? data->transient_clamp_lo : ax->clamp_lo;
		if (precenter_frac > 0.0) {
			clamp_hi = data->precenter_clamp;
			clamp_lo = -data->precenter_clamp;
		} else if (data->precenter_clamp > 0.0 && !in_hold
				&& !ax->transient_active) {
			// The normal clamp is FAST authority relative to the learned DC
			// operating point, not an absolute actuator bound. Otherwise loop
			// handoff would throw away the very bias pre-centering established.
			double dc = -ax->drift_hat / ax->K;
			if (dc > data->precenter_clamp) dc = data->precenter_clamp;
			if (dc < -data->precenter_clamp) dc = -data->precenter_clamp;
			clamp_hi = dc + ax->clamp_hi;
			clamp_lo = dc + ax->clamp_lo;
		}
		double u_unclamped = u;
		if (u > clamp_hi) u = clamp_hi;
		if (u < clamp_lo) u = clamp_lo;
		if (ax->transient_active && u != u_unclamped)
			ax->transient_saturated = true;
		// Conditional-integration anti-windup. The internal I limit is the
		// first boundary; if plant/filter dynamics still make the actuator
		// saturate and this frame's I step pushes farther into saturation,
		// roll that step back. The already bounded output remains continuous.
		if (transient_integral && u != u_unclamped
				&& (u_unclamped - u) * ax->transient_i_step > 0.0) {
			ax->transient_i = ax->transient_i_prev;
			ax->transient_i_step = 0.0;
		}
		if (ax->transient_active) {
			double ae = fabs(e), au = fabs(u);
			if (ae > data->transient_settle_error)
				ax->transient_error_t_last = now;
			if (ae > ax->transient_peak_error)
				ax->transient_peak_error = ae;
			if (au > ax->transient_peak_command)
				ax->transient_peak_command = au;
			if (fabs(ax->transient_i) > ax->transient_peak_i)
				ax->transient_peak_i = fabs(ax->transient_i);
			ax->transient_error2 += e * e;
			ax->transient_frames++;
		}

		// Apply the fractional part of the plant delay before the integer
		// command ring. H(z)=(a+z^-1)/(1+a*z^-1), with DC group delay f.
		double a_frac = (1.0 - ax->delay_frac)
			/ (1.0 + ax->delay_frac);
		double u_frac = a_frac*u + ax->frac_x1 - a_frac*ax->frac_y1;
		ax->frac_x1 = u;
		ax->frac_y1 = u_frac;
		ax->ucmd[slot_older] = u_frac;
		ax->uhead = (ax->uhead + 1) % ring_len;
		r->data[j * r->stride] = u;
	}

	data->n_seen += 1;
	if (data->t_close) data->n_closed += 1;

	// burst-guard ticker: a periodic terminal line with the activation
	// count and the fraction of closed-loop frames spent at reduced
	// authority, so a marginal run is visible while it happens
	if (data->guard_ratio > 0.0 && data->guard_tick > 0.0
			&& data->t_close) {
		if (UNLIKELY(!data->t_tick)) data->t_tick = data->t_close;
		if (now - data->t_tick >= data->guard_tick) {
			data->t_tick = now;
			size_t gy = data->axis[0].guard_events;
			size_t gx = data->axis[1].guard_events;
			double held = 100.0 * (double)(data->axis[0].guard_frames
				+ data->axis[1].guard_frames)
				/ (2.0 * (double)data->n_closed);
			log_info("fsp guard: %zu activations (y %zu / x %zu), "
				"authority shed %.2f%% of the closed loop so far",
				data->guard_events, gy, gx, held);
		}
	}

	// slow re-identification, amortized: start a new solve every
	// adapt_period, then advance any in-progress solve a few iterations per
	// frame so the Riccati never bursts (a single-frame solve stalled frame
	// delivery ~7 ms and logged as a periodic "source hiccup" -- see fsp.h).
	// Don't start a new one while the previous is still converging.
	if (data->adapt_period > 0.0 && data->t_close
			&& !data->axis[0].adapt_solving
			&& !data->axis[1].adapt_solving
			&& !data->axis[0].transient_active
			&& !data->axis[1].transient_active
			&& (now - data->t_adapt) >= data->adapt_period) {
		fsp_adapt(data);
		data->t_adapt = now;
	}
	for (int a = 0; a < 2; a++) {
		struct aylp_fsp_axis *ax = &data->axis[a];
		if (!ax->adapt_solving) continue;
		fsp_adapt_advance(ax, data, (size_t)a, AYLP_FSP_ADAPT_ITERS);
	}

	state->vector = r;
	return 0;
}


int fsp_fini(struct aylp_device *self)
{
	struct aylp_fsp_data *data = self->device_data;
	if (data->guard_ratio > 0.0 && data->n_closed)
		log_info("fsp guard final: %zu activations (y %zu / x %zu) "
			"over %zu closed-loop frames; recurring activations "
			"mean K/delay no longer match the plant -- re-run the "
			"bodes at the current coarse bias before trusting an "
			"attenuation number", data->guard_events,
			data->axis[0].guard_events, data->axis[1].guard_events,
			data->n_closed);
	if (data->gap_events)
		log_warn("fsp: %zu frame gap(s) patched this run (source "
			"stalled or hiccuped; check asi_source recovery lines)",
			data->gap_events);
	if (data->transient_sigma > 0.0)
		log_info("fsp transient final: %zu activations (y %zu / x %zu), "
			"%zu validated shadow promotion(s)",
			data->axis[0].transient_events
				+ data->axis[1].transient_events,
			data->axis[0].transient_events,
			data->axis[1].transient_events,
			data->shadow_promotions);
	if (data->transient_log_fp) {
		fclose(data->transient_log_fp);
		data->transient_log_fp = NULL;
		log_info("fsp: push-event summary written to %s",
			data->transient_log);
	}
	// a run killed during the startup hold still holds the auto-delay
	// sample buffers
	if (data->auto_dt) xfree(data->auto_dt);
	if (data->auto_compute) xfree(data->auto_compute);
	// must precede the frees below: this is the only chance to persist the
	// learned taps, and fini runs on the SIGINT path too
	fsp_save_wiener(data);
	if (data->wiener_trace_fp) {
		// final sample, so the trace ends at the same state the
		// wiener_out dump records
		fsp_trace_sample(data, data->t_last > 0.0 ? data->t_last
			: data->t0);
		fclose(data->wiener_trace_fp);
		data->wiener_trace_fp = NULL;
		log_info("fsp: observer convergence trace written to %s",
			data->wiener_trace);
	}
	for (int a = 0; a < 2; a++) {
		struct aylp_fsp_axis *ax = &data->axis[a];
		if (ax->n_dark)
			log_info("fsp: %s: %zu frame(s) the sensor could not "
				"measure (%.1f%% of %zu), %zu received a direct "
				"dark-bank prediction (%zu commanded through, the rest "
				"held) and the rest by the DC baseline; "
				"%zu frame(s) skipped training on a "
				"partly-fabricated target, %zu fabricated tap "
				"update(s) excluded",
				a == 0 ? "y" : "x", ax->n_dark,
				ax->n_total ? 1e2 * (double)ax->n_dark
					/ (double)ax->n_total : 0.0,
				ax->n_total, ax->n_dark_predicted,
				ax->n_dark_cmd,
				ax->n_train_skip_target, ax->n_syn_taps);
		xfree(data->wtrace_prev[a]);
		xfree(data->axis[a].ucmd);
		xfree(data->axis[a].broad_syn);
		xfree(data->axis[a].broad_synbuf);
		xfree(data->axis[a].broad_synbuf2);
		xfree(data->axis[a].broad_hist);
		xfree(data->axis[a].broad_hist_bl);
		xfree(data->axis[a].dark_w);
		xfree(data->axis[a].broad_w);
		xfree(data->axis[a].broad_w_next);
		xfree(data->axis[a].shadow_w);
		xfree(data->axis[a].shadow_w_next);
		xfree(data->axis[a].broad_xbuf);
		xfree(data->axis[a].broad_xbuf2);
		xfree(data->axis[a].broad_lpbuf);
	}
	if (data->res_v) xfree_type(gsl_vector, data->res_v);
	xfree(data->wiener_file);
	xfree(data->transient_log);
	xfree(data->wiener_out);
	xfree(data->wiener_trace);
	xfree(data);
	return 0;
}
