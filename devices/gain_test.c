// DC gain (and linearity) of the plant, measured with a staircase.
//
// Sits where the controller normally sits (sensor -> THIS -> DAC stage) and
// walks one output element from `low` to `high` in `step` increments, holding
// each level for `dwell` seconds. The trailing `settle_frac` of every level is
// averaged into one settled response, and a least-squares line through the
// (command, settled response) pairs gives the plant's DC gain -- the same
// quantity bode_plot reports as fit_K and latency_test as `gain`, but measured
// directly, at DC, and with the whole command range visible at once instead of
// a single small-signal amplitude.
//
// What the staircase adds over a two-level step is the shape: the per-level
// residuals show compression (the FSM's gain falling off at large command),
// `fit_range` re-fits the small-signal region the closed loop actually lives
// in, and `updown` measures hysteresis by walking back down again.
//
// Stages:
//   1. warmup   -- park at `bias` for `warmup` s so center_of_mass can acquire.
//   2. ramp in  -- glide from `bias` to `low` over `ramp` s. The staircase is
//                  entered gradually on purpose: a tracked center_of_mass
//                  window only follows motion that stays inside it, so a jump
//                  straight to `low` can lose the beam before the test starts.
//   3. sweep    -- one level per `dwell`, recording the settled mean, the
//                  within-level scatter and the sample count of each.
//   4. ramp out -- glide back to `bias` before anything else happens, so the
//                  actuator is parked whatever the analysis does.
//   5. report   -- fit, log, write the .dat and the PDF, raise AYLP_DONE.
//
// Two things make a level unusable, and both are recorded rather than fitted:
// a settled |response| past `resp_max` (the beam is running out of the sensor)
// and a settled window whose samples are all bit-identical (the centroid stage
// -- center_of_mass or fit_com -- has lost the beam and is holding its last
// coordinate, which otherwise reads as a beautifully quiet level). fit_com
// also flags this as AYLP_FRAME_REJECTED, but the bit-identical test catches
// both stages without either having to be assumed. Once the sweep has been in
// range, two unusable
// levels in a row end it early and the fit is done on what was measured.
//
// The bit-identical test only catches a level that is held THROUGHOUT: one
// differing sample clears it. A level that is partly held -- a chopped source,
// a beam that dims for a few frames -- is fitted with those held samples in the
// mean. That is harmless while the hold is short against the settled window,
// since the held value then comes from the same level, and it biases the slope
// toward zero once a hold spans a level boundary. Nothing here can tell the
// difference: watch the sensor's rejection rate.
//
// Params: see doc/devices/gain_test.md.

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

#include <gsl/gsl_fit.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_vector.h>

#include "anyloop.h"
#include "logging.h"
#include "gain_test.h"
#include "timing.h"
#include "xalloc.h"


static const char *PLOT_SCRIPT =
	"import sys\n"
	"import matplotlib\n"
	"matplotlib.use('PDF')\n"
	"import matplotlib.pyplot as plt\n"
	"import numpy as np\n"
	"d=np.loadtxt(sys.argv[1])\n"
	"if d.ndim==1: d=d[None,:]\n"
	"cmd,cv,r,sd,sem=d[:,1],d[:,2],d[:,3],d[:,4],d[:,5]\n"
	"resid,val,br=d[:,8],d[:,9].astype(int),d[:,10].astype(int)\n"
	"px=float(sys.argv[4])\n"
	"vpu=float(sys.argv[5])\n"
	"ok=val>0\n"
	"fig,(a1,a2)=plt.subplots(2,1,figsize=(9,8),sharex=True,\n"
	"    gridspec_kw={'height_ratios':[2,1]})\n"
	"fig.suptitle(sys.argv[3],fontsize=8)\n"
	"for b,mk,lb in ((0,'o','up'),(1,'s','down')):\n"
	"    m=ok&(br==b)\n"
	"    if m.any():\n"
	"        a1.errorbar(cv[m],r[m],yerr=sem[m],fmt=mk,ms=4,lw=1,\n"
	"            capsize=2,label=lb)\n"
	"        a2.errorbar(cv[m],resid[m]*px,yerr=sem[m]*px,fmt=mk,ms=4,lw=1,\n"
	"            capsize=2)\n"
	"if (~ok).any():\n"
	"    a1.plot(cv[~ok],r[~ok],'rx',ms=7,label='excluded')\n"
	"if len(sys.argv)>6 and sys.argv[6]:\n"
	"    s,i=[float(x) for x in sys.argv[6].split(',')]\n"
	"    xs=np.linspace(cmd.min(),cmd.max(),100)\n"
	"    a1.plot(xs*vpu,s*xs+i,'k--',lw=1,label='fit')\n"
	"a1.set_ylabel('settled response (normalized units)')\n"
	"a1.grid(True,ls='--',alpha=0.5)\n"
	"a1.legend(fontsize=8)\n"
	"b1=a1.twinx()\n"
	"b1.set_ylim([v*px for v in a1.get_ylim()])\n"
	"b1.set_ylabel('px')\n"
	"a2.axhline(0,color='k',lw=0.8)\n"
	"a2.set_ylabel('residual (px)')\n"
	"a2.set_xlabel('DAC output (V)')\n"
	"a2.grid(True,ls='--',alpha=0.5)\n"
	"plt.tight_layout()\n"
	"plt.savefig(sys.argv[2],format='pdf',bbox_inches='tight')\n"
	"print('Saved gain plot to '+sys.argv[2])\n";


static inline double gt_now(void)
{
	struct timespec t;
	clock_gettime(CLOCK_MONOTONIC, &t);
	return t.tv_sec + 1e-9 * t.tv_nsec;
}

/** Commanded value of level `i` of the staircase. Levels [0, n_levels) walk up
 * from `low`; with `updown`, levels [n_levels, 2*n_levels) walk back down, so
 * the top level is measured twice in a row -- a free repeatability check. */
static double gt_level_cmd(struct aylp_gain_test_data *data, size_t i)
{
	size_t n = data->n_levels;
	if (data->interleave) {
		size_t cycle = i / n;
		size_t j = i % n;
		size_t k = (j & 1) ? n - 1 - j/2 : j/2;
		if (cycle & 1) k = n - 1 - k;
		double v = data->low + (double)k * data->step;
		return v > data->high ? data->high : v;
	}
	if (i >= n) i = 2*n - 1 - i;
	double v = data->low + (double)i * data->step;
	return v > data->high ? data->high : v;
}

/** Least-squares line through the valid levels, plus the diagnostics that say
 * whether to believe it: R^2, residual rms, the largest single residual
 * (compression shows up here first), a slope from a fit that also carries a
 * linear-in-time term (beam drift over the sweep leaks straight into the plain
 * slope), an optional small-signal slope over |cmd - bias| <= fit_range, and
 * the mean up/down offset when the sweep ran both ways. */
static void gt_fit(struct aylp_gain_test_data *data)
{
	size_t n = 0;
	double *x = xmalloc(data->n_done * sizeof *x);
	double *y = xmalloc(data->n_done * sizeof *y);
	double *t = xmalloc(data->n_done * sizeof *t);
	for (size_t i = 0; i < data->n_done; i++) {
		if (!data->valid[i]) continue;
		x[n] = data->cmds[i];
		y[n] = data->resp[i];
		t[n] = data->times[i];
		n++;
	}
	// Repeated interleaved cycles exist specifically to average disturbance
	// that is independent of the commanded voltage. Fit one mean response per
	// command level rather than treating every visit as a separate point in
	// R2: the latter makes unavoidable beam motion look like a bad static
	// plant fit and does not improve when cycles are added. Raw visits remain
	// in resp[] and in the .dat output for repeatability/nonlinearity evidence.
	if (data->cycles > 1 && n > 0) {
		size_t grouped = 0;
		for (size_t i = 0; i < n; i++) {
			size_t g;
			for (g = 0; g < grouped; g++)
				if (fabs(x[g] - x[i]) <= 1e-9*fabs(data->step))
					break;
			if (g == grouped) {
				x[grouped] = x[i];
				y[grouped] = 0.0;
				t[grouped] = 0.0;
				grouped++;
			}
		}
		for (size_t g = 0; g < grouped; g++) {
			double sy = 0.0, st = 0.0;
			size_t count = 0;
			for (size_t i = 0; i < data->n_done; i++) {
				if (!data->valid[i]
					|| fabs(data->cmds[i] - x[g])
						> 1e-9*fabs(data->step))
					continue;
				sy += data->resp[i];
				st += data->times[i];
				count++;
			}
			y[g] = sy/count;
			t[g] = st/count;
		}
		log_info("gain_test: averaged %zu usable visits into %zu command "
			"levels before fitting", n, grouped);
		n = grouped;
	}
	if (n < 3) {
		log_error("gain_test: only %zu usable levels; not fitting. "
			"Did the beam leave the sensor, or is the sweep too "
			"wide?", n);
		goto out;
	}

	double c0, c1, cov00, cov01, cov11, sumsq;
	gsl_fit_linear(x, 1, y, 1, n, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
	data->intercept = c0;
	data->slope = c1;
	// The residual scatter, not the within-level noise, is what limits this
	// number: it carries the plant's own curvature and any drift, both of
	// which are real and neither of which shrinks by averaging longer.
	data->slope_err = sqrt(cov11);
	data->resid_rms = sqrt(sumsq / n);
	double ymean = 0.0, ymin = y[0], ymax = y[0];
	for (size_t i = 0; i < n; i++) {
		ymean += y[i];
		if (y[i] < ymin) ymin = y[i];
		if (y[i] > ymax) ymax = y[i];
	}
	ymean /= n;
	double tss = 0.0;
	for (size_t i = 0; i < n; i++) tss += (y[i]-ymean) * (y[i]-ymean);
	data->r2 = tss > 0.0 ? 1.0 - sumsq/tss : 0.0;
	data->span = ymax - ymin;
	data->nonlin = 0.0;
	for (size_t i = 0; i < data->n_done; i++) {
		if (!data->valid[i]) continue;
		double r = data->resp[i] - (c0 + c1 * data->cmds[i]);
		if (fabs(r) > data->nonlin) data->nonlin = fabs(r);
	}
	data->have_fit = true;

	// Drift-aware fit: y = a + b*cmd + c*t. Only possible when command and
	// time are not the same axis -- a plain rising staircase steps the
	// command in lockstep with the clock, so gain and drift are perfectly
	// collinear and no fit can tell one from the other. `updown` is what
	// breaks the tie: the down branch repeats every command at a later
	// time.
	double xm = 0.0, tm = 0.0;
	for (size_t i = 0; i < n; i++) { xm += x[i]; tm += t[i]; }
	xm /= n; tm /= n;
	double sxt = 0.0, sxx = 0.0, stt = 0.0;
	for (size_t i = 0; i < n; i++) {
		sxt += (x[i]-xm) * (t[i]-tm);
		sxx += (x[i]-xm) * (x[i]-xm);
		stt += (t[i]-tm) * (t[i]-tm);
	}
	double coll = (sxx > 0.0 && stt > 0.0)
		? fabs(sxt) / sqrt(sxx * stt) : 1.0;
	if (coll > 0.98) {
		log_info("gain_test: command and time are collinear (r = "
			"%.4f), so drift cannot be separated from gain in this "
			"sweep. Set updown to measure it: the down branch "
			"repeats each command later in the run.", coll);
	} else if (n >= 4) {
		gsl_matrix *X = gsl_matrix_alloc(n, 3);
		gsl_vector *yv = gsl_vector_alloc(n);
		gsl_vector *cv = gsl_vector_alloc(3);
		gsl_matrix *cov = gsl_matrix_alloc(3, 3);
		gsl_multifit_linear_workspace *w
			= gsl_multifit_linear_alloc(n, 3);
		for (size_t i = 0; i < n; i++) {
			gsl_matrix_set(X, i, 0, 1.0);
			gsl_matrix_set(X, i, 1, x[i]);
			gsl_matrix_set(X, i, 2, t[i]);
			gsl_vector_set(yv, i, y[i]);
		}
		double chisq;
		if (!gsl_multifit_linear(X, yv, cv, cov, &chisq, w)) {
			data->slope_drift = gsl_vector_get(cv, 1);
			data->drift_rate = gsl_vector_get(cv, 2);
		}
		gsl_multifit_linear_free(w);
		gsl_matrix_free(cov);
		gsl_vector_free(cv);
		gsl_vector_free(yv);
		gsl_matrix_free(X);
	}

	// small-signal fit
	if (data->fit_range > 0.0) {
		size_t m = 0;
		for (size_t i = 0; i < n; i++) {
			if (fabs(x[i] - data->bias) > data->fit_range) continue;
			x[m] = x[i];
			y[m] = y[i];
			m++;
		}
		if (m >= 3) {
			gsl_fit_linear(x, 1, y, 1, m, &c0, &c1,
				&cov00, &cov01, &cov11, &sumsq);
			data->slope_small = c1;
			data->slope_small_err = sqrt(cov11);
			data->n_small = m;
		} else {
			log_warn("gain_test: fit_range %G leaves only %zu "
				"levels; skipping the small-signal fit",
				data->fit_range, m);
		}
	}

	// hysteresis: mean (down - up) at matching commands
	if (data->updown) {
		double sum = 0.0;
		size_t np = 0;
		for (size_t i = 0; i < data->n_done; i++) {
			if (!data->valid[i] || data->branch[i] != GT_BRANCH_DOWN)
				continue;
			for (size_t j = 0; j < data->n_done; j++) {
				if (!data->valid[j]
					|| data->branch[j] != GT_BRANCH_UP
					|| fabs(data->cmds[j] - data->cmds[i])
						> 1e-9 * fabs(data->step))
					continue;
				sum += data->resp[i] - data->resp[j];
				np++;
				break;
			}
		}
		if (np) data->hysteresis = sum / np;
	}

out:
	xfree(x);
	xfree(y);
	xfree(t);
}

static int gt_write_dat(struct aylp_gain_test_data *data, const char *path,
	double fs)
{
	FILE *f = fopen(path, "w");
	if (!f) {
		log_error("gain_test: fopen %s: %m", path);
		return -1;
	}
	time_t tt = time(NULL);
	struct tm tm;
	localtime_r(&tt, &tm);
	double vpu = data->volts_per_unit;
	fprintf(f, "# gain_test %04d-%02d-%02d %02d:%02d:%02d\n",
		tm.tm_year+1900, tm.tm_mon+1, tm.tm_mday,
		tm.tm_hour, tm.tm_min, tm.tm_sec);
	fprintf(f, "# index_cmd %d index_err %d low %G high %G step %G "
		"bias %G dwell %G settle_frac %G updown %d\n",
		data->index_cmd, data->index_err, data->low, data->high,
		data->step, data->bias, data->dwell, data->settle_frac,
		data->updown);
	fprintf(f, "# volts_per_unit %G pixel_scale %G resp_max %G "
		"fs_measured %.2f\n", vpu, data->pixel_scale, data->resp_max,
		fs);
	if (data->have_fit) {
		fprintf(f, "# fit slope %.6G +/- %.4G intercept %.6G "
			"R2 %.6f resid_rms %.4G levels %zu/%zu\n",
			data->slope, data->slope_err, data->intercept,
			data->r2, data->resid_rms, data->n_done, data->n_points);
		fprintf(f, "# gain %.6G units/cmd", data->slope);
		if (vpu != 0.0)
			fprintf(f, " = %.6G units/V = %.4G px/V = %.4G mV/px",
				data->slope / vpu,
				data->slope * data->pixel_scale / vpu,
				1e3 * vpu
					/ (data->slope * data->pixel_scale));
		fputc('\n', f);
		fprintf(f, "# nonlinearity max_resid %.4G units = %.4G px "
			"= %.2f %% of span %.4G\n", data->nonlin,
			data->nonlin * data->pixel_scale,
			data->span > 0.0 ? 100.0*data->nonlin/data->span : 0.0,
			data->span);
		if (data->slope_drift != 0.0)
			fprintf(f, "# drift_fit slope %.6G drift %.4G units/s "
				"(%.4G px/s)\n", data->slope_drift,
				data->drift_rate,
				data->drift_rate * data->pixel_scale);
		if (data->n_small)
			fprintf(f, "# small_signal |cmd-bias|<=%G: slope %.6G "
				"+/- %.4G over %zu levels\n", data->fit_range,
				data->slope_small, data->slope_small_err,
				data->n_small);
		if (data->updown)
			fprintf(f, "# hysteresis mean(down-up) %.4G units "
				"= %.4G px\n", data->hysteresis,
				data->hysteresis * data->pixel_scale);
	} else {
		fprintf(f, "# NO FIT: too few usable levels\n");
	}
	if (data->config) fprintf(f, "# config %s\n", data->config);
	fprintf(f, "# level cmd_units cmd_V resp_mean resp_sd resp_sem "
		"n_samples t_s residual valid branch\n");
	for (size_t i = 0; i < data->n_done; i++) {
		double resid = data->have_fit
			? data->resp[i]
				- (data->intercept + data->slope*data->cmds[i])
			: 0.0;
		fprintf(f, "%zu %g %g %g %g %g %zu %g %g %d %d\n",
			i, data->cmds[i], data->cmds[i] * vpu, data->resp[i],
			data->resp_sd[i], data->resp_sem[i], data->counts[i],
			data->times[i], resid, data->valid[i], data->branch[i]);
	}
	fclose(f);
	log_info("gain_test: data saved to %s", path);
	return 0;
}

static int gt_plot(struct aylp_gain_test_data *data, const char *dat_path)
{
	char script_path[64];
	snprintf(script_path, sizeof script_path, "/tmp/gain_test_%d.py",
		getpid());
	FILE *script = fopen(script_path, "w");
	if (!script) {
		log_error("gain_test: fopen %s: %m", script_path);
		return -1;
	}
	fputs(PLOT_SCRIPT, script);
	fclose(script);

	char title[512];
	double vpu = data->volts_per_unit;
	if (data->have_fit && vpu != 0.0)
		snprintf(title, sizeof title,
			"K=%.4G units/V (%.4G px/V, %.4G mV/px)  R2=%.5f  "
			"resid rms %.3G px  max resid %.2f%% of span%s",
			data->slope / vpu,
			data->slope * data->pixel_scale / vpu,
			1e3 * vpu / (data->slope * data->pixel_scale),
			data->r2, data->resid_rms * data->pixel_scale,
			data->span > 0.0 ? 100.0*data->nonlin/data->span : 0.0,
			data->updown ? "  (up+down)" : "");
	else
		snprintf(title, sizeof title, "gain_test (no fit)");

	char fitarg[64] = "";
	if (data->have_fit)
		snprintf(fitarg, sizeof fitarg, "%.10G,%.10G",
			data->slope, data->intercept);

	char cmd[1024];
	snprintf(cmd, sizeof cmd, "python3 '%s' '%s' '%s' '%s' %G %G '%s'",
		script_path, dat_path, data->output_file, title,
		data->pixel_scale, data->volts_per_unit, fitarg);
	int ret = system(cmd);
	int ok = (ret != -1) && WIFEXITED(ret) && !WEXITSTATUS(ret);
	if (!ok)
		log_error("gain_test: plot script failed (exit status %d); is "
			"python3 with numpy and matplotlib installed?",
			(ret != -1 && WIFEXITED(ret)) ? WEXITSTATUS(ret) : ret);
	else
		log_info("gain_test: plot saved to %s", data->output_file);
	unlink(script_path);
	return ok ? 0 : -1;
}

/** Fit, log the headline numbers, and write the .dat and the PDF. */
/** Install the measured gain into a run config's fsp stage.
 *
 * Writes the SMALL-SIGNAL slope, not the full-span one: the closed loop lives
 * near the bias, and a mirror that compresses at large command would otherwise
 * hand the controller an averaged gain it never operates at. Magnitude only --
 * K is positive in these configs and the command->voltage sign lives in the DAC
 * stage's `scale`, where a sign error is a wiring question rather than a
 * silently inverted loop.
 *
 * The gates are the calibration suite's: a fit at all, a small-signal fit over
 * enough levels, R2 >= 0.90, uncertainty within 15%, and small-signal within
 * 25% of full-span. Failing any of them leaves the config untouched and says
 * why. A wrong K is a wrong loop gain into a fixed delay, so refusing to write
 * is always the cheaper mistake.
 */
static void gt_autowrite(struct aylp_gain_test_data *data)
{
	const char *path = data->write_config;
	const char *axis = data->write_axis ? "x" : "y";
	const char *why = 0;
	if (!data->have_fit) why = "there is no fit";
	else if (data->n_small < 5) why = "the small-signal fit has under 5 levels";
	else if (!(data->r2 >= 0.90)) why = "R2 is below 0.90";
	else if (!(fabs(data->slope_small) > 0.0)) why = "the slope is zero";
	else if (fabs(data->slope_small_err) > 0.15 * fabs(data->slope_small))
		why = "the small-signal uncertainty is over 15%";
	else if (fabs(fabs(data->slope_small) - fabs(data->slope))
			> 0.25 * fabs(data->slope_small))
		why = "small-signal and full-span gains differ by over 25%, "
			"which is too nonlinear to install";
	if (why) {
		log_error("gain_test: NOT writing %s to %s: %s. The run's .dat "
			"still has everything; install by hand if you disagree.",
			axis, path, why);
		return;
	}
	double K = fabs(data->slope_small);

	struct json_object *cfg = json_object_from_file(path);
	if (!cfg) {
		log_error("gain_test: could not read %s to write %s.K into",
			path, axis);
		return;
	}
	struct json_object *pipeline, *fsp = 0;
	if (json_object_object_get_ex(cfg, "pipeline", &pipeline)
			&& json_object_is_type(pipeline, json_type_array)) {
		size_t n = json_object_array_length(pipeline);
		for (size_t i = 0; i < n; i++) {
			struct json_object *dev, *uri;
			dev = json_object_array_get_idx(pipeline, i);
			if (!json_object_object_get_ex(dev, "uri", &uri))
				continue;
			if (strcmp(json_object_get_string(uri), "anyloop:fsp"))
				continue;
			if (fsp) {
				// two controllers means we cannot tell which one
				// this measurement belongs to
				log_error("gain_test: %s has more than one "
					"anyloop:fsp; not writing", path);
				json_object_put(cfg);
				return;
			}
			json_object_object_get_ex(dev, "params", &fsp);
		}
	}
	struct json_object *ax = 0;
	if (!fsp || !json_object_object_get_ex(fsp, axis, &ax)) {
		log_error("gain_test: %s has no anyloop:fsp stage with a \"%s\" "
			"axis; not writing", path, axis);
		json_object_put(cfg);
		return;
	}
	double was = 0.0;
	struct json_object *old;
	if (json_object_object_get_ex(ax, "K", &old))
		was = json_object_get_double(old);
	json_object_object_add(ax, "K", json_object_new_double(K));

	char note[512];
	snprintf(note, sizeof note,
		"K %.6G written by gain_test from the small-signal slope "
		"(%.6G +/- %.4G over %zu levels, full-span %.6G, R2 %.4f). "
		"Was %.6G. Re-run the gain test after any ROI, exposure or "
		"optics change: K is in units normalized over the frame.",
		K, data->slope_small, data->slope_small_err, data->n_small,
		data->slope, data->r2, was);
	json_object_object_add(ax, "_auto_gain_write",
		json_object_new_string(note));

	// write through a temporary so an interrupted write cannot leave a
	// half-written controller config behind
	char *tmp = xmalloc(strlen(path) + 5);
	sprintf(tmp, "%s.tmp", path);
	int err = json_object_to_file_ext(tmp, cfg,
		JSON_C_TO_STRING_PRETTY | JSON_C_TO_STRING_PRETTY_TAB);
	if (!err && rename(tmp, path))
		err = -1;
	if (err) {
		log_error("gain_test: could not write %s (%m)", path);
		remove(tmp);
	} else {
		log_info("gain_test: wrote %s.K = %.6G (was %.6G) into %s",
			axis, K, was, path);
	}
	xfree(tmp);
	json_object_put(cfg);
}


static void gt_report(struct aylp_gain_test_data *data)
{
	gt_fit(data);
	double fs = data->dt_cnt ? data->dt_cnt / data->dt_sum : 0.0;
	double vpu = data->volts_per_unit;

	log_info("gain_test: ================ RESULTS ================");
	log_info("gain_test: %zu of %zu levels usable; loop ran at %.1f Hz",
		data->n_done, data->n_points, fs);
	if (data->have_fit) {
		log_info("gain_test: slope %.6G +/- %.4G response units per "
			"command unit (intercept %.4G, R2 %.6f, residual rms "
			"%.4G units = %.4G px)", data->slope, data->slope_err,
			data->intercept, data->r2, data->resid_rms,
			data->resid_rms * data->pixel_scale);
		if (vpu != 0.0)
			log_info("gain_test: K = %.6G units/V = %.4G px/V = "
				"%.4G mV/px", data->slope / vpu,
				data->slope * data->pixel_scale / vpu,
				1e3 * vpu
					/ (data->slope * data->pixel_scale));
		log_info("gain_test: largest deviation from the line %.4G "
			"units = %.4G px = %.2f %% of the response span -- "
			"this is the plant's compression, not noise",
			data->nonlin, data->nonlin * data->pixel_scale,
			data->span > 0.0 ? 100.0*data->nonlin/data->span : 0.0);
		if (data->n_small)
			log_info("gain_test: small-signal slope over "
				"|cmd-bias| <= %G: %.6G +/- %.4G over %zu "
				"levels -- this is the number the closed loop "
				"lives at, if it differs from the full-range "
				"slope the plant is compressing",
				data->fit_range, data->slope_small,
				data->slope_small_err, data->n_small);
		if (data->slope_drift != 0.0)
			log_info("gain_test: with a linear drift term the "
				"slope is %.6G (drift %.4G px/s); a large gap "
				"to %.6G means the beam moved during the sweep",
				data->slope_drift,
				data->drift_rate * data->pixel_scale,
				data->slope);
		if (data->updown)
			log_info("gain_test: hysteresis (mean down - up) "
				"%.4G units = %.4G px", data->hysteresis,
				data->hysteresis * data->pixel_scale);
	}

	if (data->min_r2 > 0.0 &&
			(!data->have_fit || data->r2 < data->min_r2)) {
		log_error("gain_test: R2 %.6f is below min_r2 %.6f; refusing "
			"to publish a calibration report",
			data->have_fit ? data->r2 : 0.0, data->min_r2);
		return;
	}
	if (data->write_config) gt_autowrite(data);
	if (!data->output_file) return;
	const char *out = data->output_file;
	const char *dot = strrchr(out, '.');
	size_t base = dot ? (size_t)(dot - out) : strlen(out);
	char *dat_path = xmalloc(base + 5);
	memcpy(dat_path, out, base);
	memcpy(dat_path + base, ".dat", 5);
	if (!gt_write_dat(data, dat_path, fs)) gt_plot(data, dat_path);
	xfree(dat_path);
}


int gain_test_init(struct aylp_device *self)
{
	self->proc = &gain_test_proc;
	self->fini = &gain_test_fini;
	struct aylp_gain_test_data *data = xcalloc(1, sizeof *data);
	self->device_data = data;

	// defaults
	data->index_cmd = 1;
	data->index_err = 1;
	data->out_size = 2;
	data->low = -0.1;
	data->high = 0.1;
	data->step = 0.01;
	data->bias = 0.0;
	data->dwell = 1.0;
	data->settle_frac = 0.5;
	data->warmup = 5.0;
	data->ramp = 2.0;
	data->cycles = 1;
	data->volts_per_unit = 1.0;
	data->pixel_scale = 1.0;
	data->resp_max = 0.9;
	data->write_axis = -1;		// -1 = follow index_err
	data->output_file = xstrdup("gain_test.pdf");

	if (self->params) { json_object_object_foreach(self->params, key, val) {
		if (key[0] == '_') {
		} else if (!strcmp(key, "index_cmd")) {
			data->index_cmd = (int)json_object_get_int64(val);
		} else if (!strcmp(key, "index_err")) {
			data->index_err = (int)json_object_get_int64(val);
		} else if (!strcmp(key, "out_size")) {
			data->out_size = (size_t)json_object_get_uint64(val);
		} else if (!strcmp(key, "low")) {
			data->low = json_object_get_double(val);
		} else if (!strcmp(key, "high")) {
			data->high = json_object_get_double(val);
		} else if (!strcmp(key, "step")) {
			data->step = json_object_get_double(val);
		} else if (!strcmp(key, "bias")) {
			data->bias = json_object_get_double(val);
		} else if (!strcmp(key, "dwell")) {
			data->dwell = json_object_get_double(val);
		} else if (!strcmp(key, "settle_frac")) {
			data->settle_frac = json_object_get_double(val);
		} else if (!strcmp(key, "warmup")) {
			data->warmup = json_object_get_double(val);
		} else if (!strcmp(key, "ramp")) {
			data->ramp = json_object_get_double(val);
		} else if (!strcmp(key, "updown")) {
			data->updown = json_object_get_boolean(val);
		} else if (!strcmp(key, "interleave")) {
			data->interleave = json_object_get_boolean(val);
		} else if (!strcmp(key, "cycles")) {
			data->cycles = json_object_get_uint64(val);
		} else if (!strcmp(key, "use_rejected")) {
			data->use_rejected = json_object_get_boolean(val);
		} else if (!strcmp(key, "volts_per_unit")) {
			data->volts_per_unit = json_object_get_double(val);
		} else if (!strcmp(key, "pixel_scale")) {
			// "auto" = take it from the frame the source is actually
			// delivering, which is the only way to be right when the
			// source sizes its own ROI
			if (json_object_is_type(val, json_type_string)
					&& !strcmp(json_object_get_string(val),
						"auto"))
				data->pixel_scale_auto = true;
			else data->pixel_scale = json_object_get_double(val);
		} else if (!strcmp(key, "resp_max")) {
			data->resp_max = json_object_get_double(val);
		} else if (!strcmp(key, "fit_range")) {
			data->fit_range = json_object_get_double(val);
		} else if (!strcmp(key, "write_config")) {
			const char *v = json_object_get_string(val);
			xfree(data->write_config);
			if (v && *v) data->write_config = xstrdup(v);
		} else if (!strcmp(key, "write_axis")) {
			const char *v = json_object_get_string(val);
			if (v && (!strcmp(v, "x") || !strcmp(v, "y")))
				data->write_axis = !strcmp(v, "x");
			else {
				log_error("gain_test: write_axis must be "
					"\"x\" or \"y\"");
				return -1;
			}
		} else if (!strcmp(key, "min_r2")) {
			data->min_r2 = json_object_get_double(val);
		} else if (!strcmp(key, "output_file")) {
			// an empty name means "log the numbers, write nothing"
			const char *v = json_object_get_string(val);
			xfree(data->write_config);
	xfree(data->output_file);
			if (v && *v) data->output_file = xstrdup(v);
		} else if (!strcmp(key, "config")) {
			data->config = xstrdup(json_object_get_string(val));
		} else {
			log_warn("gain_test: unknown param \"%s\"", key);
		}
	} }

	if ((size_t)data->index_cmd >= data->out_size) {
		log_error("gain_test: index_cmd %d not < out_size %zu",
			data->index_cmd, data->out_size);
		return -1;
	}
	if (data->step <= 0.0 || data->high <= data->low) {
		log_error("gain_test: need low < high and step > 0 (got "
			"low %G high %G step %G)",
			data->low, data->high, data->step);
		return -1;
	}
	if (data->dwell <= 0.0
			|| data->settle_frac <= 0.0 || data->settle_frac >= 1.0) {
		log_error("gain_test: bad dwell/settle_frac");
		return -1;
	}
	double nspan = (data->high - data->low) / data->step;
	data->n_levels = (size_t)floor(nspan + 0.5) + 1;
	if (fabs(nspan - floor(nspan + 0.5)) > 1e-6)
		log_warn("gain_test: (high - low) is not a whole number of "
			"steps; the last level is clamped to high %G",
			data->high);
	if (data->n_levels < 3) {
		log_error("gain_test: %zu levels is not a line; reduce step",
			data->n_levels);
		return -1;
	}
	if (data->n_levels > 4096) {
		log_error("gain_test: %zu levels is absurd; raise step",
			data->n_levels);
		return -1;
	}
	if (!data->cycles) {
		log_error("gain_test: cycles must be nonzero");
		return -1;
	}
	if (!isfinite(data->min_r2) || data->min_r2 < 0.0 || data->min_r2 > 1.0) {
		log_error("gain_test: min_r2 must be in [0,1]");
		return -1;
	}
	data->n_points = data->interleave ? data->cycles * data->n_levels
		: (data->updown ? 2 * data->n_levels : data->n_levels);

	data->out = xcalloc_type(gsl_vector, data->out_size);
	data->cmds = xcalloc(data->n_points, sizeof *data->cmds);
	data->resp = xcalloc(data->n_points, sizeof *data->resp);
	data->resp_sd = xcalloc(data->n_points, sizeof *data->resp_sd);
	data->resp_sem = xcalloc(data->n_points, sizeof *data->resp_sem);
	data->times = xcalloc(data->n_points, sizeof *data->times);
	data->counts = xcalloc(data->n_points, sizeof *data->counts);
	data->valid = xcalloc(data->n_points, sizeof *data->valid);
	data->branch = xcalloc(data->n_points, sizeof *data->branch);
	data->stage = GT_WARMUP;
	data->cmd = data->bias;

	log_info("gain_test: %zu levels from %G to %G in steps of %G, %G s "
		"each%s; %.1f s of sweep after a %G s warmup and a %G s ramp",
		data->n_levels, data->low, data->high, data->step, data->dwell,
		data->interleave ? " (interleaved)" :
			(data->updown ? " (up then down)" : ""),
		data->n_points * data->dwell, data->warmup, data->ramp);
	if (data->write_axis < 0) data->write_axis = data->index_err ? 1 : 0;
	if (data->pixel_scale_auto) {
		// index_err 0 is y, 1 is x, matching center_of_mass's output
		struct aylp_frame_geometry g;
		if (aylp_frame_geometry_get(&g)) {
			size_t dim = data->index_err ? g.width : g.height;
			data->pixel_scale = 0.5 * (double)(dim - 1);
			log_info("gain_test: auto pixel_scale %G px/unit, from "
				"the source's %zux%zu frame on the %s axis",
				data->pixel_scale, g.height, g.width,
				data->index_err ? "x" : "y");
		} else {
			log_warn("gain_test: pixel_scale is \"auto\" but no "
				"source published a frame size; reporting in "
				"response units (px columns will read as "
				"units)");
		}
	}
	if (data->volts_per_unit != 1.0)
		log_info("gain_test: command unit is %G V at the DAC, so the "
			"sweep runs %G V to %G V in %G V steps",
			data->volts_per_unit,
			data->low * data->volts_per_unit,
			data->high * data->volts_per_unit,
			fabs(data->step * data->volts_per_unit));

	self->type_in = AYLP_T_VECTOR;
	self->units_in = AYLP_U_ANY;
	self->type_out = AYLP_T_VECTOR;
	self->units_out = AYLP_U_MINMAX;
	return 0;
}


int gain_test_proc(struct aylp_device *self, struct aylp_state *state)
{
	struct aylp_gain_test_data *data = self->device_data;
	double t = gt_now();
	int ret = 0;

	gsl_vector *s = state->vector;
	if (UNLIKELY((size_t)data->index_err >= s->size)) {
		log_error("gain_test: index_err %d out of range (input size "
			"%zu)", data->index_err, s->size);
		ret = -1;
		goto fatal;
	}
	double err = s->data[data->index_err * s->stride];

	if (UNLIKELY(!data->t0)) {
		data->t0 = t;
		data->stage_t0 = t;
		data->t_prev = t;
	}

	switch (data->stage) {
	case GT_WARMUP:
		data->cmd = data->bias;
		if (t - data->stage_t0 < data->warmup) break;
		data->stage = GT_RAMP_IN;
		data->stage_t0 = t;
		data->ramp_from = data->bias;
		break;

	case GT_RAMP_IN: {
		double target = gt_level_cmd(data, 0);
		double f = data->ramp > 0.0
			? (t - data->stage_t0) / data->ramp : 1.0;
		if (f >= 1.0) {
			data->cmd = target;
			data->stage = GT_SWEEP;
			data->stage_t0 = t;
			data->point = 0;
			data->wmean = data->wm2 = 0.0;
			data->wcnt = 0;
			log_info("gain_test: sweeping");
		} else {
			data->cmd = data->ramp_from
				+ f * (target - data->ramp_from);
		}
		break;
	}

	case GT_SWEEP: {
		data->dt_sum += t - data->t_prev;
		data->dt_cnt++;
		data->cmd = gt_level_cmd(data, data->point);
		// settled window: the tail of the level
		bool sensor_live = !(state->header.status
			& (AYLP_FRAME_REJECTED | AYLP_BEAM_LOST));
		if (t - data->stage_t0
				>= data->dwell * (1.0 - data->settle_frac)
				&& (data->use_rejected || sensor_live)) {
			if (!data->wcnt) {
				data->wfirst = err;
				data->wmoved = false;
			} else if (err != data->wfirst) {
				data->wmoved = true;
			}
			data->wcnt++;
			double d = err - data->wmean;
			data->wmean += d / data->wcnt;
			data->wm2 += d * (err - data->wmean);
		}
		if (t - data->stage_t0 < data->dwell) break;

		// level over: fold it into the results
		size_t i = data->n_done;
		if (data->wcnt < 8) {
			log_error("gain_test: only %zu settled samples at "
				"level %zu; dwell %G s is too short for this "
				"loop rate", data->wcnt, data->point,
				data->dwell);
			ret = -1;
			goto fatal;
		}
		double mean = data->wmean;
		double var = data->wcnt > 1 ? data->wm2 / (data->wcnt - 1) : 0.0;
		data->cmds[i] = data->cmd;
		data->resp[i] = mean;
		data->resp_sd[i] = sqrt(var);
		data->resp_sem[i] = sqrt(var / data->wcnt);
		data->counts[i] = data->wcnt;
		data->times[i] = t - data->dwell*data->settle_frac/2.0 - data->t0;
		data->branch[i] = data->interleave
			? ((data->point / data->n_levels) & 1)
			: (data->point < data->n_levels
				? GT_BRANCH_UP : GT_BRANCH_DOWN);
		// A level is unusable if the beam is running out of the sensor,
		// or if the sensor stopped moving altogether -- a tracked
		// centroid stage that has lost the beam (center_of_mass past
		// min_peak, fit_com on a failed fit) holds its last coordinate,
		// which reads as a perfectly quiet level rather than as an
		// error.
		bool edge = fabs(mean) > data->resp_max;
		// exact equality, not a small variance: a held coordinate is
		// bit-identical every frame, while a real sensor that happens
		// to be very quiet is not
		bool frozen = !data->wmoved;
		data->valid[i] = !edge && !frozen;
		if (edge)
			log_warn("gain_test: level %zu (cmd %G): |response| "
				"%.4G is past resp_max %G -- the beam is "
				"leaving the sensor; excluded from the fit",
				data->point, data->cmd, fabs(mean),
				data->resp_max);
		else if (frozen)
			log_warn("gain_test: level %zu (cmd %G): response is "
				"perfectly constant -- the sensor has lost the "
				"beam and is holding its last value; excluded",
				data->point, data->cmd);
		else
			log_info("gain_test: level %zu/%zu: cmd %G -> %.5G "
				"+/- %.2G (%zu samples)", data->point + 1,
				data->n_points, data->cmd, mean,
				data->resp_sem[i], data->wcnt);
		data->n_bad_run = data->valid[i] ? 0 : data->n_bad_run + 1;
		if (data->valid[i]) data->seen_valid = true;
		data->n_done++;

		data->wmean = data->wm2 = 0.0;
		data->wcnt = 0;
		data->stage_t0 = t;
		data->point++;
		// Two unusable levels in a row end the sweep -- but only once
		// the sweep has been in range at least once. A staircase that
		// starts wider than the sensor walks INTO range, and aborting
		// on its first levels would throw away the good middle; one
		// that leaves range at the far end is walking the beam further
		// off with every remaining level, and should stop.
		if (data->n_bad_run >= 2 && data->seen_valid) {
			log_warn("gain_test: two unusable levels in a row after "
				"the sweep had been in range; ending early and "
				"fitting the %zu levels measured so far",
				data->n_done);
			data->point = data->n_points;
		}
		if (data->point >= data->n_points) {
			data->stage = GT_RAMP_OUT;
			data->ramp_from = data->cmd;
		} else {
			// publish the new level in the same iteration its dwell
			// starts, rather than one frame into it
			data->cmd = gt_level_cmd(data, data->point);
		}
		break;
	}

	case GT_RAMP_OUT: {
		double f = data->ramp > 0.0
			? (t - data->stage_t0) / data->ramp : 1.0;
		if (f >= 1.0) {
			data->cmd = data->bias;
			// park first, analyse second: the report writes files
			// and runs a plot script, and the actuator should
			// already be sitting at bias while that happens
			data->stage = GT_DONE;
			gt_report(data);
			state->header.status |= AYLP_DONE;
		} else {
			data->cmd = data->ramp_from
				+ f * (data->bias - data->ramp_from);
		}
		break;
	}

	case GT_DONE:
		data->cmd = data->bias;
		break;
	}
	data->t_prev = t;
	goto publish;

fatal:
	// proc errors are recoverable as far as the main loop is concerned, so
	// a fatal test condition has to end the loop itself -- and still
	// publish a parked command, or the DAC stage would be handed the
	// sensor's vector and left holding the last level
	data->stage = GT_DONE;
	data->cmd = data->bias;
	if (data->n_done >= 3) gt_report(data);
	state->header.status |= AYLP_DONE;

publish:
	gsl_vector_set_zero(data->out);
	data->out->data[data->index_cmd * data->out->stride] = data->cmd;
	state->vector = data->out;
	state->header.type = self->type_out;
	state->header.units = self->units_out;
	state->header.log_dim.y = data->out->size;
	state->header.log_dim.x = 1;
	return ret;
}


int gain_test_fini(struct aylp_device *self)
{
	struct aylp_gain_test_data *data = self->device_data;
	xfree(data->output_file);
	xfree(data->config);
	xfree_type(gsl_vector, data->out);
	xfree(data->cmds);
	xfree(data->resp);
	xfree(data->resp_sd);
	xfree(data->resp_sem);
	xfree(data->times);
	xfree(data->counts);
	xfree(data->valid);
	xfree(data->branch);
	xfree(data);
	return 0;
}
