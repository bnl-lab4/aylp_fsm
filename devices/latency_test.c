// End-to-end actuation latency test.
//
// Sits where the pid normally sits (sensor -> THIS -> DAC stage) and toggles
// one output element between `low` and `high` every `period` seconds, while
// watching one input element (the CoM error) for the response. The command
// edge is timestamped in the same iteration that hands the new value to the
// DAC stage, and the response is timestamped at the proc that first carries
// the moved centroid, so the measured delay covers the whole real path:
// serial write -> DAC -> actuator motion -> exposure -> camera readout/USB ->
// centroid computation. Exactly the loop's feedback-path delay.
//
// Stages:
//   1. warmup  -- hold `low` for `warmup` s so center_of_mass can acquire.
//   2. cal     -- run `cal_cycles` low/high cycles; the last `settle_frac` of
//                 each half-cycle is averaged into a settled window mean.
//                 The step is the mean difference of consecutive high/low
//                 window means (differencing cancels slow open-loop beam
//                 drift) and the noise floor is the within-window scatter
//                 only, so drift doesn't count as noise. The test refuses to
//                 continue if the step is < 8 sigma of that noise.
//   3. measure -- for `n_steps` edges, record two latencies per edge:
//                 departure: first of 2 consecutive samples past
//                   max(4 sigma, 2% of step) away from the baseline toward
//                   the target (transport delay + first motion);
//                 50% crossing: first sample past baseline + step/2.
//                 The baseline is the settled mean of the half-cycle the
//                 edge departs from, re-measured every edge, so slow beam
//                 drift between edges cannot stale it.
//   4. report  -- log median/mean/min/max of both, the mean sample interval,
//                 and the latency in frames; then park the command at 0 and
//                 set AYLP_DONE.
//
// Every measured edge is also logged as it happens, so a wedged test is
// visible immediately.
//
// Params:
//   index_cmd (int)    -- output element to drive (default 0 = x_fine)
//   index_err (int)    -- input element to watch (default 1 = x, CoM is [y,x])
//   out_size (int)     -- output vector length (default 4)
//   low, high (float)  -- command levels in minmax units (default -+0.05).
//                         Keep the resulting beam step well inside the
//                         center_of_mass window or tracking will be lost.
//   period (float)     -- seconds per half-cycle (default 0.25)
//   warmup (float)     -- seconds to hold `low` first (default 6; must exceed
//                         center_of_mass's acquire_seconds)
//   cal_cycles (int)   -- calibration cycles (default 3)
//   n_steps (int)      -- edges to measure (default 40)
//   settle_frac (float)-- settled fraction of each half-cycle (default 0.4)
//   results_file (str) -- append the full results (run header, calibration,
//                         every edge, summary stats) to this file. Appended,
//                         not truncated, so one file accumulates a history of
//                         runs; each run starts with a wall-clock timestamp.
//   label (str)        -- free-text label written into the run header in
//                         results_file, to tell runs apart.

#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "anyloop.h"
#include "logging.h"
#include "xalloc.h"
#include "latency_test.h"


static double lt_monotonic_s(void)
{
	struct timespec t;
	clock_gettime(CLOCK_MONOTONIC, &t);
	return t.tv_sec + 1e-9 * t.tv_nsec;
}

static int lt_cmp_dbl(const void *a, const void *b)
{
	double d = *(const double *)a - *(const double *)b;
	return (d > 0) - (d < 0);
}

/** Append a line to the results file, if one is open. Flushed immediately so
 * a wedged or killed test still leaves everything measured so far on disk. */
static void lt_fout(struct aylp_latency_test_data *data, const char *fmt, ...)
{
	if (!data->rf) return;
	va_list ap;
	va_start(ap, fmt);
	vfprintf(data->rf, fmt, ap);
	va_end(ap);
	fflush(data->rf);
}

/** Log median/mean/min/max of n latencies (s). Sorts arr in place. */
static void lt_report_one(struct aylp_latency_test_data *data,
	const char *name, const char *key, double *arr, size_t n)
{
	if (!n) {
		log_warn("latency_test: %s: no edges detected!", name);
		lt_fout(data, "# %s: no edges detected\n", key);
		return;
	}
	qsort(arr, n, sizeof *arr, lt_cmp_dbl);
	double sum = 0.0;
	for (size_t i = 0; i < n; i++) sum += arr[i];
	log_info("latency_test: %s: median %.2f ms | mean %.2f ms | "
		"min %.2f ms | max %.2f ms | n=%zu", name,
		1e3 * arr[n/2], 1e3 * sum / n,
		1e3 * arr[0], 1e3 * arr[n-1], n);
	lt_fout(data, "# %s median %.3f mean %.3f min %.3f max %.3f n %zu\n",
		key, 1e3 * arr[n/2], 1e3 * sum / n,
		1e3 * arr[0], 1e3 * arr[n-1], n);
}


int latency_test_init(struct aylp_device *self)
{
	self->proc = &latency_test_proc;
	self->fini = &latency_test_fini;
	struct aylp_latency_test_data *data = xcalloc(1, sizeof *data);
	self->device_data = data;

	// defaults
	data->index_cmd = 0;
	data->index_err = 1;
	data->out_size = 4;
	data->low = -0.05;
	data->high = 0.05;
	data->period = 0.25;
	data->warmup = 6.0;
	data->settle_frac = 0.4;
	data->cal_cycles = 3;
	data->n_steps = 40;

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
		} else if (!strcmp(key, "period")) {
			data->period = json_object_get_double(val);
		} else if (!strcmp(key, "warmup")) {
			data->warmup = json_object_get_double(val);
		} else if (!strcmp(key, "cal_cycles")) {
			data->cal_cycles = (size_t)json_object_get_uint64(val);
		} else if (!strcmp(key, "n_steps")) {
			data->n_steps = (size_t)json_object_get_uint64(val);
		} else if (!strcmp(key, "settle_frac")) {
			data->settle_frac = json_object_get_double(val);
		} else if (!strcmp(key, "results_file")) {
			data->results_file = strdup(json_object_get_string(val));
		} else if (!strcmp(key, "label")) {
			data->label = strdup(json_object_get_string(val));
		} else {
			log_warn("latency_test: unknown param \"%s\"", key);
		}
	} }

	if ((size_t)data->index_cmd >= data->out_size) {
		log_error("latency_test: index_cmd %d not < out_size %zu",
			data->index_cmd, data->out_size);
		return -1;
	}
	if (data->period <= 0.0 || data->settle_frac <= 0.0
			|| data->settle_frac >= 1.0
			|| !data->cal_cycles || !data->n_steps) {
		log_error("latency_test: bad period/settle_frac/cal_cycles/"
			"n_steps");
		return -1;
	}
	if (data->low == data->high) {
		log_error("latency_test: low == high; no step to measure");
		return -1;
	}

	data->out = xcalloc_type(gsl_vector, data->out_size);
	data->dep = xcalloc(data->n_steps, sizeof *data->dep);
	data->half = xcalloc(data->n_steps, sizeof *data->half);
	data->stage = LT_WARMUP;
	data->dep_pend = -1.0;

	if (data->results_file) {
		data->rf = fopen(data->results_file, "a");
		if (!data->rf) {
			log_error("latency_test: cannot open results_file "
				"\"%s\"", data->results_file);
			return -1;
		}
		time_t now = time(0);
		char ts[32];
		strftime(ts, sizeof ts, "%Y-%m-%dT%H:%M:%S", localtime(&now));
		lt_fout(data, "# ===== latency_test run %s%s%s =====\n", ts,
			data->label ? " label " : "",
			data->label ? data->label : "");
		lt_fout(data, "# index_cmd %d index_err %d low %G high %G "
			"period %G warmup %G cal_cycles %zu n_steps %zu "
			"settle_frac %G\n",
			data->index_cmd, data->index_err, data->low,
			data->high, data->period, data->warmup,
			data->cal_cycles, data->n_steps, data->settle_frac);
		lt_fout(data, "# columns: edge direction departure_ms half_ms "
			"(-1 = not detected)\n");
	}

	log_info("latency_test: stepping output %d between %G and %G every "
		"%G s; watching input %d; %zu cal cycles then %zu edges",
		data->index_cmd, data->low, data->high, data->period,
		data->index_err, data->cal_cycles, data->n_steps);

	self->type_in = AYLP_T_VECTOR;
	self->units_in = AYLP_U_ANY;
	self->type_out = AYLP_T_VECTOR;
	self->units_out = AYLP_U_MINMAX;
	return 0;
}


int latency_test_proc(struct aylp_device *self, struct aylp_state *state)
{
	struct aylp_latency_test_data *data = self->device_data;
	double t = lt_monotonic_s();
	int ret = 0;

	gsl_vector *s = state->vector;
	if (UNLIKELY((size_t)data->index_err >= s->size)) {
		log_error("latency_test: index_err %d out of range (input "
			"size %zu)", data->index_err, s->size);
		ret = -1;
		goto fatal;
	}
	double err = s->data[data->index_err * s->stride];

	if (UNLIKELY(!data->t0)) data->t0 = t;

	switch (data->stage) {
	case LT_WARMUP:
		// hold low; the sensor side (acquisition, tracking) settles
		if (t - data->t0 >= data->warmup) {
			data->stage = LT_CAL;
			data->level_high = true;	// first toggle: low -> high
			data->half_t0 = t;
			data->halfcycles = 0;
			log_info("latency_test: warmup done; calibrating");
		}
		break;

	case LT_CAL: {
		// settled window: the tail of the half-cycle
		if (t - data->half_t0 >= data->period * (1.0 - data->settle_frac)) {
			data->wsum += err;
			data->wsum2 += err * err;
			data->wcnt++;
		}
		if (t - data->half_t0 < data->period) break;
		// half-cycle over: fold its settled window into the calibration
		if (data->wcnt < 8) {
			log_error("latency_test: only %zu settled samples in "
				"a half-cycle; period too short for this loop "
				"rate?", data->wcnt);
			ret = -1;
			goto fatal;
		}
		int lv = data->level_high;
		double wmean = data->wsum / data->wcnt;
		double wvar = data->wsum2 / data->wcnt - wmean * wmean;
		if (wvar < 0.0) wvar = 0.0;
		data->var_sum += wvar * data->wcnt;
		data->var_cnt += data->wcnt;
		data->msum[lv] += wmean;
		data->mcnt[lv]++;
		if (data->win_prev_valid) {
			// consecutive windows alternate level; sign as high - low
			data->step_sum += lv ? wmean - data->win_mean_prev
				: data->win_mean_prev - wmean;
			data->step_pairs++;
		}
		data->win_mean_prev = wmean;
		data->win_prev_valid = true;
		data->wsum = data->wsum2 = 0.0;
		data->wcnt = 0;
		// toggle
		data->level_high = !data->level_high;
		data->half_t0 = t;
		data->halfcycles++;
		if (data->halfcycles < 2 * data->cal_cycles) break;
		// calibration complete. Differencing consecutive windows for the
		// step cancels slow beam drift; taking the noise from within-
		// window scatter only keeps that same drift out of sigma.
		data->sigma = sqrt(data->var_sum / data->var_cnt);
		data->step = data->step_sum / data->step_pairs;
		// plant DC gain: settled error response per unit minmax command.
		// The step is the response to a (high - low) command change, so
		// dividing by that gives error-units/command-unit -- the same
		// quantity (and units) bode_plot reports as fit_K, from the same
		// point on the plant, so a latency run doubles as a DC-gain check.
		data->gain = data->step / (data->high - data->low);
		log_info("latency_test: calibrated: settled error %G (low) / "
			"%G (high), step %G, fast noise sigma %G, step/sigma "
			"%.1f, gain (px/cmd) %G", data->msum[0] / data->mcnt[0],
			data->msum[1] / data->mcnt[1], data->step, data->sigma,
			data->sigma > 0
				? fabs(data->step) / data->sigma : INFINITY,
			data->gain);
		lt_fout(data, "# calibrated: settled_low %G settled_high %G "
			"step %G sigma %G step/sigma %.1f gain %G\n",
			data->msum[0] / data->mcnt[0],
			data->msum[1] / data->mcnt[1], data->step, data->sigma,
			data->sigma > 0
				? fabs(data->step) / data->sigma : INFINITY,
			data->gain);
		if (fabs(data->step) <= 0.0
				|| fabs(data->step) < 8.0 * data->sigma) {
			log_error("latency_test: step is < 8 sigma of noise; "
				"raise |high - low| or fix the beam/tracking. "
				"Refusing to measure garbage.");
			ret = -1;
			goto fatal;
		}
		// this toggle is the first measured edge; it departs from the
		// settled level of the window that just closed
		data->from_level = wmean;
		data->stage = LT_MEAS;
		data->dep_found = data->half_found = false;
		data->dep_pend = -1.0;
		data->t_prev = t;
		log_info("latency_test: measuring %zu edges...", data->n_steps);
		break;
	}

	case LT_MEAS: {
		data->dt_sum += t - data->t_prev;
		data->dt_cnt++;
		data->t_prev = t;

		// signed expected response for this edge: the command went
		// toward `to`, so the error should move by delta from the
		// baseline (the settled level of the previous half-cycle)
		int to = data->level_high;
		double delta = to ? data->step : -data->step;
		double dir = delta > 0.0 ? 1.0 : -1.0;
		double step = fabs(delta);

		if (!data->dep_found) {
			// first movement: past max(4 sigma, 2% of step), on 2
			// consecutive samples so a lone noise spike can't fire
			double thr = 4.0 * data->sigma;
			if (0.02 * step > thr) thr = 0.02 * step;
			if ((err - data->from_level) * dir > thr) {
				if (data->dep_pend >= 0.0) {
					data->dep_found = true;
					data->dep[data->n_dep++] =
						data->dep_pend - data->half_t0;
				} else {
					data->dep_pend = t;
				}
			} else {
				data->dep_pend = -1.0;
			}
		}
		if (!data->half_found
				&& (err - data->from_level) * dir > step / 2.0) {
			data->half_found = true;
			data->half[data->n_half++] = t - data->half_t0;
		}

		// settled tail of this half-cycle: the next edge's baseline
		if (t - data->half_t0 >= data->period * (1.0 - data->settle_frac)) {
			data->wsum += err;
			data->wsum2 += err * err;
			data->wcnt++;
		}

		if (t - data->half_t0 < data->period) break;
		// half-cycle over: close out this edge
		data->edges_done++;
		if (data->dep_found && data->half_found) {
			log_info("latency_test: edge %zu/%zu (%s): departure "
				"%.2f ms, 50%% crossing %.2f ms",
				data->edges_done, data->n_steps,
				to ? "low->high" : "high->low",
				1e3 * data->dep[data->n_dep-1],
				1e3 * data->half[data->n_half-1]);
		} else {
			data->edges_missed++;
			log_warn("latency_test: edge %zu/%zu: response not "
				"detected (departure %s, 50%% %s)",
				data->edges_done, data->n_steps,
				data->dep_found ? "yes" : "NO",
				data->half_found ? "yes" : "NO");
		}
		lt_fout(data, "%zu %s %.3f %.3f\n",
			data->edges_done, to ? "low->high" : "high->low",
			data->dep_found ? 1e3 * data->dep[data->n_dep-1] : -1.0,
			data->half_found
				? 1e3 * data->half[data->n_half-1] : -1.0);
		if (data->edges_done >= data->n_steps) {
			// report
			double mean_dt = data->dt_cnt
				? data->dt_sum / data->dt_cnt : 0.0;
			log_info("latency_test: ================ RESULTS "
				"================");
			log_info("latency_test: %zu edges measured, %zu "
				"missed", data->edges_done, data->edges_missed);
			log_info("latency_test: plant DC gain = %G px per "
				"unit minmax command (step %G / cmd %G); "
				"compare to bode fit_K", data->gain,
				data->step, data->high - data->low);
			lt_fout(data, "# edges_measured %zu edges_missed %zu\n",
				data->edges_done, data->edges_missed);
			lt_fout(data, "# gain %G step %G cmd %G\n",
				data->gain, data->step,
				data->high - data->low);
			lt_report_one(data, "departure (first motion)",
				"departure_ms", data->dep, data->n_dep);
			lt_report_one(data, "50% crossing",
				"half_ms", data->half, data->n_half);
			if (mean_dt > 0.0) {
				log_info("latency_test: sample interval %.3f "
					"ms (%.0f Hz loop)",
					1e3 * mean_dt, 1.0 / mean_dt);
				lt_fout(data, "# sample_interval_ms %.3f "
					"loop_hz %.0f\n",
					1e3 * mean_dt, 1.0 / mean_dt);
				if (data->n_dep) {
					log_info("latency_test: median "
						"departure = %.1f frames",
						data->dep[data->n_dep/2]
						/ mean_dt);
					lt_fout(data, "# median_departure_"
						"frames %.1f\n",
						data->dep[data->n_dep/2]
						/ mean_dt);
				}
			}
			// park the command at bias and stop the loop
			gsl_vector_set_zero(data->out);
			data->stage = LT_DONE;
			state->header.status |= AYLP_DONE;
		} else {
			// re-baseline on the settled level actually reached; if
			// the loop hiccupped and the window is empty, predict it
			// from the old baseline instead
			if (data->wcnt >= 4)
				data->from_level = data->wsum / data->wcnt;
			else
				data->from_level += delta;
			data->wsum = data->wsum2 = 0.0;
			data->wcnt = 0;
			// next edge
			data->level_high = !data->level_high;
			data->half_t0 = t;
			data->dep_found = data->half_found = false;
			data->dep_pend = -1.0;
		}
		break;
	}

	case LT_DONE:
		break;
	}
	goto publish;

fatal:
	// proc errors are treated as recoverable by the main loop, so a fatal
	// test condition must end the loop itself -- and still publish a
	// parked output vector, or the DAC stage would be handed the sensor's
	// vector and be left holding the last stepped command
	lt_fout(data, "# FATAL: test aborted before completion (see log)\n");
	data->stage = LT_DONE;
	gsl_vector_set_zero(data->out);
	state->header.status |= AYLP_DONE;

publish:
	if (data->stage != LT_DONE)
		data->out->data[data->index_cmd * data->out->stride] =
			data->level_high ? data->high : data->low;

	state->vector = data->out;
	state->header.type = self->type_out;
	state->header.units = self->units_out;
	state->header.log_dim.y = data->out->size;
	state->header.log_dim.x = 1;
	return ret;
}


int latency_test_fini(struct aylp_device *self)
{
	struct aylp_latency_test_data *data = self->device_data;
	if (data->rf) fclose(data->rf);
	free(data->results_file);
	free(data->label);
	xfree_type(gsl_vector, data->out);
	xfree(data->dep);
	xfree(data->half);
	xfree(data);
	return 0;
}
