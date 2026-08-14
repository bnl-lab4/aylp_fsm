// End-to-end loop delay from a burst of pseudorandom binary drive.
//
// Sits where the controller normally sits (sensor -> THIS -> DAC stage) and
// drives one output element with a maximum-length sequence: +-`amplitude`
// about `bias`, one chip every `chip_frames` loop iterations, `2^order - 1`
// chips per burst, `n_bursts` bursts separated by `quiet_frames` of rest. The
// response of one input element is recorded frame by frame, and the delay
// comes out of the cross-correlation between the command that was emitted and
// the response that came back.
//
// Why a burst rather than one long sequence: every burst is driven from the
// same LFSR seed, so all of them are the *same* stimulus and their responses
// can be ensemble-averaged before anything else happens. Room tone, the 5-80 Hz
// disturbance lines and centroid noise are uncorrelated with the sequence, so
// they fall as 1/sqrt(n_bursts) while the response does not. The quiet tail
// also lets the actuator return to rest, so each burst starts from the same
// state, and it gives the correlogram somewhere to show its own noise floor.
//
// Three delays are reported, because they are three different things and the
// gap between them is the plant's rise, not disagreement:
//   - onset: the first lag where the correlation leaves the noise. This is the
//     transport delay -- DAC write, actuator dead time, exposure, readout,
//     centroid.
//   - peak: the lag of the largest correlation, refined to sub-frame by
//     parabolic interpolation. Transport delay plus the plant's rise.
//   - phase slope: the group delay from a weighted straight-line fit to the
//     phase of the cross-spectrum over [phase_f_lo, phase_f_hi]. This is the
//     same estimator (and the same band) that bode_plot's fit_tau_ms uses, so
//     it is the number to compare against a sweep, and the one to put in an
//     fsp `delay`/`delay_frac` pair.
// Per-burst estimates are made as well; their median absolute deviation is the
// honest uncertainty on the ensemble numbers.
//
// Only frames the sensor says it measured are correlated. A centroid stage
// that cannot fit a frame -- a chopped or blinking source, a beam that dims --
// re-publishes its previous coordinate and raises AYLP_FRAME_REJECTED, so the
// stream contains samples carrying an older frame's value. Those are not
// missing data but wrong data, and correlating them pulls the apparent lag
// toward the age of what was held: at a 25% hold rate with holds averaging 7
// frames old, the lag reads roughly 1.5 frames long. Every frame's flag is
// recorded and the held ones are dropped from the ensemble mean and from every
// correlation, so no duty cycle is assumed and none has to be periodic. What
// makes this cheap is that the chop is not locked to the burst period: a
// window position dark in one burst is live in others, so the ensemble mean
// fills in. `use_rejected` restores the old behaviour for comparison.
//
// A lag of L frames means the response appeared L iterations after the command
// was emitted. L = 0 is physically impossible -- the sensor value read in an
// iteration was captured before that iteration's command went out -- so a
// measured lag below 1 frame means the correlation found noise, not the plant.
// There is no hidden frame to add back the way there is in bode_plot: the
// command and the response are correlated exactly as the loop sees them.
//
// Params: see doc/devices/prbs_test.md.

#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

#include <gsl/gsl_fft_real.h>

#include "anyloop.h"
#include "logging.h"
#include "prbs_test.h"
#include "xalloc.h"


static const char *PLOT_SCRIPT =
	"import sys\n"
	"import matplotlib\n"
	"matplotlib.use('PDF')\n"
	"import matplotlib.pyplot as plt\n"
	"import numpy as np\n"
	"c=np.loadtxt(sys.argv[1])\n"
	"tr=np.loadtxt(sys.argv[2])\n"
	"onset,peak,px=float(sys.argv[5]),float(sys.argv[6]),float(sys.argv[7])\n"
	"fig,(a1,a2)=plt.subplots(2,1,figsize=(10,8))\n"
	"fig.suptitle(sys.argv[4],fontsize=8)\n"
	"n=min(len(tr),int(sys.argv[8]))\n"
	"a1.step(tr[:n,1],tr[:n,2],where='post',lw=1,color='tab:gray',\n"
	"    label='command')\n"
	"a1.set_xlabel('time into burst (ms)')\n"
	"a1.set_ylabel('command (units)')\n"
	"b1=a1.twinx()\n"
	"b1.errorbar(tr[:n,1],tr[:n,3]*px,yerr=tr[:n,4]*px,fmt='-',lw=1,\n"
	"    color='tab:blue',elinewidth=0.5,label='response')\n"
	"b1.set_ylabel('ensemble-mean response (px)',color='tab:blue')\n"
	"a1.grid(True,ls='--',alpha=0.5)\n"
	"a1.set_title('stimulus and averaged response',fontsize=9)\n"
	"a2.axhline(0,color='k',lw=0.8)\n"
	"a2.plot(c[:,0],c[:,2],'.-',lw=1,ms=3)\n"
	"a2.axvline(onset,color='tab:green',ls='--',lw=1,\n"
	"    label='onset %.2f fr'%onset)\n"
	"a2.axvline(peak,color='tab:red',ls='--',lw=1,label='peak %.2f fr'%peak)\n"
	"a2.axvspan(c[0,0],0,color='0.9',label='acausal (noise floor)')\n"
	"a2.set_xlabel('lag (frames)')\n"
	"a2.set_ylabel('normalized cross-correlation')\n"
	"a2.grid(True,ls='--',alpha=0.5)\n"
	"a2.legend(fontsize=8)\n"
	"a2.set_title('command-response cross-correlation',fontsize=9)\n"
	"plt.tight_layout()\n"
	"plt.savefig(sys.argv[3],format='pdf',bbox_inches='tight')\n"
	"print('Saved PRBS plot to '+sys.argv[3])\n";


static inline double pt_now(void)
{
	struct timespec t;
	clock_gettime(CLOCK_MONOTONIC, &t);
	return t.tv_sec + 1e-9 * t.tv_nsec;
}

/** Galois LFSR feedback masks giving a maximal-length sequence, i.e. one that
 * visits all 2^order - 1 nonzero states. Verified at init by walking the whole
 * cycle, so a wrong entry here fails loudly instead of quietly producing a
 * short sequence with a periodic autocorrelation. */
static unsigned pt_lfsr_mask(unsigned order)
{
	switch (order) {
	case  5: return 0x0014u;
	case  6: return 0x0030u;
	case  7: return 0x0060u;
	case  8: return 0x00B8u;
	case  9: return 0x0110u;
	case 10: return 0x0240u;
	case 11: return 0x0500u;
	case 12: return 0x0829u;
	case 13: return 0x100Du;
	case 14: return 0x2015u;
	case 15: return 0x6000u;
	case 16: return 0xD008u;
	default: return 0;
	}
}

/** Fill data->chips with +-1 and check the sequence really is maximal-length. */
static int pt_make_chips(struct aylp_prbs_test_data *data)
{
	unsigned mask = pt_lfsr_mask(data->order);
	if (!mask) {
		log_error("prbs_test: order %u is not in the tap table "
			"(5 to 16)", data->order);
		return -1;
	}
	unsigned state = 1;
	for (size_t i = 0; i < data->n_chips; i++) {
		data->chips[i] = (state & 1) ? 1 : -1;
		unsigned lsb = state & 1;
		state >>= 1;
		if (lsb) state ^= mask;
		if (state == 1 && i + 1 < data->n_chips) {
			log_error("prbs_test: order %u tap mask 0x%X closes "
				"after %zu chips, not %zu -- not a "
				"maximum-length sequence", data->order, mask,
				i + 1, data->n_chips);
			return -1;
		}
	}
	if (state != 1) {
		log_error("prbs_test: order %u tap mask 0x%X did not return to "
			"its seed after %zu chips", data->order, mask,
			data->n_chips);
		return -1;
	}
	return 0;
}

static int pt_cmp_dbl(const void *a, const void *b)
{
	double d = *(const double *)a - *(const double *)b;
	return (d > 0) - (d < 0);
}

/** Median of arr (sorted in place), and the median absolute deviation. */
static void pt_median_mad(double *arr, size_t n, double *med, double *mad)
{
	*med = *mad = 0.0;
	if (!n) return;
	qsort(arr, n, sizeof *arr, pt_cmp_dbl);
	*med = arr[n/2];
	double *d = xmalloc(n * sizeof *d);
	for (size_t i = 0; i < n; i++) d[i] = fabs(arr[i] - *med);
	qsort(d, n, sizeof *d, pt_cmp_dbl);
	*mad = d[n/2];
	xfree(d);
}

/** Cross-correlate one response window against the command, over lags
 * [-neg_lags, max_lag]. Every lag is evaluated over the same set of command
 * samples, so the correlogram's shape is exact rather than tapering at the
 * ends. The acausal lags cannot carry any real response, so their scatter is
 * the noise floor everything else is judged against.
 *
 * `live` marks which response samples the sensor actually measured; the rest
 * are values it held through a dark frame and are skipped. A held sample is
 * not missing data, it is WRONG data -- it carries a value from some earlier
 * frame, so correlating it drags the apparent lag toward that frame's age.
 * Skipping means each lag sees a different subset of the response, so each is
 * normalized over its own subset (a plain Pearson r) rather than all sharing
 * one. With every sample live that is the same correlogram as before, up to
 * the per-lag mean removal. `live` may be null, meaning every sample counts.
 *
 * Returns 0 if a peak was found, -1 if the correlation never leaves the noise.
 * rho may be null if only the summary numbers are wanted. */
static int pt_correlate(struct aylp_prbs_test_data *data, const double *r,
	const unsigned char *live, double *rho_out, double *peak_rho,
	double *lag_peak, double *lag_onset, double *lag_cent,
	double *noise_out, size_t *n_thin)
{
	size_t nl = data->neg_lags + data->max_lag + 1;
	size_t k0 = data->neg_lags, k1 = data->win - data->max_lag;
	const double *c = data->cmd_win;

	// a lag fit over a handful of surviving samples is noise dressed as a
	// correlation, so it is dropped rather than allowed to win the peak
	size_t span = k1 - k0;
	size_t need = (size_t)(data->min_pair_frac * (double)span);
	if (need < 8) need = 8;

	double *rho = rho_out ? rho_out : xmalloc(nl * sizeof *rho);
	unsigned char *fit = xcalloc(nl, sizeof *fit);
	size_t n_ok = 0, thin = 0;
	for (size_t i = 0; i < nl; i++) {
		long lag = (long)i - (long)data->neg_lags;
		double sc = 0, sr = 0, scc = 0, srr = 0, scr = 0;
		size_t n = 0;
		for (size_t k = k0; k < k1; k++) {
			size_t j = (size_t)((long)k + lag);
			if (live && !live[j]) continue;
			double cv = c[k], rv = r[j];
			sc += cv; sr += rv;
			scc += cv*cv; srr += rv*rv; scr += cv*rv;
			n++;
		}
		rho[i] = 0.0;
		if (n < need) { thin++; continue; }
		double dn = (double)n;
		double cov = scr - sc*sr/dn;
		double vc = scc - sc*sc/dn, vr = srr - sr*sr/dn;
		if (vc > 0.0 && vr > 0.0) {
			rho[i] = cov / sqrt(vc * vr);
			fit[i] = 1;
			n_ok++;
		} else {
			thin++;
		}
	}
	if (n_thin) *n_thin = thin;
	if (!n_ok) {
		xfree(fit);
		if (!rho_out) xfree(rho);
		return -1;
	}

	// noise floor: the acausal half, which by causality holds no response.
	// Lags that were dropped for want of live pairs are left out of it --
	// their rho is 0 by convention, not by measurement, and averaging those
	// zeros in would report a quieter floor than the run actually had.
	double nsum = 0.0;
	size_t nn = 0;
	for (size_t i = 0; i < data->neg_lags; i++) {
		if (!fit[i]) continue;
		nsum += rho[i]*rho[i];
		nn++;
	}
	double noise = nn ? sqrt(nsum / nn) : 0.0;

	// peak of |rho| over the causal lags that were actually fit
	size_t zero = data->neg_lags, ipk = nl;
	for (size_t i = zero; i < nl; i++)
		if (fit[i] && (ipk == nl || fabs(rho[i]) > fabs(rho[ipk])))
			ipk = i;
	if (ipk == nl) {
		xfree(fit);
		if (!rho_out) xfree(rho);
		return -1;
	}
	double sign = rho[ipk] < 0.0 ? -1.0 : 1.0;
	double peak = rho[ipk];
	int ret = -1;
	if (fabs(peak) > 4.0 * noise && fabs(peak) > 1e-9) {
		ret = 0;
		// parabolic refinement on the three points around the peak,
		// but only if both neighbours were fit -- a dropped lag reads
		// 0 and would pull the refinement toward itself
		double dl = 0.0;
		if (ipk > 0 && ipk + 1 < nl && fit[ipk-1] && fit[ipk+1]) {
			double ym = sign*rho[ipk-1], y0 = sign*rho[ipk];
			double yp = sign*rho[ipk+1];
			double den = ym - 2.0*y0 + yp;
			if (fabs(den) > 1e-15) dl = 0.5 * (ym - yp) / den;
			if (dl > 1.0 || dl < -1.0) dl = 0.0;
		}
		*lag_peak = (double)ipk - (double)zero + dl;

		// onset: first crossing of max(onset_frac*peak, 4 sigma),
		// interpolated between the two straddling lags
		double thr = data->onset_frac * fabs(peak);
		if (4.0 * noise > thr) thr = 4.0 * noise;
		*lag_onset = *lag_peak;
		for (size_t i = zero; i <= ipk; i++) {
			if (sign*rho[i] < thr) continue;
			double lag = (double)i - (double)zero;
			if (i > zero) {
				double y0 = sign*rho[i-1], y1 = sign*rho[i];
				if (y1 > y0) lag -= (y1 - thr) / (y1 - y0);
			}
			*lag_onset = lag;
			break;
		}

		// centroid of the correlation lobe around the peak
		double wsum = 0.0, lsum = 0.0, cut = 0.1 * fabs(peak);
		for (size_t i = ipk; i < nl && sign*rho[i] > cut; i++) {
			wsum += sign*rho[i];
			lsum += sign*rho[i] * ((double)i - (double)zero);
		}
		for (size_t i = ipk; i-- > zero && sign*rho[i] > cut; ) {
			wsum += sign*rho[i];
			lsum += sign*rho[i] * ((double)i - (double)zero);
		}
		*lag_cent = wsum > 0.0 ? lsum / wsum : *lag_peak;
	}
	*peak_rho = peak;
	*noise_out = noise;
	xfree(fit);
	if (!rho_out) xfree(rho);
	return ret;
}

/** Group delay from the slope of the cross-spectrum's phase, weighted by
 * cross-spectral magnitude, over [phase_f_lo, phase_f_hi]. Same estimator and
 * same band as bode_plot's fit_tau_ms, so the two are directly comparable. */
static void pt_phase_slope(struct aylp_prbs_test_data *data, double sign)
{
	size_t nfft = 1;
	while (nfft < data->win) nfft <<= 1;
	if (!(data->fs > 0.0)) return;

	double *cf = xcalloc(nfft, sizeof *cf);
	double *rf = xcalloc(nfft, sizeof *rf);
	double rmean = 0.0;
	for (size_t k = 0; k < data->win; k++) rmean += data->rbar[k];
	rmean /= data->win;
	for (size_t k = 0; k < data->win; k++) {
		cf[k] = data->cmd_win[k];
		rf[k] = data->rbar[k] - rmean;
	}
	gsl_fft_real_radix2_transform(cf, 1, nfft);
	gsl_fft_real_radix2_transform(rf, 1, nfft);

	// halfcomplex layout: re(k) = a[k], im(k) = a[nfft-k]
	double sw = 0, swf = 0, swf2 = 0, swp = 0, swfp = 0;
	size_t n = 0;
	// The correlation peak is an independent, time-domain delay estimate.
	// Use it to select the phase branch at every frequency independently.
	// Sequential unwrapping lets one corrupted bin (notably a PWM-lighting
	// harmonic) cycle-slip every later bin by 2*pi.
	double tau_seed = data->lag_peak / data->fs;
	double *fs_ = xmalloc((nfft/2) * sizeof *fs_);
	double *ph_ = xmalloc((nfft/2) * sizeof *ph_);
	double *w_ = xmalloc((nfft/2) * sizeof *w_);
	for (size_t k = 1; k < nfft/2; k++) {
		double f = k * data->fs / nfft;
		if (f < data->phase_f_lo || f > data->phase_f_hi) continue;
		if (data->phase_notch_hz > 0.0 &&
			data->phase_notch_width_hz > 0.0) {
			double harmonic = round(f / data->phase_notch_hz);
			if (harmonic >= 1.0 && fabs(f - harmonic *
				data->phase_notch_hz) <= data->phase_notch_width_hz)
				continue;
		}
		double cr = cf[k], ci = -cf[nfft-k];	// conj(C)
		double rr = rf[k], ri = rf[nfft-k];
		double sr = rr*cr - ri*ci;
		double si = rr*ci + ri*cr;
		double mag = hypot(sr, si);
		if (!(mag > 0.0)) continue;
		double phi = atan2(si, sr);
		// a plant with negative DC gain sits at pi; take it out before
		// unwrapping or the fit chases a constant offset
		if (sign < 0.0) phi += M_PI;
		double expected = -2.0*M_PI * f * tau_seed;
		phi += 2.0*M_PI * round((expected - phi) / (2.0*M_PI));
		fs_[n] = f;
		ph_[n] = phi;
		w_[n] = mag;
		sw += mag; swf += mag*f; swf2 += mag*f*f;
		swp += mag*phi; swfp += mag*f*phi;
		n++;
	}
	double den = sw*swf2 - swf*swf;
	if (n >= 4 && fabs(den) > 0.0) {
		double slope = (sw*swfp - swf*swp) / den;
		double icept = (swf2*swp - swf*swfp) / den;
		// phi = -2 pi f tau, so tau is -slope / 2 pi
		data->tau_phase_ms = -1e3 * slope / (2.0*M_PI);
		double ss = 0.0;
		for (size_t i = 0; i < n; i++) {
			double r = ph_[i] - (icept + slope*fs_[i]);
			ss += w_[i] * r * r;
		}
		data->phase_resid_deg = 180.0/M_PI * sqrt(ss / sw);
		data->phase_n = n;
	}
	xfree(fs_);
	xfree(ph_);
	xfree(w_);
	xfree(cf);
	xfree(rf);
}

static int pt_write_dat(struct aylp_prbs_test_data *data, const char *path,
	const char *traces_path)
{
	FILE *f = fopen(path, "w");
	if (!f) {
		log_error("prbs_test: fopen %s: %m", path);
		return -1;
	}
	time_t tt = time(NULL);
	struct tm tm;
	localtime_r(&tt, &tm);
	double dt_ms = data->fs > 0.0 ? 1e3 / data->fs : 0.0;
	fprintf(f, "# prbs_test %04d-%02d-%02d %02d:%02d:%02d\n",
		tm.tm_year+1900, tm.tm_mon+1, tm.tm_mday,
		tm.tm_hour, tm.tm_min, tm.tm_sec);
	fprintf(f, "# index_cmd %d index_err %d amplitude %G bias %G "
		"order %u chips %zu chip_frames %zu bursts %zu quiet %zu\n",
		data->index_cmd, data->index_err, data->amplitude, data->bias,
		data->order, data->n_chips, data->chip_frames, data->n_bursts,
		data->quiet_frames);
	fprintf(f, "# fs_measured %.2f dt_ms %.4f jitter_rms_ms %.4f "
		"volts_per_unit %G pixel_scale %G\n", data->fs, dt_ms,
		1e3 * (data->dt_cnt > 1
			? sqrt(fmax(0.0, data->dt_sum2/data->dt_cnt
				- (data->dt_sum/data->dt_cnt)
					* (data->dt_sum/data->dt_cnt)))
			: 0.0),
		data->volts_per_unit, data->pixel_scale);
	fprintf(f, "# sensor_live %.4f (%zu of %zu frames) holes %zu "
		"min_live_bursts %zu thin_lags %zu use_rejected %d\n",
		data->n_frames
			? (double)data->n_live / (double)data->n_frames : 0.0,
		data->n_live, data->n_frames, data->n_holes, data->worst_cnt,
		data->n_thin_lags, data->use_rejected ? 1 : 0);
	if (data->have_result) {
		fprintf(f, "# onset %.3f fr = %.4f ms\n",
			data->lag_onset, data->lag_onset * dt_ms);
		fprintf(f, "# peak %.3f fr = %.4f ms (rho %.4f, noise rms "
			"%.5f, peak/noise %.1f)\n", data->lag_peak,
			data->lag_peak * dt_ms, data->rho_peak,
			data->rho_noise, data->rho_noise > 0.0
				? fabs(data->rho_peak)/data->rho_noise : 0.0);
		fprintf(f, "# centroid %.3f fr = %.4f ms\n",
			data->lag_centroid, data->lag_centroid * dt_ms);
		fprintf(f, "# recommended_delay %.3f fr source %s\n",
			data->use_correlation_delay ? data->lag_centroid :
				(data->tau_phase_ms / dt_ms),
			data->use_correlation_delay ? "correlation_centroid" :
				"phase_slope");
		if (data->phase_n)
			fprintf(f, "# phase_slope tau %.4f ms = %.3f fr over "
				"%G-%G Hz (%zu bins, residual %.2f deg)\n",
				data->tau_phase_ms,
				dt_ms > 0.0 ? data->tau_phase_ms/dt_ms : 0.0,
				data->phase_f_lo, data->phase_f_hi,
				data->phase_n, data->phase_resid_deg);
		fprintf(f, "# per_burst onset %.3f +/- %.3f fr  peak %.3f "
			"+/- %.3f fr (median +/- MAD over %zu of %zu bursts)\n",
			data->lag_onset_med, data->lag_onset_mad,
			data->lag_peak_med, data->lag_peak_mad,
			data->n_burst_ok, data->n_bursts);
		fprintf(f, "# quality %s%s%s\n",
			data->quality_ok ? "PASS" : "FAIL",
			data->quality_ok ? "" : ": ",
			data->quality_ok ? "" : data->quality_reason);
	} else {
		fprintf(f, "# NO RESULT: correlation never left the noise\n");
	}
	if (data->config) fprintf(f, "# config %s\n", data->config);
	fprintf(f, "# lag_frames lag_ms rho\n");
	size_t nl = data->neg_lags + data->max_lag + 1;
	for (size_t i = 0; i < nl; i++) {
		double lag = (double)i - (double)data->neg_lags;
		fprintf(f, "%g %g %g\n", lag, lag * dt_ms, data->rho[i]);
	}
	fclose(f);
	log_info("prbs_test: correlogram saved to %s", path);

	f = fopen(traces_path, "w");
	if (!f) {
		log_error("prbs_test: fopen %s: %m", traces_path);
		return -1;
	}
	fprintf(f, "# prbs_test ensemble-averaged traces over %zu bursts\n",
		data->n_bursts);
	fprintf(f, "# resp_n is how many bursts were live at that position; 0 "
		"means resp_mean is interpolated, not measured\n");
	fprintf(f, "# frame t_ms command resp_mean resp_sem resp_n\n");
	for (size_t k = 0; k < data->win; k++)
		fprintf(f, "%zu %g %g %g %g %zu\n", k, k * dt_ms,
			data->cmd_win[k], data->rbar[k], data->rsem[k],
			data->rcnt[k]);
	fclose(f);
	log_info("prbs_test: averaged traces saved to %s", traces_path);
	return 0;
}

static int pt_plot(struct aylp_prbs_test_data *data, const char *dat_path,
	const char *traces_path)
{
	char script_path[64];
	snprintf(script_path, sizeof script_path, "/tmp/prbs_test_%d.py",
		getpid());
	FILE *script = fopen(script_path, "w");
	if (!script) {
		log_error("prbs_test: fopen %s: %m", script_path);
		return -1;
	}
	fputs(PLOT_SCRIPT, script);
	fclose(script);

	double dt_ms = data->fs > 0.0 ? 1e3 / data->fs : 0.0;
	char title[512];
	snprintf(title, sizeof title,
		"onset %.2f fr (%.3f ms)  peak %.2f fr (%.3f ms)  phase-slope "
		"%.3f ms (%.2f fr)  rho %.3f, peak/noise %.0f  fs %.1f Hz",
		data->lag_onset, data->lag_onset * dt_ms, data->lag_peak,
		data->lag_peak * dt_ms, data->tau_phase_ms,
		dt_ms > 0.0 ? data->tau_phase_ms/dt_ms : 0.0, data->rho_peak,
		data->rho_noise > 0.0
			? fabs(data->rho_peak)/data->rho_noise : 0.0, data->fs);

	// show the burst plus a little of the quiet tail, not the whole window
	size_t nshow = data->burst_frames + data->max_lag + 8;
	if (nshow > data->win) nshow = data->win;

	char cmd[1024];
	snprintf(cmd, sizeof cmd, "python3 '%s' '%s' '%s' '%s' '%s' %G %G %G %zu",
		script_path, dat_path, traces_path, data->output_file, title,
		data->lag_onset, data->lag_peak, data->pixel_scale, nshow);
	int ret = system(cmd);
	int ok = (ret != -1) && WIFEXITED(ret) && !WEXITSTATUS(ret);
	if (!ok)
		log_error("prbs_test: plot script failed (exit status %d); is "
			"python3 with numpy and matplotlib installed?",
			(ret != -1 && WIFEXITED(ret)) ? WEXITSTATUS(ret) : ret);
	else
		log_info("prbs_test: plot saved to %s", data->output_file);
	unlink(script_path);
	return ok ? 0 : -1;
}

/** Fill window positions no burst was live at, by linear interpolation between
 * the nearest live positions on either side (edge-held at the ends). Only the
 * phase slope and the plot use these: `rlive` still marks them dead, so the
 * correlogram never sees an invented sample. */
static void pt_fill_holes(struct aylp_prbs_test_data *data)
{
	size_t win = data->win;
	size_t prev = win;	// last live position seen, win = none yet
	for (size_t k = 0; k < win; k++) {
		if (data->rlive[k]) { prev = k; continue; }
		size_t next = k;
		while (next < win && !data->rlive[next]) next++;
		double a = prev < win ? data->rbar[prev] : 0.0;
		double b = next < win ? data->rbar[next] : a;
		if (prev >= win) a = b;
		for (size_t j = k; j < next; j++) {
			double f = next > prev && prev < win
				? (double)(j - prev) / (double)(next - prev)
				: 0.0;
			data->rbar[j] = a + f * (b - a);
		}
		k = next - 1;
	}
}

static void pt_quality_fail(struct aylp_prbs_test_data *data, size_t *used,
	const char *fmt, ...)
{
	size_t cap = sizeof data->quality_reason;
	if (*used >= cap - 1) return;
	if (*used) {
		data->quality_reason[(*used)++] = ';';
		if (*used < cap - 1) data->quality_reason[(*used)++] = ' ';
	}
	va_list ap;
	va_start(ap, fmt);
	int nw = vsnprintf(data->quality_reason + *used, cap - *used, fmt, ap);
	va_end(ap);
	if (nw <= 0) return;
	size_t add = (size_t)nw;
	*used += add < cap - *used ? add : cap - *used - 1;
}

/** Average the bursts, correlate, fit, log, and write the files. */
static void pt_analyze(struct aylp_prbs_test_data *data)
{
	size_t win = data->win, nb = data->burst;
	if (!nb) {
		log_error("prbs_test: no bursts recorded");
		return;
	}
	data->fs = data->dt_cnt ? data->dt_cnt / data->dt_sum : 0.0;

	// Ensemble mean and its standard error, frame by frame, over the bursts
	// in which the sensor was live at that window position. This is where
	// most of the repair happens: the chop is not locked to the burst
	// period, so a position that was dark in one burst is nearly always
	// live in others, and the mean is built from real samples rather than
	// from a mixture of real ones and stale ones.
	data->worst_cnt = nb;
	data->n_holes = 0;
	for (size_t k = 0; k < win; k++) {
		double s = 0.0, s2 = 0.0;
		size_t n = 0;
		for (size_t b = 0; b < nb; b++) {
			if (!data->use_rejected && !data->live[b*win + k])
				continue;
			double v = data->resp[b*win + k];
			s += v;
			s2 += v * v;
			n++;
		}
		data->rcnt[k] = n;
		data->rlive[k] = n > 0;
		if (n < data->worst_cnt) data->worst_cnt = n;
		if (!n) {
			data->n_holes++;
			data->rbar[k] = 0.0;
			data->rsem[k] = 0.0;
			continue;
		}
		double m = s / n;
		double var = s2/n - m*m;
		if (var < 0.0) var = 0.0;
		data->rbar[k] = m;
		data->rsem[k] = n > 1 ? sqrt(var / (n - 1)) : 0.0;
	}
	// A position no burst ever saw is a hole. The correlator skips it, but
	// the phase slope is an FFT and cannot, so fill it by interpolating its
	// live neighbours -- and say how many, because an interpolated sample
	// is a guess about the very quantity being measured.
	if (data->n_holes) pt_fill_holes(data);

	double lag_cent = 0.0;
	int ok = pt_correlate(data, data->rbar, data->rlive, data->rho,
		&data->rho_peak, &data->lag_peak, &data->lag_onset, &lag_cent,
		&data->rho_noise, &data->n_thin_lags);
	data->lag_centroid = lag_cent;
	data->have_result = !ok;
	data->quality_ok = false;
	data->quality_reason[0] = '\0';

	if (!data->have_result) {
		log_error("prbs_test: the correlation never left its noise "
			"floor (peak rho %.5G against noise rms %.5G). Raise "
			"amplitude or n_bursts, or check that the actuator is "
			"actually moving the beam.", data->rho_peak,
			data->rho_noise);
	} else {
		pt_phase_slope(data, data->rho_peak < 0.0 ? -1.0 : 1.0);

		// per-burst spread
		double *pk = xmalloc(nb * sizeof *pk);
		double *on = xmalloc(nb * sizeof *on);
		size_t n = 0;
		for (size_t b = 0; b < nb; b++) {
			double p, o, c, rp, nz;
			const unsigned char *lv = data->use_rejected
				? NULL : data->live + b*win;
			if (pt_correlate(data, data->resp + b*win, lv, NULL,
					&rp, &p, &o, &c, &nz, NULL))
				continue;
			pk[n] = p;
			on[n] = o;
			n++;
		}
		data->n_burst_ok = n;
		pt_median_mad(pk, n, &data->lag_peak_med, &data->lag_peak_mad);
		pt_median_mad(on, n, &data->lag_onset_med,
			&data->lag_onset_mad);
		xfree(pk);
		xfree(on);

		// A correlation peak is only a candidate result. Calibration runs can
		// demand enough SNR, phase linearity, independent-burst support and
		// live coverage to make that candidate safe to copy into a controller.
		// Keep every failed condition in one message so a bench run says what
		// needs fixing instead of making the operator play whack-a-mole.
		size_t used = 0;
#define PT_QUALITY_FAIL(...) \
	pt_quality_fail(data, &used, __VA_ARGS__)
		double snr = data->rho_noise > 0.0
			? fabs(data->rho_peak) / data->rho_noise : INFINITY;
		double live_frac = data->n_frames
			? (double)data->n_live / data->n_frames : 0.0;
		double burst_frac = nb ? (double)data->n_burst_ok / nb : 0.0;
		double phase_frames = data->fs > 0.0
			? data->tau_phase_ms * data->fs / 1e3 : 0.0;
		if (data->min_peak_snr > 0.0 && snr < data->min_peak_snr)
			PT_QUALITY_FAIL("peak/noise %.2f < %.2f", snr,
				data->min_peak_snr);
		if (!data->phase_n && !data->use_correlation_delay)
			PT_QUALITY_FAIL("no phase-slope fit");
		else if (!data->use_correlation_delay) {
			if (data->max_phase_resid_deg > 0.0
					&& data->phase_resid_deg
					> data->max_phase_resid_deg)
				PT_QUALITY_FAIL("phase residual %.2f deg > %.2f deg",
					data->phase_resid_deg,
					data->max_phase_resid_deg);
			if (data->min_phase_bins
					&& data->phase_n < data->min_phase_bins)
				PT_QUALITY_FAIL("phase bins %zu < %zu", data->phase_n,
					data->min_phase_bins);
			if (data->max_phase_peak_delta > 0.0
					&& fabs(phase_frames - data->lag_peak)
					> data->max_phase_peak_delta)
				PT_QUALITY_FAIL("phase/peak delta %.3f fr > %.3f fr",
					fabs(phase_frames - data->lag_peak),
					data->max_phase_peak_delta);
		}
		if (data->min_burst_frac > 0.0
				&& burst_frac < data->min_burst_frac)
			PT_QUALITY_FAIL("usable bursts %.1f%% < %.1f%%",
				100.0 * burst_frac, 100.0 * data->min_burst_frac);
		if (data->max_peak_mad > 0.0
				&& data->lag_peak_mad > data->max_peak_mad)
			PT_QUALITY_FAIL("peak MAD %.3f fr > %.3f fr",
				data->lag_peak_mad, data->max_peak_mad);
		if (data->min_live_frac > 0.0 && live_frac < data->min_live_frac)
			PT_QUALITY_FAIL("live frames %.1f%% < %.1f%%",
				100.0 * live_frac, 100.0 * data->min_live_frac);
		if (data->n_holes > data->max_holes)
			PT_QUALITY_FAIL("holes %zu > %zu", data->n_holes,
				data->max_holes);
		data->quality_ok = !used;
#undef PT_QUALITY_FAIL
	}

	double dt_ms = data->fs > 0.0 ? 1e3 / data->fs : 0.0;
	double jit = data->dt_cnt > 1
		? sqrt(fmax(0.0, data->dt_sum2/data->dt_cnt
			- (data->dt_sum/data->dt_cnt)
				* (data->dt_sum/data->dt_cnt))) : 0.0;
	log_info("prbs_test: ================ RESULTS ================");
	log_info("prbs_test: %zu bursts of %zu chips at %zu frame(s)/chip; "
		"loop ran at %.1f Hz (%.4f ms/frame, jitter rms %.4f ms)",
		nb, data->n_chips, data->chip_frames, data->fs, dt_ms,
		1e3 * jit);
	if (jit > 0.2 / (data->fs > 0.0 ? data->fs : 1.0))
		log_warn("prbs_test: frame interval jitter is over 20%% of the "
			"period; lags in frames are still exact but their "
			"conversion to ms is smeared");
	double lf = data->n_frames
		? (double)data->n_live / (double)data->n_frames : 0.0;
	log_info("prbs_test: sensor live on %.1f%% of frames (%zu of %zu); "
		"%zu window position(s) had no live sample in any burst, "
		"fewest live bursts at any position %zu of %zu",
		1e2 * lf, data->n_live, data->n_frames, data->n_holes,
		data->worst_cnt, nb);
	if (data->use_rejected)
		log_warn("prbs_test: use_rejected is set, so samples the sensor "
			"flagged were correlated as if they were data. A held "
			"sample carries an older frame's value, which biases "
			"the lag LONG; this run is for comparison, not for a "
			"number to configure a loop with");
	else if (lf < 0.999)
		log_info("prbs_test: those %zu held samples were excluded from "
			"the ensemble mean and from every correlation, so the "
			"lags below are fit only where the beam was actually "
			"on", data->n_frames - data->n_live);
	if (!data->use_rejected && lf < 0.5)
		log_warn("prbs_test: the sensor was live on under half the "
			"frames -- the ensemble mean is thin and the noise "
			"floor is correspondingly higher. Fix the source "
			"before trusting this");
	if (data->n_holes)
		log_warn("prbs_test: %zu window position(s) were dark in EVERY "
			"burst, so their ensemble mean is interpolated from "
			"neighbours. The correlogram skips them, but the phase "
			"slope cannot -- treat tau as approximate. A chop "
			"locked to the burst period does exactly this; change "
			"quiet_frames to walk the two apart",
			data->n_holes);
	if (data->n_thin_lags)
		log_warn("prbs_test: %zu lag(s) had too few live sample pairs "
			"to fit (min_pair_frac %G) and were dropped from the "
			"correlogram", data->n_thin_lags, data->min_pair_frac);
	if (!data->have_result) goto write;
	if (!data->quality_ok) {
		log_error("prbs_test: QUALITY GATES FAILED: %s. Do not copy this "
			"delay into a controller.", data->quality_reason);
		goto write;
	}
	log_info("prbs_test: QUALITY GATES PASSED");
	log_info("prbs_test: RECOMMENDED DELAY %.3f frames from %s",
		data->use_correlation_delay ? data->lag_centroid :
			(data->tau_phase_ms / dt_ms),
		data->use_correlation_delay ? "dropout-safe correlation centroid" :
			"phase slope");
	log_info("prbs_test: ONSET (transport delay) %.3f frames = %.4f ms",
		data->lag_onset, data->lag_onset * dt_ms);
	log_info("prbs_test: PEAK (delay + plant rise) %.3f frames = %.4f ms; "
		"correlation %.4f against a %.5f noise floor (%.0fx)",
		data->lag_peak, data->lag_peak * dt_ms, data->rho_peak,
		data->rho_noise, data->rho_noise > 0.0
			? fabs(data->rho_peak)/data->rho_noise : 0.0);
	log_info("prbs_test: centroid of the lobe %.3f frames = %.4f ms",
		data->lag_centroid, data->lag_centroid * dt_ms);
	if (data->phase_n)
		log_info("prbs_test: PHASE SLOPE tau %.4f ms = %.3f frames "
			"over %G-%G Hz (%zu bins, residual %.2f deg) -- this "
			"is the bode_plot-comparable number, and the one for "
			"an fsp delay/delay_frac pair", data->tau_phase_ms,
			dt_ms > 0.0 ? data->tau_phase_ms/dt_ms : 0.0,
			data->phase_f_lo, data->phase_f_hi, data->phase_n,
			data->phase_resid_deg);
	log_info("prbs_test: per-burst scatter: onset %.3f +/- %.3f fr, peak "
		"%.3f +/- %.3f fr (median +/- MAD over %zu of %zu bursts)",
		data->lag_onset_med, data->lag_onset_mad, data->lag_peak_med,
		data->lag_peak_mad, data->n_burst_ok, nb);
	if (data->lag_peak < 1.0)
		log_warn("prbs_test: a peak below 1 frame is not physical -- "
			"the sensor value read in an iteration was captured "
			"before that iteration's command went out. Treat this "
			"run as a failure, not as a fast loop.");

write:
	if (!data->output_file) return;
	const char *out = data->output_file;
	const char *dot = strrchr(out, '.');
	size_t base = dot ? (size_t)(dot - out) : strlen(out);
	char *dat_path = xmalloc(base + 5);
	memcpy(dat_path, out, base);
	memcpy(dat_path + base, ".dat", 5);
	char *traces_path = xmalloc(base + 12);
	memcpy(traces_path, out, base);
	memcpy(traces_path + base, "_traces.dat", 12);
	if (!pt_write_dat(data, dat_path, traces_path)
			&& data->have_result && data->quality_ok)
		pt_plot(data, dat_path, traces_path);
	xfree(dat_path);
	xfree(traces_path);
}


int prbs_test_init(struct aylp_device *self)
{
	self->proc = &prbs_test_proc;
	self->fini = &prbs_test_fini;
	struct aylp_prbs_test_data *data = xcalloc(1, sizeof *data);
	self->device_data = data;

	// defaults
	data->index_cmd = 1;
	data->index_err = 1;
	data->out_size = 2;
	data->amplitude = 0.05;
	data->order = 7;
	data->chip_frames = 1;
	data->n_bursts = 32;
	data->quiet_frames = 128;
	data->warmup = 5.0;
	data->max_lag = 48;
	data->neg_lags = 24;
	data->onset_frac = 0.2;
	data->phase_f_lo = 20.0;
	data->phase_f_hi = 200.0;
	data->volts_per_unit = 1.0;
	data->pixel_scale = 1.0;
	data->min_pair_frac = 0.25;
	data->max_holes = SIZE_MAX;
	data->output_file = xstrdup("prbs_test.pdf");

	if (self->params) { json_object_object_foreach(self->params, key, val) {
		if (key[0] == '_') {
		} else if (!strcmp(key, "index_cmd")) {
			data->index_cmd = (int)json_object_get_int64(val);
		} else if (!strcmp(key, "index_err")) {
			data->index_err = (int)json_object_get_int64(val);
		} else if (!strcmp(key, "out_size")) {
			data->out_size = (size_t)json_object_get_uint64(val);
		} else if (!strcmp(key, "amplitude")) {
			data->amplitude = json_object_get_double(val);
		} else if (!strcmp(key, "bias")) {
			data->bias = json_object_get_double(val);
		} else if (!strcmp(key, "order")) {
			data->order = (unsigned)json_object_get_uint64(val);
		} else if (!strcmp(key, "chip_frames")) {
			data->chip_frames = (size_t)json_object_get_uint64(val);
		} else if (!strcmp(key, "n_bursts")) {
			data->n_bursts = (size_t)json_object_get_uint64(val);
		} else if (!strcmp(key, "quiet_frames")) {
			data->quiet_frames = (size_t)json_object_get_uint64(val);
		} else if (!strcmp(key, "warmup")) {
			data->warmup = json_object_get_double(val);
		} else if (!strcmp(key, "max_lag")) {
			data->max_lag = (size_t)json_object_get_uint64(val);
		} else if (!strcmp(key, "neg_lags")) {
			data->neg_lags = (size_t)json_object_get_uint64(val);
		} else if (!strcmp(key, "onset_frac")) {
			data->onset_frac = json_object_get_double(val);
		} else if (!strcmp(key, "phase_f_lo")) {
			data->phase_f_lo = json_object_get_double(val);
		} else if (!strcmp(key, "phase_f_hi")) {
			data->phase_f_hi = json_object_get_double(val);
		} else if (!strcmp(key, "phase_notch_hz")) {
			data->phase_notch_hz = json_object_get_double(val);
		} else if (!strcmp(key, "phase_notch_width_hz")) {
			data->phase_notch_width_hz = json_object_get_double(val);
		} else if (!strcmp(key, "use_correlation_delay")) {
			data->use_correlation_delay = json_object_get_boolean(val);
		} else if (!strcmp(key, "volts_per_unit")) {
			data->volts_per_unit = json_object_get_double(val);
		} else if (!strcmp(key, "pixel_scale")) {
			data->pixel_scale = json_object_get_double(val);
		} else if (!strcmp(key, "use_rejected")) {
			data->use_rejected = json_object_get_boolean(val);
		} else if (!strcmp(key, "min_pair_frac")) {
			data->min_pair_frac = json_object_get_double(val);
		} else if (!strcmp(key, "min_peak_snr")) {
			data->min_peak_snr = json_object_get_double(val);
		} else if (!strcmp(key, "max_phase_resid_deg")) {
			data->max_phase_resid_deg = json_object_get_double(val);
		} else if (!strcmp(key, "min_phase_bins")) {
			data->min_phase_bins = json_object_get_uint64(val);
		} else if (!strcmp(key, "min_burst_frac")) {
			data->min_burst_frac = json_object_get_double(val);
		} else if (!strcmp(key, "max_peak_mad")) {
			data->max_peak_mad = json_object_get_double(val);
		} else if (!strcmp(key, "min_live_frac")) {
			data->min_live_frac = json_object_get_double(val);
		} else if (!strcmp(key, "max_holes")) {
			data->max_holes = json_object_get_uint64(val);
		} else if (!strcmp(key, "max_phase_peak_delta")) {
			data->max_phase_peak_delta = json_object_get_double(val);
		} else if (!strcmp(key, "output_file")) {
			// an empty name means "log the numbers, write nothing"
			const char *v = json_object_get_string(val);
			xfree(data->output_file);
			if (v && *v) data->output_file = xstrdup(v);
		} else if (!strcmp(key, "config")) {
			data->config = xstrdup(json_object_get_string(val));
		} else {
			log_warn("prbs_test: unknown param \"%s\"", key);
		}
	} }

	if ((size_t)data->index_cmd >= data->out_size) {
		log_error("prbs_test: index_cmd %d not < out_size %zu",
			data->index_cmd, data->out_size);
		return -1;
	}
	if (data->amplitude == 0.0 || !data->chip_frames || !data->n_bursts) {
		log_error("prbs_test: need a nonzero amplitude and at least "
			"one frame per chip and one burst");
		return -1;
	}
	if (data->onset_frac <= 0.0 || data->onset_frac >= 1.0) {
		log_error("prbs_test: onset_frac must be in (0, 1)");
		return -1;
	}
	if (data->min_pair_frac <= 0.0 || data->min_pair_frac > 1.0) {
		log_error("prbs_test: min_pair_frac must be in (0, 1]");
		return -1;
	}
	if (data->min_peak_snr < 0.0 || data->max_phase_resid_deg < 0.0
			|| data->min_burst_frac < 0.0 || data->min_burst_frac > 1.0
			|| data->max_peak_mad < 0.0 || data->min_live_frac < 0.0
			|| data->min_live_frac > 1.0
			|| data->max_phase_peak_delta < 0.0) {
		log_error("prbs_test: invalid quality-gate parameter");
		return -1;
	}
	if (data->order < 5 || data->order > 16) {
		log_error("prbs_test: order %u out of range (5 to 16)",
			data->order);
		return -1;
	}

	data->n_chips = ((size_t)1 << data->order) - 1;
	data->burst_frames = data->n_chips * data->chip_frames;
	data->win = data->burst_frames + data->quiet_frames;
	// every lag has to be evaluated over the same stretch of command, so
	// the window has to be longer than the lags it is asked about
	if (data->win < data->neg_lags + data->max_lag + 32) {
		log_error("prbs_test: window of %zu frames is too short for "
			"lags -%zu..%zu; lengthen the burst (order, "
			"chip_frames) or quiet_frames", data->win,
			data->neg_lags, data->max_lag);
		return -1;
	}
	if (data->quiet_frames < data->max_lag)
		log_warn("prbs_test: quiet_frames %zu is shorter than max_lag "
			"%zu, so the tail of one burst's response overlaps the "
			"start of the next", data->quiet_frames, data->max_lag);

	data->chips = xmalloc(data->n_chips * sizeof *data->chips);
	if (pt_make_chips(data)) return -1;

	data->out = xcalloc_type(gsl_vector, data->out_size);
	data->cmd_win = xcalloc(data->win, sizeof *data->cmd_win);
	data->resp = xcalloc(data->n_bursts * data->win, sizeof *data->resp);
	data->live = xcalloc(data->n_bursts * data->win, sizeof *data->live);
	data->rbar = xcalloc(data->win, sizeof *data->rbar);
	data->rsem = xcalloc(data->win, sizeof *data->rsem);
	data->rcnt = xcalloc(data->win, sizeof *data->rcnt);
	data->rlive = xcalloc(data->win, sizeof *data->rlive);
	data->rho = xcalloc(data->neg_lags + data->max_lag + 1,
		sizeof *data->rho);

	// the command of one whole window, with its mean removed once here so
	// the correlator never has to think about the burst's own DC
	double mean = 0.0;
	for (size_t k = 0; k < data->burst_frames; k++) {
		data->cmd_win[k] =
			data->amplitude * data->chips[k / data->chip_frames];
		mean += data->cmd_win[k];
	}
	mean /= data->win;
	for (size_t k = 0; k < data->win; k++) data->cmd_win[k] -= mean;

	data->stage = PT_WARMUP;
	data->cmd = data->bias;

	log_info("prbs_test: order %u = %zu chips at %zu frame(s)/chip = %zu "
		"frames per burst, +%zu quiet; %zu bursts, +-%G about %G",
		data->order, data->n_chips, data->chip_frames,
		data->burst_frames, data->quiet_frames, data->n_bursts,
		fabs(data->amplitude), data->bias);
	if (data->volts_per_unit != 1.0)
		log_info("prbs_test: command unit is %G V at the DAC, so the "
			"drive is +-%G V about %G V", data->volts_per_unit,
			fabs(data->amplitude * data->volts_per_unit),
			data->bias * data->volts_per_unit + 0.0);

	self->type_in = AYLP_T_VECTOR;
	self->units_in = AYLP_U_ANY;
	self->type_out = AYLP_T_VECTOR;
	self->units_out = AYLP_U_MINMAX;
	return 0;
}


int prbs_test_proc(struct aylp_device *self, struct aylp_state *state)
{
	struct aylp_prbs_test_data *data = self->device_data;
	double t = pt_now();
	int ret = 0;

	gsl_vector *s = state->vector;
	if (UNLIKELY((size_t)data->index_err >= s->size)) {
		log_error("prbs_test: index_err %d out of range (input size "
			"%zu)", data->index_err, s->size);
		ret = -1;
		goto fatal;
	}
	double err = s->data[data->index_err * s->stride];

	if (UNLIKELY(!data->t0)) {
		data->t0 = t;
		data->t_prev = t;
	}

	switch (data->stage) {
	case PT_WARMUP:
		data->cmd = data->bias;
		if (t - data->t0 < data->warmup) break;
		data->stage = PT_BURST;
		data->frame = 0;
		data->burst = 0;
		log_info("prbs_test: driving %zu bursts", data->n_bursts);
		break;

	case PT_BURST:
	case PT_QUIET: {
		// The response recorded at window-frame k is the one the sensor
		// delivered to THIS iteration, and the command chosen just
		// below is the one this same iteration hands to the DAC. So a
		// correlation lag of L frames is L whole iterations from
		// command to response, with nothing hidden in between.
		size_t k = data->frame;
		size_t idx = data->burst * data->win + k;
		data->resp[idx] = err;
		// Whether the beam was on for THIS frame, taken from the
		// sensor's own flags rather than from any assumed duty cycle.
		// fit_com, wfs_com and center_of_mass all re-publish their
		// previous coordinate when they cannot measure one, so without
		// this the correlator cannot tell a held value from a fresh
		// one -- and a held value is an older frame's, which drags the
		// apparent lag long.
		bool live = !(state->header.status
			& (AYLP_FRAME_REJECTED | AYLP_BEAM_LOST));
		data->live[idx] = live;
		data->n_frames++;
		data->n_live += live;
		double dt = t - data->t_prev;
		data->dt_sum += dt;
		data->dt_sum2 += dt * dt;
		data->dt_cnt++;

		if (k < data->burst_frames) {
			data->stage = PT_BURST;
			data->cmd = data->bias + data->amplitude
				* data->chips[k / data->chip_frames];
		} else {
			data->stage = PT_QUIET;
			data->cmd = data->bias;
		}

		data->frame++;
		if (data->frame >= data->win) {
			data->burst++;
			data->frame = 0;
			data->stage = PT_BURST;
			if (!(data->burst % 8) || data->burst == data->n_bursts)
				log_info("prbs_test: %zu/%zu bursts done",
					data->burst, data->n_bursts);
			if (data->burst >= data->n_bursts) {
				// park before analysing: writing files and
				// running a plot script takes long enough that
				// the actuator should already be at bias
				data->cmd = data->bias;
				data->stage = PT_DONE;
				pt_analyze(data);
				state->header.status |= AYLP_DONE;
			}
		}
		break;
	}

	case PT_DONE:
		data->cmd = data->bias;
		break;
	}
	data->t_prev = t;
	goto publish;

fatal:
	data->stage = PT_DONE;
	data->cmd = data->bias;
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


int prbs_test_fini(struct aylp_device *self)
{
	struct aylp_prbs_test_data *data = self->device_data;
	xfree(data->output_file);
	xfree(data->config);
	xfree_type(gsl_vector, data->out);
	xfree(data->chips);
	xfree(data->cmd_win);
	xfree(data->resp);
	xfree(data->live);
	xfree(data->rbar);
	xfree(data->rsem);
	xfree(data->rcnt);
	xfree(data->rlive);
	xfree(data->rho);
	xfree(data);
	return 0;
}
