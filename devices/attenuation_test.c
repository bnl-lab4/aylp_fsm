#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/wait.h>
#include <gsl/gsl_fft_real.h>

#include "anyloop.h"
#include "logging.h"
#include "attenuation_test.h"
#include "xalloc.h"


static const char *PLOT_SCRIPT =
	"import sys\n"
	"import matplotlib\n"
	"matplotlib.use('PDF')\n"
	"import matplotlib.pyplot as plt\n"
	"from matplotlib.backends.backend_pdf import PdfPages\n"
	"import numpy as np\n"
	"datpath=sys.argv[1]\n"
	"data=np.loadtxt(datpath)\n"
	"labels=sys.argv[3].split(',') if len(sys.argv)>3 and sys.argv[3] else []\n"
	"# metadata written as # comment lines by the C side\n"
	"meta={}; stats={}; cfg=[]\n"
	"for line in open(datpath):\n"
	"    if not line.startswith('#'): continue\n"
	"    b=line[1:].strip()\n"
	"    if b.startswith('start='): meta['start']=b[6:]\n"
	"    elif b.startswith('pixel_scale='): meta['ps']=float(b.split('=',1)[1])\n"
	"    elif b.startswith('stat '):\n"
	"        p=b.split(); stats[int(p[1])]=[float(x) for x in p[2:6]]\n"
	"    elif b.startswith('config '): cfg.append(b[7:])\n"
	"ps=meta.get('ps',1.0)\n"
	"ulab='px' if ps!=1.0 else 'unit'\n"
	"freq=data[:,0]\n"
	"n_elem=(data.shape[1]-1)//3\n"
	"m=freq>0\n"
	"f=freq[m]\n"
	"pdf=PdfPages(sys.argv[2])\n"
	"fig,axes=plt.subplots(3,n_elem,figsize=(6.5*n_elem,13),sharex='row',squeeze=False)\n"
	"for e in range(n_elem):\n"
	"    po=data[:,1+3*e][m]; pc=data[:,2+3*e][m]; att=data[:,3+3*e][m]\n"
	"    lab=labels[e] if e<len(labels) else 'element %d'%e\n"
	"    ax=axes[0][e]\n"
	"    ax.loglog(f,po,label='open loop',alpha=0.8)\n"
	"    ax.loglog(f,pc,label='closed loop',alpha=0.8)\n"
	"    ax.set_title(lab)\n"
	"    ax.set_ylabel('PSD (unit$^2$/Hz)')\n"
	"    ax.legend()\n"
	"    ax.grid(True,which='both',ls='--',alpha=0.5)\n"
	"    ax=axes[1][e]\n"
	"    ax.semilogx(f,att,lw=0.8)\n"
	"    ax.axhline(0,color='k',lw=0.8)\n"
	"    ax.fill_between(f,att,0,where=att<0,color='tab:green',alpha=0.25,label='attenuated')\n"
	"    ax.fill_between(f,att,0,where=att>0,color='tab:red',alpha=0.25,label='amplified')\n"
	"    # annotate the strongest changes on a lightly smoothed curve,\n"
	"    # keeping picks a few Hz apart so one wide dip isn't labelled thrice\n"
	"    k=np.ones(5)/5\n"
	"    s=np.convolve(att,k,mode='same')\n"
	"    def pick(order,n,thresh):\n"
	"        out=[]\n"
	"        for i in order:\n"
	"            if thresh(s[i]) and all(abs(f[i]-f[j])>5 for j in out):\n"
	"                out.append(i)\n"
	"            if len(out)>=n: break\n"
	"        return out\n"
	"    for i in pick(np.argsort(s),3,lambda v:v<-3):\n"
	"        ax.annotate('%.1f Hz\\n%+.1f dB'%(f[i],att[i]),(f[i],att[i]),\n"
	"                    fontsize=7,ha='center',va='top',color='tab:green')\n"
	"    for i in pick(np.argsort(s)[::-1],2,lambda v:v>3):\n"
	"        ax.annotate('%.1f Hz\\n%+.1f dB'%(f[i],att[i]),(f[i],att[i]),\n"
	"                    fontsize=7,ha='center',va='bottom',color='tab:red')\n"
	"    ax.set_ylabel('closed/open (dB)')\n"
	"    ax.set_xlabel('Frequency (Hz)')\n"
	"    ax.legend(loc='upper left',fontsize=8)\n"
	"    ax.grid(True,which='both',ls='--',alpha=0.5)\n"
	"    # error variance contributed by each 5 Hz band (integral of the\n"
	"    # PSD over the band); top closed-loop bands tagged with their\n"
	"    # share of the total closed-loop error variance\n"
	"    ax=axes[2][e]\n"
	"    df=f[1]-f[0] if len(f)>1 else 1.0\n"
	"    band=(f/5).astype(int)\n"
	"    vo=np.bincount(band,po*df)\n"
	"    vc=np.bincount(band,pc*df)\n"
	"    ctr=(np.arange(len(vo))+0.5)*5\n"
	"    ax.bar(ctr-1.15,vo,width=2.3,label='open loop',alpha=0.8)\n"
	"    ax.bar(ctr+1.15,vc,width=2.3,label='closed loop',alpha=0.8)\n"
	"    ax.set_yscale('log')\n"
	"    for i in np.argsort(vc)[::-1][:3]:\n"
	"        ax.annotate('%.0f%%'%(100*vc[i]/vc.sum()),(ctr[i],vc[i]),\n"
	"                    fontsize=7,ha='center',va='bottom')\n"
	"    ax.set_ylabel('error variance per 5 Hz band (unit$^2$)')\n"
	"    ax.set_xlabel('Frequency (Hz)')\n"
	"    ax.legend(loc='upper right',fontsize=8)\n"
	"    ax.grid(True,which='both',ls='--',alpha=0.5)\n"
	"fig.suptitle('Closed-loop attenuation')\n"
	"plt.tight_layout()\n"
	"pdf.savefig(fig,bbox_inches='tight')\n"
	"plt.close(fig)\n"
	"# second page: run config, mean/jitter summary, and the RMS that each\n"
	"# 10 Hz band contributes (open vs closed), all scaled to pixels\n"
	"df=f[1]-f[0] if len(f)>1 else 1.0\n"
	"edges=np.arange(0,310,10)\n"
	"def bandrms(psd):\n"
	"    return [ps*np.sqrt(np.sum(psd[(f>=lo)&(f<lo+10)])*df) for lo in edges[:-1]]\n"
	"bo={e:bandrms(data[:,1+3*e][m]) for e in range(n_elem)}\n"
	"bc={e:bandrms(data[:,2+3*e][m]) for e in range(n_elem)}\n"
	"fig2=plt.figure(figsize=(6.5*max(n_elem,1),13))\n"
	"fig2.suptitle('Run configuration and per-band RMS')\n"
	"info=[]\n"
	"if 'start' in meta: info.append('Started: '+meta['start'])\n"
	"info.append('pixel_scale: %g %s per output unit'%(ps,ulab))\n"
	"for e in range(n_elem):\n"
	"    lab=labels[e] if e<len(labels) else 'elem %d'%e\n"
	"    if e in stats:\n"
	"        mo,so,mc,sc=stats[e]\n"
	"        info.append('%s: mean %+.3f -> %+.3f %s;  RMS about mean %.3f -> %.3f %s   (open -> closed)'%(lab,ps*mo,ps*mc,ulab,ps*so,ps*sc,ulab))\n"
	"txt='\\n'.join(info)\n"
	"if cfg: txt+='\\n\\nConfig used:\\n'+'\\n'.join(cfg)\n"
	"fig2.text(0.03,0.97,txt,va='top',ha='left',family='monospace',fontsize=8)\n"
	"col=['band (Hz)']\n"
	"for e in range(n_elem):\n"
	"    lab=labels[e] if e<len(labels) else 'elem %d'%e\n"
	"    col+=['%s open'%lab,'%s closed'%lab]\n"
	"rows=[]\n"
	"for i,lo in enumerate(edges[:-1]):\n"
	"    r=['%d-%d'%(lo,lo+10)]\n"
	"    for e in range(n_elem): r+=['%.4f'%bo[e][i],'%.4f'%bc[e][i]]\n"
	"    rows.append(r)\n"
	"tot=['total 0-300']\n"
	"for e in range(n_elem):\n"
	"    tot+=['%.4f'%np.sqrt(np.sum(np.array(bo[e])**2)),'%.4f'%np.sqrt(np.sum(np.array(bc[e])**2))]\n"
	"rows.append(tot)\n"
	"axt=fig2.add_axes([0.03,0.02,0.94,0.70]); axt.axis('off')\n"
	"axt.set_title('RMS contributed by each 10 Hz band (%s)'%ulab,fontsize=9)\n"
	"tb=axt.table(cellText=rows,colLabels=col,loc='upper center',cellLoc='center')\n"
	"tb.auto_set_font_size(False); tb.set_fontsize(7); tb.scale(1,1.15)\n"
	"for j in range(len(col)): tb[0,j].set_facecolor('#dddddd')\n"
	"for j in range(len(col)): tb[len(rows),j].set_facecolor('#eeeeee')\n"
	"pdf.savefig(fig2,bbox_inches='tight')\n"
	"plt.close(fig2)\n"
	"pdf.close()\n"
	"print('Saved attenuation plot to '+sys.argv[2])\n";


static inline double get_time_seconds(void)
{
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return ts.tv_sec + ts.tv_nsec * 1e-9;
}


// append one sample vector to an interleaved buffer, growing it as needed
static void record_sample(const gsl_vector *v, size_t n_elem,
	double **buf, size_t *n, size_t *cap)
{
	if (*n >= *cap) {
		*cap = *cap ? *cap * 2 : 4096;
		*buf = xrealloc(*buf, *cap * n_elem * sizeof(double));
	}
	for (size_t e = 0; e < n_elem; e++)
		(*buf)[*n * n_elem + e] = gsl_vector_get(v, e);
	*n += 1;
}


/** One-sided Welch PSD with Hann window, mean removal, and 50% overlap.
 * x has n samples; segments are nfft samples long (any length <= n).
 * nfft_fft (a power of 2, >= nfft) zero-pads each segment before the FFT,
 * sampling the same windowed spectrum on a grid nfft_fft/nfft times finer;
 * psd must hold nfft_fft/2+1 bins, spaced fs/nfft_fft apart. */
static void welch_psd(const double *x, size_t n, size_t nfft,
	size_t nfft_fft, double fs, double *psd)
{
	size_t nbins = nfft_fft/2 + 1;
	memset(psd, 0, nbins * sizeof(double));
	double *w = xmalloc(nfft * sizeof(double));
	double u = 0.0;	// window power, sum(w^2)
	for (size_t i = 0; i < nfft; i++) {
		w[i] = 0.5 * (1.0 - cos(2.0*M_PI*i / (nfft-1)));
		u += w[i]*w[i];
	}
	size_t hop = nfft/2;
	size_t nseg = 1 + (n - nfft)/hop;
	double *seg = xmalloc(nfft_fft * sizeof(double));
	for (size_t s = 0; s < nseg; s++) {
		const double *xs = x + s*hop;
		double mean = 0.0;
		for (size_t i = 0; i < nfft; i++) mean += xs[i];
		mean /= nfft;
		for (size_t i = 0; i < nfft; i++)
			seg[i] = (xs[i] - mean) * w[i];
		memset(seg + nfft, 0, (nfft_fft - nfft) * sizeof(double));
		gsl_fft_real_radix2_transform(seg, 1, nfft_fft);
		// halfcomplex layout: re(k)=seg[k], im(k)=seg[nfft_fft-k]
		psd[0] += seg[0]*seg[0];
		for (size_t k = 1; k < nfft_fft/2; k++)
			psd[k] += seg[k]*seg[k]
				+ seg[nfft_fft-k]*seg[nfft_fft-k];
		psd[nfft_fft/2] += seg[nfft_fft/2]*seg[nfft_fft/2];
	}
	double scale = 1.0 / (fs * u * nseg);
	psd[0] *= scale;
	psd[nfft_fft/2] *= scale;
	for (size_t k = 1; k < nfft_fft/2; k++)
		psd[k] *= 2.0 * scale;
	xfree(w);
	xfree(seg);
}


// extract one element's time series from an interleaved buffer
static void extract_element(const double *buf, size_t n, size_t n_elem,
	size_t e, double *out)
{
	for (size_t i = 0; i < n; i++)
		out[i] = buf[i*n_elem + e];
}


// mean and RMS-about-mean (population std) of one element of a phase buffer,
// in raw output units; the report scales both by pixel_scale
static void phase_stats(const double *buf, size_t n, size_t n_elem, size_t e,
	double *mean, double *rms)
{
	double s = 0.0, s2 = 0.0;
	for (size_t i = 0; i < n; i++) {
		double v = buf[i*n_elem + e];
		s += v;
		s2 += v*v;
	}
	double m = n ? s/n : 0.0;
	double var = n ? s2/n - m*m : 0.0;
	*mean = m;
	*rms = var > 0.0 ? sqrt(var) : 0.0;
}


static int analyze_and_plot(struct aylp_attenuation_test_data *data)
{
	if (data->open_n < 2 || data->closed_n < 2) {
		log_error("attenuation_test: not enough samples "
			"(open %zu, closed %zu)", data->open_n, data->closed_n);
		return -1;
	}
	double fs_open = (data->open_n - 1)
		/ (data->open_t_last - data->open_t_first);
	double fs_closed = (data->closed_n - 1)
		/ (data->closed_t_last - data->closed_t_first);
	// closing the loop usually slows the iteration (DAC writes), so the
	// phases legitimately run at different rates; each PSD is computed on
	// its own grid and the open one is re-gridded onto the closed axis
	// below, so differing rates are handled, not just tolerated
	if (fs_open != fs_closed)
		log_info("attenuation_test: sample rate %g Hz open, "
			"%g Hz closed", fs_open, fs_closed);
	if (fs_closed > fs_open)
		log_warn("attenuation_test: closed phase ran faster than "
			"open; the attenuation curve above the open Nyquist "
			"(%g Hz) compares against a clamped open PSD",
			fs_open / 2.0);
	double fs = fs_closed;

	// clamp nfft to a power of two that fits the shorter recording
	size_t n_min = data->open_n < data->closed_n
		? data->open_n : data->closed_n;
	size_t nfft = data->nfft;
	while (nfft > n_min) nfft /= 2;
	if (nfft < 64) {
		log_error("attenuation_test: recordings too short for a PSD "
			"(min %zu samples)", n_min);
		return -1;
	}
	if (nfft != data->nfft)
		log_warn("attenuation_test: nfft reduced %zu -> %zu to fit "
			"the recording", data->nfft, nfft);

	size_t nbins = nfft/2 + 1;
	size_t n_elem = data->n_elem;
	double *psd_open = xmalloc(n_elem * nbins * sizeof(double));
	double *psd_closed = xmalloc(n_elem * nbins * sizeof(double));
	size_t n_max = data->open_n > data->closed_n
		? data->open_n : data->closed_n;
	double *series = xmalloc(n_max * sizeof(double));
	// the output frequency axis is the closed phase's (k*fs_closed/nfft),
	// but the open PSD lives on its own grid. Two things make the
	// regridded dB ratio exact when the rates differ: (1) the open
	// segment spans the SAME TIME as a closed segment (nfft*fs_open/
	// fs_closed samples -- only the FFT length must be a power of 2, not
	// the window), so both windows have the same resolution bandwidth in
	// Hz and an off-bin line scallops identically in both PSDs, cancelling
	// in the ratio; (2) the open FFTs are zero-padded to >= 8x the segment
	// length, which samples the exact windowed spectrum finely enough that
	// the linear interpolation onto the closed axis errs < 0.01 dB
	// (whole-bin interpolation flattened line peaks by a few tenths of a
	// dB). Bins past the open Nyquist clamp to the last open bin.
	size_t seg_open = nfft;
	size_t nfft_open = nfft;
	if (fs_open != fs_closed) {
		seg_open = (size_t)lround(nfft * fs_open / fs_closed);
		if (seg_open > data->open_n)
			seg_open = data->open_n;
		nfft_open = 1;
		while (nfft_open < 8*seg_open) nfft_open *= 2;
	}
	size_t nbins_fine = nfft_open/2 + 1;
	double *fine = fs_open != fs_closed
		? xmalloc(nbins_fine * sizeof(double)) : NULL;
	for (size_t e = 0; e < n_elem; e++) {
		extract_element(data->open_buf, data->open_n, n_elem, e,
			series);
		if (fine) {
			welch_psd(series, data->open_n, seg_open,
				nfft_open, fs_open, fine);
			double *po = psd_open + e*nbins;
			for (size_t k = 0; k < nbins; k++) {
				// closed-axis freq in open fine-grid bins
				// (spacing fs_open/nfft_open)
				double pos = (double)k * fs_closed * nfft_open
					/ (fs_open * nfft);
				size_t i0 = (size_t)pos;
				if (i0 >= nbins_fine - 1) {
					po[k] = fine[nbins_fine - 1];
				} else {
					double frac = pos - i0;
					po[k] = fine[i0] * (1.0 - frac)
						+ fine[i0 + 1] * frac;
				}
			}
		} else {
			welch_psd(series, data->open_n, nfft, nfft, fs_open,
				psd_open + e*nbins);
		}
		extract_element(data->closed_buf, data->closed_n, n_elem, e,
			series);
		welch_psd(series, data->closed_n, nfft, nfft, fs_closed,
			psd_closed + e*nbins);
	}
	xfree(fine);
	xfree(series);

	// data file next to the PDF, like bode_plot does
	const char *out = data->output_file;
	const char *dot = strrchr(out, '.');
	size_t base = dot ? (size_t)(dot - out) : strlen(out);
	char *dat_path = xmalloc(base + 5);
	memcpy(dat_path, out, base);
	memcpy(dat_path + base, ".dat", 5);
	FILE *f = fopen(dat_path, "w");
	if (!f) {
		log_error("attenuation_test: fopen %s: %m", dat_path);
		xfree(dat_path);
		xfree(psd_open);
		xfree(psd_closed);
		return -1;
	}
	fprintf(f, "# fs_open=%g fs_closed=%g nfft=%zu open_s=%g closed_s=%g"
		" labels=%s\n", fs_open, fs_closed, nfft,
		data->open_t_last - data->open_t_first,
		data->closed_t_last - data->closed_t_first,
		data->labels ? data->labels : "");
	// report metadata for the PDF: start time, unit->px scale, per-element
	// mean and RMS-about-mean (raw units; the plot scales by pixel_scale),
	// and the caller-supplied config summary (one line each)
	char start_str[32] = "";
	struct tm tmv;
	if (localtime_r(&data->start_wall, &tmv))
		strftime(start_str, sizeof start_str, "%Y-%m-%d %H:%M:%S", &tmv);
	fprintf(f, "# start=%s\n", start_str);
	fprintf(f, "# pixel_scale=%g\n", data->pixel_scale);
	for (size_t e = 0; e < n_elem; e++) {
		double mo, so, mc, sc;
		phase_stats(data->open_buf, data->open_n, n_elem, e, &mo, &so);
		phase_stats(data->closed_buf, data->closed_n, n_elem, e,
			&mc, &sc);
		fprintf(f, "# stat %zu %g %g %g %g\n", e, mo, so, mc, sc);
	}
	if (data->config) {
		const char *p = data->config;
		while (*p) {
			const char *nl = strchr(p, '\n');
			int len = nl ? (int)(nl - p) : (int)strlen(p);
			fprintf(f, "# config %.*s\n", len, p);
			if (!nl) break;
			p = nl + 1;
		}
	}
	fprintf(f, "# freq_Hz");
	for (size_t e = 0; e < n_elem; e++)
		fprintf(f, " psd_open_%zu psd_closed_%zu atten_dB_%zu",
			e, e, e);
	fprintf(f, "\n");
	for (size_t k = 0; k < nbins; k++) {
		fprintf(f, "%g", k * fs / nfft);
		for (size_t e = 0; e < n_elem; e++) {
			double po = psd_open[e*nbins + k];
			double pc = psd_closed[e*nbins + k];
			fprintf(f, " %g %g %g", po, pc,
				10.0 * log10((pc + 1e-300)/(po + 1e-300)));
		}
		fprintf(f, "\n");
	}
	fclose(f);
	log_info("Attenuation data saved to %s", dat_path);
	xfree(psd_open);
	xfree(psd_closed);

	char script_path[64];
	snprintf(script_path, sizeof(script_path),
		"/tmp/attenuation_plot_%d.py", getpid());
	FILE *script = fopen(script_path, "w");
	if (!script) {
		log_error("attenuation_test: fopen %s: %m", script_path);
		xfree(dat_path);
		return -1;
	}
	fputs(PLOT_SCRIPT, script);
	fclose(script);

	char cmd[512];
	snprintf(cmd, sizeof(cmd), "python3 '%s' '%s' '%s' '%s'", script_path,
		dat_path, data->output_file,
		data->labels ? data->labels : "");
	int ret = system(cmd);
	int ok = (ret != -1) && WIFEXITED(ret) && (WEXITSTATUS(ret) == 0);
	if (!ok)
		log_error("Attenuation plot script failed (exit status %d); "
			"is python3 with numpy and matplotlib installed?",
			(ret != -1 && WIFEXITED(ret)) ? WEXITSTATUS(ret) : ret);
	else
		log_info("Attenuation plot saved to %s", data->output_file);
	unlink(script_path);
	xfree(dat_path);
	return ok ? 0 : -1;
}


int attenuation_test_init(struct aylp_device *self)
{
	self->proc = &attenuation_test_proc;
	self->fini = &attenuation_test_fini;
	self->device_data =
		xcalloc(1, sizeof(struct aylp_attenuation_test_data));
	struct aylp_attenuation_test_data *data = self->device_data;

	// defaults
	data->start_delay = 0.0;
	data->settle_time = 5.0;
	data->nfft        = 4096;
	data->output_file = xstrdup("attenuation.pdf");
	data->pixel_scale = 1.0;

	if (!self->params) {
		log_error("No params object found.");
		return -1;
	}
	json_object_object_foreach(self->params, key, val) {
		if (key[0] == '_') {
			// comment
		} else if (!strcmp(key, "start_delay")) {
			data->start_delay = json_object_get_double(val);
			log_trace("start_delay = %G s", data->start_delay);
		} else if (!strcmp(key, "open_time")) {
			data->open_time = json_object_get_double(val);
			log_trace("open_time = %G s", data->open_time);
		} else if (!strcmp(key, "settle_time")) {
			data->settle_time = json_object_get_double(val);
			log_trace("settle_time = %G s", data->settle_time);
		} else if (!strcmp(key, "closed_time")) {
			data->closed_time = json_object_get_double(val);
			log_trace("closed_time = %G s", data->closed_time);
		} else if (!strcmp(key, "nfft")) {
			data->nfft = json_object_get_uint64(val);
			log_trace("nfft = %zu", data->nfft);
		} else if (!strcmp(key, "output_file")) {
			xfree(data->output_file);
			data->output_file =
				xstrdup(json_object_get_string(val));
			log_trace("output_file = %s", data->output_file);
		} else if (!strcmp(key, "labels")) {
			xfree(data->labels);
			data->labels = xstrdup(json_object_get_string(val));
			log_trace("labels = %s", data->labels);
		} else if (!strcmp(key, "pixel_scale")) {
			data->pixel_scale = json_object_get_double(val);
			log_trace("pixel_scale = %G", data->pixel_scale);
		} else if (!strcmp(key, "config")) {
			xfree(data->config);
			data->config = xstrdup(json_object_get_string(val));
			log_trace("config = %s", data->config);
		} else if (!strcmp(key, "pass_open")) {
			data->pass_open = json_object_get_boolean(val);
			log_trace("pass_open = %d", data->pass_open);
		} else {
			log_warn("Unknown parameter \"%s\"", key);
		}
	}

	if (data->open_time <= 0.0 || data->closed_time <= 0.0) {
		log_error("open_time and closed_time are required "
			"and must be > 0.");
		return -1;
	}
	if (data->nfft < 64 || (data->nfft & (data->nfft - 1))) {
		log_error("nfft must be a power of 2 and >= 64.");
		return -1;
	}
	if (data->pixel_scale <= 0.0) {
		log_error("pixel_scale must be > 0.");
		return -1;
	}
	if (data->labels && strchr(data->labels, '\'')) {
		// the label string is single-quoted into a shell command
		log_error("labels must not contain '");
		return -1;
	}
	if (strchr(data->output_file, '\'')) {
		// The output and adjacent .dat path are single-quoted when the plot
		// helper is launched, allowing destination directories with spaces.
		log_error("output_file must not contain '");
		return -1;
	}

	log_info("attenuation_test: %g s start delay, %g s open, "
		"%g s settle, %g s closed", data->start_delay,
		data->open_time, data->settle_time, data->closed_time);

	self->type_in = AYLP_T_VECTOR;
	self->units_in = AYLP_U_ANY;
	self->type_out = AYLP_T_UNCHANGED;
	self->units_out = AYLP_U_UNCHANGED;
	return 0;
}


int attenuation_test_proc(struct aylp_device *self, struct aylp_state *state)
{
	struct aylp_attenuation_test_data *data = self->device_data;
	double now = get_time_seconds();

	if (!data->t0) {
		data->t0 = now;
		data->start_wall = time(NULL);
		data->n_elem = state->vector->size;
	}
	double elapsed = now - data->t0;

	// phase boundaries, from loop start
	double t_open   = data->start_delay;
	double t_settle = t_open + data->open_time;
	double t_closed = t_settle + data->settle_time;
	double t_done   = t_closed + data->closed_time;

	unsigned phase;
	if      (elapsed < t_open)   phase = ATTEN_PHASE_START;
	else if (elapsed < t_settle) phase = ATTEN_PHASE_OPEN;
	else if (elapsed < t_closed) phase = ATTEN_PHASE_SETTLE;
	else if (elapsed < t_done)   phase = ATTEN_PHASE_CLOSED;
	else                         phase = ATTEN_PHASE_DONE;

	if (phase != data->phase) {
		switch (phase) {
		case ATTEN_PHASE_OPEN:
			log_info("attenuation_test: recording open loop "
				"for %g s", data->open_time);
			break;
		case ATTEN_PHASE_SETTLE:
			log_info("attenuation_test: loop closed; settling "
				"for %g s", data->settle_time);
			break;
		case ATTEN_PHASE_CLOSED:
			log_info("attenuation_test: recording closed loop "
				"for %g s", data->closed_time);
			break;
		}
		data->phase = phase;
	}

	switch (phase) {
	case ATTEN_PHASE_START:
		gsl_vector_set_zero(state->vector);
		return 0;
	case ATTEN_PHASE_OPEN:
		record_sample(state->vector, data->n_elem,
			&data->open_buf, &data->open_n, &data->open_cap);
		if (data->open_n == 1) data->open_t_first = now;
		data->open_t_last = now;
		// Normally hold the downstream error at zero. A predictive
		// controller may instead request the raw measurement so it can
		// identify while its own matched startup hold keeps command at zero.
		if (!data->pass_open) gsl_vector_set_zero(state->vector);
		return 0;
	case ATTEN_PHASE_SETTLE:
		return 0;
	case ATTEN_PHASE_CLOSED:
		record_sample(state->vector, data->n_elem,
			&data->closed_buf, &data->closed_n, &data->closed_cap);
		if (data->closed_n == 1) data->closed_t_first = now;
		data->closed_t_last = now;
		return 0;
	default: {	// ATTEN_PHASE_DONE
		int ret = 0;
		if (!data->analyzed) {
			data->analyzed = 1;
			ret = analyze_and_plot(data);
		}
		state->header.status |= AYLP_DONE;
		return ret;
	}
	}
}


int attenuation_test_fini(struct aylp_device *self)
{
	struct aylp_attenuation_test_data *data = self->device_data;
	xfree(data->open_buf);
	xfree(data->closed_buf);
	xfree(data->output_file);
	xfree(data->labels);
	xfree(data->config);
	xfree(data);
	return 0;
}
