#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <termios.h>
#include <unistd.h>
#include <sys/wait.h>
#include <time.h>

#include "anyloop.h"
#include "logging.h"
#include "bode_plot.h"
#include "xalloc.h"

static int fg_serial_open(const char *port)
{
	int fd = open(port, O_RDWR | O_NOCTTY);
	if (fd < 0) return -1;
	struct termios tty;
	if (tcgetattr(fd, &tty)) { close(fd); return -1; }
	cfmakeraw(&tty);
	cfsetispeed(&tty, B115200);
	cfsetospeed(&tty, B115200);
	tty.c_cc[VMIN]  = 0;
	tty.c_cc[VTIME] = 5;
	if (tcsetattr(fd, TCSANOW, &tty)) { close(fd); return -1; }
	return fd;
}

static void fg_drain(int fd)
{
	char c;
	while (read(fd, &c, 1) == 1 && c != '\n')
		;
}

static int fg_send(int fd, const char *cmd)
{
	size_t n = strlen(cmd);
	if (write(fd, cmd, n) != (ssize_t)n) return -1;
	fg_drain(fd);
	return 0;
}


static const char *PLOT_SCRIPT =
	"import sys\n"
	"import matplotlib\n"
	"matplotlib.use('PDF')\n"
	"import matplotlib.pyplot as plt\n"
	"import numpy as np\n"
	"data=np.loadtxt(sys.argv[1])\n"
	"f,m,p=data[:,0],data[:,1],data[:,2]\n"
	"ms,ps,q=data[:,3],data[:,4],data[:,5]\n"
	"ok=m>1e-6\n"
	"f,m,p,ms,ps,q=f[ok],m[ok],p[ok],ms[ok],ps[ok],q[ok]\n"
	"pd=np.degrees(np.unwrap(p))\n"
	"fig,(ax1,ax2)=plt.subplots(2,1,figsize=(10,8),sharex=True)\n"
	"fig.suptitle(sys.argv[3],fontsize=9)\n"
	"mdb=20*np.log10(m+1e-300)\n"
	"mdb_err=8.6859*ms/np.maximum(m,1e-12)\n"
	"ax1.errorbar(f,mdb,yerr=mdb_err,fmt='.-',lw=1,ms=3,capsize=2)\n"
	"bad=q<0.9\n"
	"ax1.plot(f[bad],mdb[bad],'rx',ms=6,label='quality<0.9')\n"
	"if bad.any(): ax1.legend(fontsize=8)\n"
	"ax1.set_xscale('log')\n"
	"ax1.set_ylabel('Magnitude (dB)')\n"
	"ax1.grid(True,which='both',ls='--',alpha=0.5)\n"
	"ax2.errorbar(f,pd,yerr=np.degrees(ps),fmt='.-',lw=1,ms=3,capsize=2)\n"
	"ax2.plot(f[bad],pd[bad],'rx',ms=6)\n"
	"if data.shape[1]>=8:\n"
	"    m2,p2=data[:,6][ok],data[:,7][ok]\n"
	"    ax1.plot(f,20*np.log10(m2+1e-300),'--',color='tab:orange',lw=1,label='amplitude2')\n"
	"    ax2.plot(f,np.degrees(np.unwrap(p2)),'--',color='tab:orange',lw=1)\n"
	"    ax1.legend(fontsize=8)\n"
	"ax2.set_xscale('log')\n"
	"ax2.set_ylabel('Phase (deg)')\n"
	"ax2.set_xlabel('Frequency (Hz)')\n"
	"ax2.grid(True,which='both',ls='--',alpha=0.5)\n"
	"plt.tight_layout()\n"
	"plt.savefig(sys.argv[2],format='pdf',bbox_inches='tight')\n"
	"print('Saved Bode plot to '+sys.argv[2])\n";


static inline double get_time_seconds(void)
{
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return ts.tv_sec + ts.tv_nsec * 1e-9;
}

// flat-gain + pure-delay fit over the pass-0 points: K = error-weighted mean
// |H| below k_max_freq; tau from the error-weighted LS slope of the unwrapped
// phase. Points with quality < 0.5 are excluded. The measured phase EXCLUDES
// one frame (the correlator references phi[n-1]); fit_frames adds it back
// using the MEASURED mean sample rate.
static void bode_fit(struct aylp_bode_plot_data *data)
{
	size_t n = data->n_freqs;
	double ph_u_prev = 0.0;
	double sw = 0, swf = 0, swf2 = 0, swp = 0, swfp = 0;
	double kw = 0, kwm = 0;
	size_t n_used = 0, n_k = 0;
	double fs_meas = data->dt_n ? (double)data->dt_n / data->dt_sum
		: data->sample_rate;

	// first pass: weighted phase line fit + weighted K
	for (size_t i = 0; i < n; i++) {
		if (data->mags[i] <= 1e-9 || data->quality[i] < 0.5) continue;
		double ph = data->phases[i];
		if (n_used)	// unwrap toward the previous used point
			ph += 2.0*M_PI * round((ph_u_prev - ph) / (2.0*M_PI));
		ph_u_prev = ph;
		double sem = data->ph_sem[i] > 1e-3 ? data->ph_sem[i] : 1e-3;
		double w = 1.0 / (sem * sem);
		double f = data->freqs[i];
		sw += w; swf += w*f; swf2 += w*f*f; swp += w*ph; swfp += w*f*ph;
		n_used++;
		if (f <= data->k_max_freq) {
			double msem = data->mag_sem[i] > 1e-6
				? data->mag_sem[i] : 1e-6;
			double mw = 1.0 / (msem * msem);
			kw += mw;
			kwm += mw * data->mags[i];
			n_k++;
		}
	}
	if (n_used < 3 || kw <= 0.0) {
		log_warn("bode_plot: too few clean points for the K/tau fit "
			"(%zu used)", n_used);
		return;
	}
	double det = sw*swf2 - swf*swf;
	if (fabs(det) < 1e-30) return;
	double b = (sw*swfp - swf*swp) / det;	// rad per Hz
	double a = (swf2*swp - swf*swfp) / det;
	data->fit_tau_ms = -b / (2.0*M_PI) * 1e3;
	data->fit_frames = data->fit_tau_ms * 1e-3 * fs_meas + 1.0;
	data->fit_K = kwm / kw;
	data->fit_K_err = sqrt(1.0 / kw) * sqrt((double)(n_k > 1 ? n_k : 1));

	// residuals (unweighted rms over the used points)
	double pr2 = 0, mr2 = 0;
	size_t n2 = 0;
	ph_u_prev = 0.0;
	bool first = true;
	for (size_t i = 0; i < n; i++) {
		if (data->mags[i] <= 1e-9 || data->quality[i] < 0.5) continue;
		double ph = data->phases[i];
		if (!first)
			ph += 2.0*M_PI * round((ph_u_prev - ph) / (2.0*M_PI));
		first = false;
		ph_u_prev = ph;
		double r = ph - (a + b*data->freqs[i]);
		pr2 += r*r;
		double mr = (data->mags[i] - data->fit_K) / data->fit_K;
		mr2 += mr*mr;
		n2++;
	}
	data->fit_ph_resid_deg = sqrt(pr2 / n2) * 180.0 / M_PI;
	data->fit_mag_flat_rms = sqrt(mr2 / n2);

	log_info("bode_plot FIT: K(<=%G Hz) = %.4f +/- %.4f, tau = %.4f ms "
		"(+1 hidden frame => %.3f frames at measured fs %.1f Hz); "
		"phase resid %.2f deg rms, |H| flatness %.1f%% rms over %zu pts",
		data->k_max_freq, data->fit_K, data->fit_K_err,
		data->fit_tau_ms, data->fit_frames, fs_meas,
		data->fit_ph_resid_deg, 100.0*data->fit_mag_flat_rms, n_used);
	if (data->amplitude2 > 0.0) {
		// Robust linearity: a real amplitude nonlinearity moves MANY
		// points together, so judge on the MEDIAN |ratio-1| and the
		// count over threshold -- a lone contaminated segment (one bad
		// amp2 point) must not condemn K. Only points clean in the main
		// pass (quality >= 0.9) count.
		double dev[512];
		size_t nd = 0, n_over = 0;
		double worst = 0.0, wf = 0.0;
		for (size_t i = 0; i < n && nd < 512; i++) {
			if (data->mags[i] <= 1e-9 || data->mags2[i] <= 1e-9
					|| data->quality[i] < 0.9)
				continue;
			double d = fabs(data->mags2[i]/data->mags[i] - 1.0);
			dev[nd++] = d;
			if (d > 0.10) n_over++;
			if (d > worst) { worst = d; wf = data->freqs[i]; }
		}
		double med = 0.0;
		if (nd) {
			for (size_t a = 0; a < nd; a++)	// tiny insertion sort
				for (size_t b = a+1; b < nd; b++)
					if (dev[b] < dev[a]) {
						double t = dev[a];
						dev[a] = dev[b]; dev[b] = t;
					}
			med = dev[nd/2];
		}
		bool nonlin = med > 0.05 || n_over > nd/10 + 1;
		log_info("bode_plot LINEARITY: median |H| change between "
			"amplitudes %G and %G is %.1f%% (%zu/%zu pts over 10%%, "
			"worst %.1f%% at %.1f Hz) -- %s", data->amplitude,
			data->amplitude2, 100.0*med, n_over, nd, 100.0*worst, wf,
			nonlin ? "NONLINEAR/CLIPPING, do not trust K"
			: "plant is linear, K trustworthy");
	}
}

// write the annotated .dat (header + per-point columns) to path
static int write_dat(struct aylp_bode_plot_data *data, const char *path)
{
	FILE *f = fopen(path, "w");
	if (!f) {
		log_error("bode_plot: fopen %s: %m", path);
		return -1;
	}
	time_t tt = time(NULL);
	struct tm tm;
	localtime_r(&tt, &tm);
	double fs_meas = data->dt_n ? (double)data->dt_n / data->dt_sum
		: data->sample_rate;
	fprintf(f, "# bode_plot v2 %04d-%02d-%02d %02d:%02d:%02d\n",
		tm.tm_year+1900, tm.tm_mon+1, tm.tm_mday,
		tm.tm_hour, tm.tm_min, tm.tm_sec);
	fprintf(f, "# element=%zu inject=%s output=%s amplitude=%G amplitude2=%G\n",
		data->element,
		data->inject_mode == BODE_INJECT_ADD ? "add" : "replace",
		data->output_mode == BODE_OUTPUT_RAW ? "raw" : "diff",
		data->amplitude, data->amplitude2);
	fprintf(f, "# n_freqs=%zu cycles=%G (%zu segments x %G cycles) "
		"n_settle=%zu+%Gcyc\n", data->n_freqs, data->n_cycles,
		data->n_segments, data->cycles_per_seg,
		data->n_settle, data->n_settle_cycles);
	fprintf(f, "# fs_nominal=%G fs_measured=%.2f\n",
		data->sample_rate, fs_meas);
	if (data->fit_K != 0.0) {
		fprintf(f, "# fit K(<=%GHz)=%.4f +/- %.4f  tau_meas_ms=%.4f  "
			"loop_frames=%.3f (tau*fs_meas + 1 hidden frame)  "
			"ph_resid_deg=%.2f  mag_flat_rms=%.3f\n",
			data->k_max_freq, data->fit_K, data->fit_K_err,
			data->fit_tau_ms, data->fit_frames,
			data->fit_ph_resid_deg, data->fit_mag_flat_rms);
	}
	fprintf(f, "# note: phase EXCLUDES one frame (correlator references "
		"phi[n-1]); loop_frames above already adds it back\n");
	if (data->config_note)
		fprintf(f, "# config %s\n", data->config_note);
	fprintf(f, "# freq_Hz mag phase_rad mag_sem phase_sem_rad quality%s\n",
		data->amplitude2 > 0.0 ? " mag2 phase2_rad" : "");
	for (size_t k = 0; k < data->n_freqs; k++) {
		fprintf(f, "%g %g %g %g %g %g", data->freqs[k], data->mags[k],
			data->phases[k], data->mag_sem[k], data->ph_sem[k],
			data->quality[k]);
		if (data->amplitude2 > 0.0)
			fprintf(f, " %g %g", data->mags2[k], data->phases2[k]);
		fputc('\n', f);
	}
	fclose(f);
	log_info("Bode data saved to %s", path);
	return 0;
}

static int generate_plot(struct aylp_bode_plot_data *data)
{
	// .dat path: output_file with the extension swapped for .dat
	const char *out = data->output_file;
	const char *dot = strrchr(out, '.');
	size_t base = dot ? (size_t)(dot - out) : strlen(out);
	char *dat_path = xmalloc(base + 5);
	memcpy(dat_path, out, base);
	memcpy(dat_path + base, ".dat", 5);
	if (write_dat(data, dat_path)) {
		xfree(dat_path);
		return -1;
	}

	char script_path[64];
	snprintf(script_path, sizeof(script_path), "/tmp/bode_plot_%d.py",
		getpid());
	FILE *script = fopen(script_path, "w");
	if (!script) {
		log_error("fopen for plot script failed: %m");
		xfree(dat_path);
		return -1;
	}
	fputs(PLOT_SCRIPT, script);
	fclose(script);

	char title[256];
	if (data->fit_K != 0.0)
		snprintf(title, sizeof title,
			"K(<=%gHz)=%.4f+/-%.4f  tau=%.4fms => %.3f frames "
			"(incl 1 hidden)  ph resid %.2fdeg  flat %.1f%%",
			data->k_max_freq, data->fit_K, data->fit_K_err,
			data->fit_tau_ms, data->fit_frames,
			data->fit_ph_resid_deg, 100.0*data->fit_mag_flat_rms);
	else
		snprintf(title, sizeof title, "Bode Plot");

	char cmd[1024];
	snprintf(cmd, sizeof(cmd), "python3 %s %s %s '%s'", script_path,
		dat_path, data->output_file, title);
	int ret = system(cmd);
	int ok = (ret != -1) && WIFEXITED(ret) && (WEXITSTATUS(ret) == 0);
	if (!ok)
		log_error("Bode plot script failed (exit status %d); "
			"is python3 with numpy and matplotlib installed?",
			(ret != -1 && WIFEXITED(ret)) ? WEXITSTATUS(ret) : ret);
	else
		log_info("Bode plot saved to %s", data->output_file);
	unlink(script_path);
	xfree(dat_path);
	return ok ? 0 : -1;
}

// get the scalar of interest from the current pipeline state
static double get_element(struct aylp_state *state, size_t elem)
{
	switch (state->header.type) {
	case AYLP_T_VECTOR:
		return gsl_vector_get(state->vector, elem);
	case AYLP_T_MATRIX:
		return gsl_matrix_get(state->matrix,
			elem / state->matrix->size2,
			elem % state->matrix->size2);
	default:
		return 0.0;
	}
}

// write to element in the current pipeline state (add or replace)
static void write_element(struct aylp_state *state, size_t elem,
	double val, unsigned inject_mode)
{
	double out = val;
	switch (state->header.type) {
	case AYLP_T_VECTOR:
		if (inject_mode == BODE_INJECT_ADD)
			out += gsl_vector_get(state->vector, elem);
		gsl_vector_set(state->vector, elem, out);
		break;
	case AYLP_T_MATRIX: {
		size_t row = elem / state->matrix->size2;
		size_t col = elem % state->matrix->size2;
		if (inject_mode == BODE_INJECT_ADD)
			out += gsl_matrix_get(state->matrix, row, col);
		gsl_matrix_set(state->matrix, row, col, out);
		break;
	}
	default:
		break;
	}
}

// reset per-frequency measurement state
static void advance_frequency(struct aylp_bode_plot_data *data)
{
	double now = get_time_seconds();
	data->start_time = now;
	data->prev_time  = now;
	data->integrating = false;
	data->seg_idx = 0;
	data->seg_Yc = data->seg_Ys = data->seg_T = 0.0;
	data->seg_start_phi = 0.0;
}


int bode_plot_init(struct aylp_device *self)
{
	self->proc = &bode_plot_proc;
	self->fini = &bode_plot_fini;
	self->device_data = xcalloc(1, sizeof(struct aylp_bode_plot_data));
	struct aylp_bode_plot_data *data = self->device_data;

	// defaults
	data->inject_mode  = BODE_INJECT_ADD;
	data->freq_start   = 1.0;
	data->freq_end     = 0.0;   // resolved to Nyquist below if unset
	data->n_freqs      = 30;
	data->n_settle     = 50;
	data->n_settle_cycles = 2.0;
	data->n_cycles     = 10.0;
	data->n_segments   = 7;
	data->amplitude    = 0.1;
	data->amplitude2   = 0.0;
	data->avoid_df     = 0.75;
	data->k_max_freq   = 30.0;
	data->output_mode  = BODE_OUTPUT_RAW;
	data->element      = 0;
	data->output_file  = xstrdup("bode_plot.pdf");
	data->fg_fd        = -1;
	data->fg_board     = 0;
	data->fg_channel   = 1;
	data->fg_level     = 1;
	data->fg_offset    = 0.0;

	if (!self->params) {
		log_error("No params object found.");
		return -1;
	}
	json_object_object_foreach(self->params, key, val) {
		if (key[0] == '_') {
			// comment
		} else if (!strcmp(key, "sample_rate")) {
			data->sample_rate = json_object_get_double(val);
			log_trace("sample_rate = %G Hz", data->sample_rate);
		} else if (!strcmp(key, "freq_start")) {
			data->freq_start = json_object_get_double(val);
			log_trace("freq_start = %G Hz", data->freq_start);
		} else if (!strcmp(key, "freq_end")) {
			data->freq_end = json_object_get_double(val);
			log_trace("freq_end = %G Hz", data->freq_end);
		} else if (!strcmp(key, "n_freqs")) {
			data->n_freqs = json_object_get_uint64(val);
			log_trace("n_freqs = %zu", data->n_freqs);
		} else if (!strcmp(key, "n_settle")) {
			data->n_settle = json_object_get_uint64(val);
			log_trace("n_settle = %zu", data->n_settle);
		} else if (!strcmp(key, "n_settle_cycles")) {
			data->n_settle_cycles = json_object_get_double(val);
			log_trace("n_settle_cycles = %G", data->n_settle_cycles);
		} else if (!strcmp(key, "n_cycles")) {
			data->n_cycles = json_object_get_double(val);
			log_trace("n_cycles = %G", data->n_cycles);
		} else if (!strcmp(key, "n_segments")) {
			data->n_segments = json_object_get_uint64(val);
			log_trace("n_segments = %zu", data->n_segments);
		} else if (!strcmp(key, "amplitude")) {
			data->amplitude = json_object_get_double(val);
			log_trace("amplitude = %G", data->amplitude);
		} else if (!strcmp(key, "amplitude2")) {
			data->amplitude2 = json_object_get_double(val);
			log_trace("amplitude2 = %G", data->amplitude2);
		} else if (!strcmp(key, "avoid_freqs")) {
			if (json_object_is_type(val, json_type_array)) {
				size_t na = json_object_array_length(val);
				if (na > BODE_MAX_AVOID) na = BODE_MAX_AVOID;
				for (size_t i = 0; i < na; i++)
					data->avoid_freqs[i] =
						json_object_get_double(
						json_object_array_get_idx(
						val, i));
				data->n_avoid = na;
			}
		} else if (!strcmp(key, "avoid_df")) {
			data->avoid_df = fabs(json_object_get_double(val));
		} else if (!strcmp(key, "k_max_freq")) {
			data->k_max_freq = json_object_get_double(val);
		} else if (!strcmp(key, "config")) {
			xfree(data->config_note);
			data->config_note = xstrdup(
				json_object_get_string(val));
		} else if (!strcmp(key, "start_delay")) {
			data->start_delay = json_object_get_double(val);
			log_trace("start_delay = %G s", data->start_delay);
		} else if (!strcmp(key, "inject")) {
			const char *s = json_object_get_string(val);
			if (!strcmp(s, "add"))
				data->inject_mode = BODE_INJECT_ADD;
			else if (!strcmp(s, "replace"))
				data->inject_mode = BODE_INJECT_REPLACE;
			else
				log_error("Unknown inject mode \"%s\"; use \"add\" or \"replace\"", s);
			log_trace("inject = %s", s);
		} else if (!strcmp(key, "output")) {
			const char *s = json_object_get_string(val);
			if (!strcmp(s, "raw"))
				data->output_mode = BODE_OUTPUT_RAW;
			else if (!strcmp(s, "diff"))
				data->output_mode = BODE_OUTPUT_DIFF;
			else
				log_error("Unknown output mode \"%s\"; use \"raw\" or \"diff\"", s);
			log_trace("output = %s", s);
		} else if (!strcmp(key, "element")) {
			data->element = json_object_get_uint64(val);
			log_trace("element = %zu", data->element);
		} else if (!strcmp(key, "freq_linear")) {
			data->freq_linear = json_object_get_boolean(val);
			log_trace("freq_linear = %s", data->freq_linear ? "true" : "false");
		} else if (!strcmp(key, "output_file")) {
			xfree(data->output_file);
			data->output_file = xstrdup(json_object_get_string(val));
			log_trace("output_file = %s", data->output_file);
		} else if (!strcmp(key, "fg_port")) {
			xfree(data->fg_port);
			data->fg_port = xstrdup(json_object_get_string(val));
			log_trace("fg_port = %s", data->fg_port);
		} else if (!strcmp(key, "fg_board")) {
			data->fg_board = (int)json_object_get_int64(val);
			log_trace("fg_board = %d", data->fg_board);
		} else if (!strcmp(key, "fg_channel")) {
			data->fg_channel = (int)json_object_get_int64(val);
			log_trace("fg_channel = %d", data->fg_channel);
		} else if (!strcmp(key, "fg_level")) {
			data->fg_level = (int)json_object_get_int64(val);
			log_trace("fg_level = %d", data->fg_level);
		} else if (!strcmp(key, "fg_offset")) {
			data->fg_offset = json_object_get_double(val);
			log_trace("fg_offset = %g V", data->fg_offset);
		} else {
			log_warn("Unknown parameter \"%s\"", key);
		}
	}

	if (!data->sample_rate) {
		log_error("sample_rate is required.");
		return -1;
	}
	if (!data->n_freqs) {
		log_error("n_freqs must be > 0.");
		return -1;
	}
	if (!data->freq_end)
		data->freq_end = data->sample_rate / 2.0;
	if (data->freq_start <= 0.0) {
		log_error("freq_start must be > 0.");
		return -1;
	}
	if (data->n_segments < 3) {
		log_warn("n_segments %zu < 3 gives meaningless error bars; "
			"using 3", data->n_segments);
		data->n_segments = 3;
	}
	if (data->n_segments > BODE_MAX_SEGMENTS)
		data->n_segments = BODE_MAX_SEGMENTS;
	// integer cycles per segment: boundaries land on excitation-phase
	// crossings so each segment integrates whole cycles (no partial-cycle
	// DC/line leakage), and n_segments phasors give the error bars
	data->cycles_per_seg = round(data->n_cycles / (double)data->n_segments);
	if (data->cycles_per_seg < 1.0) data->cycles_per_seg = 1.0;
	double actual = data->cycles_per_seg * (double)data->n_segments;
	if (fabs(actual - data->n_cycles) > 0.5)
		log_info("bode_plot: n_cycles %G rounded to %G (%zu segments "
			"x %G cycles)", data->n_cycles, actual,
			data->n_segments, data->cycles_per_seg);
	data->n_cycles = actual;
	if (data->fg_port && data->amplitude2 > 0.0) {
		log_warn("bode_plot: amplitude2 is a software-injection "
			"feature; hardware FG amplitude is quantized -- "
			"disabling the linearity pass");
		data->amplitude2 = 0.0;
	}

	data->freqs   = xcalloc(data->n_freqs, sizeof(double));
	data->mags    = xcalloc(data->n_freqs, sizeof(double));
	data->phases  = xcalloc(data->n_freqs, sizeof(double));
	data->mag_sem = xcalloc(data->n_freqs, sizeof(double));
	data->ph_sem  = xcalloc(data->n_freqs, sizeof(double));
	data->quality = xcalloc(data->n_freqs, sizeof(double));
	data->mags2   = xcalloc(data->n_freqs, sizeof(double));
	data->phases2 = xcalloc(data->n_freqs, sizeof(double));

	// frequency points: linear or log spaced
	for (size_t k = 0; k < data->n_freqs; k++) {
		double t = (data->n_freqs > 1) ? (double)k / (data->n_freqs - 1) : 0.0;
		if (data->freq_linear) {
			data->freqs[k] = data->freq_start + t * (data->freq_end - data->freq_start);
		} else {
			double log_start = log10(data->freq_start);
			double log_end   = log10(data->freq_end);
			data->freqs[k] = pow(10.0, log_start + t * (log_end - log_start));
		}
	}
	// nudge grid points off known disturbance lines: a line inside the
	// lock-in's leakage skirt biases the point coherently (segments do not
	// average it out), so move the injection avoid_df away instead
	for (size_t k = 0; k < data->n_freqs; k++) {
		for (size_t i = 0; i < data->n_avoid; i++) {
			double d = data->freqs[k] - data->avoid_freqs[i];
			if (fabs(d) < data->avoid_df) {
				double moved = data->avoid_freqs[i]
					+ (d >= 0.0 ? data->avoid_df
						: -data->avoid_df);
				log_info("bode_plot: moved sweep point "
					"%.2f -> %.2f Hz (%.2f Hz line)",
					data->freqs[k], moved,
					data->avoid_freqs[i]);
				data->freqs[k] = moved;
			}
		}
	}

	data->freq_idx   = 0;
	data->pass       = 0;
	data->step_count = 0;
	data->phase      = 0.0;
	data->prev_value = 0.0;
	advance_frequency(data);

	log_info("Sweeping %zu frequencies from %G to %G Hz (%zu segments x "
		"%G cycles each%s)", data->n_freqs, data->freq_start,
		data->freq_end, data->n_segments, data->cycles_per_seg,
		data->amplitude2 > 0.0 ? "; second linearity sweep" : "");

	if (data->fg_port) {
		data->fg_fd = fg_serial_open(data->fg_port);
		if (data->fg_fd < 0) {
			log_error("bode_plot: open FG port %s: %s",
				data->fg_port, strerror(errno));
			return -1;
		}
		char cmd[64];
		snprintf(cmd, sizeof cmd, "DAQC2.fgTYPE(%d, %d, 1)\n",
			data->fg_board, data->fg_channel);
		if (fg_send(data->fg_fd, cmd)) {
			log_error("bode_plot: FG init failed: %s", strerror(errno));
			return -1;
		}
		snprintf(cmd, sizeof cmd, "DAQC2.fgLEVEL(%d, %d, %d)\n",
			data->fg_board, data->fg_channel, data->fg_level);
		fg_send(data->fg_fd, cmd);
		if (data->fg_offset > 0.0) {
			// setDAC on the corresponding DAC channel sets the DC center;
			// the FG oscillates around that value when fgON is called.
			// FG channel N maps to DAC channel N-1.
			snprintf(cmd, sizeof cmd, "DAQC2.setDAC(%d, %d, %.4f)\n",
				data->fg_board, data->fg_channel - 1, data->fg_offset);
			fg_send(data->fg_fd, cmd);
		}
		log_info("bode_plot: hardware FG on %s board=%d channel=%d level=%d offset=%.3fV",
			data->fg_port, data->fg_board, data->fg_channel, data->fg_level,
			data->fg_offset > 0.0 ? data->fg_offset : 2.048);
	}

	self->type_in  = AYLP_T_VECTOR | AYLP_T_MATRIX;
	self->units_in = AYLP_U_ANY;
	self->type_out = AYLP_T_UNCHANGED;
	self->units_out = AYLP_U_UNCHANGED;

	return 0;
}


int bode_plot_proc(struct aylp_device *self, struct aylp_state *state)
{
	struct aylp_bode_plot_data *data = self->device_data;

	if (data->freq_idx >= data->n_freqs) {
		aylp_set_done(&state->header);
		return 0;
	}

	// 0. Hold before the sweep: inject 0 for start_delay seconds so a
	// downstream DAC bridge parks the driven channel at its offset (replace
	// mode) while e.g. coarse channels stabilize. The sweep clock starts
	// only after the hold expires.
	if (UNLIKELY(data->start_delay > 0.0 && !data->delay_done)) {
		double now = get_time_seconds();
		if (data->delay_start == 0.0) {
			data->delay_start = now;
			log_info("Holding injection at 0 for %g s before sweep",
				data->start_delay);
		}
		if (now - data->delay_start < data->start_delay) {
			if (data->fg_fd < 0)
				write_element(state, data->element, 0.0,
					data->inject_mode);
			return 0;
		}
		data->delay_done = true;
	}

	// 1. Initialize clock baseline cleanly on the very first sample of each frequency.
	// In FG mode, program the hardware FG first so start_time is captured right
	// after fgON completes — the FG is running from that moment and elapsed is the
	// correct phase reference.
	if (data->step_count == 0) {
		if (data->fg_fd >= 0) {
			double f = data->freqs[data->freq_idx];
			int f_int = (int)round(f);
			if (f_int < 10)    f_int = 10;
			if (f_int > 10000) f_int = 10000;
			char cmd[64];
			snprintf(cmd, sizeof cmd, "DAQC2.fgFREQ(%d, %d, %d)\n",
				data->fg_board, data->fg_channel, f_int);
			fg_send(data->fg_fd, cmd);
			snprintf(cmd, sizeof cmd, "DAQC2.fgON(%d, %d)\n",
				data->fg_board, data->fg_channel);
			fg_send(data->fg_fd, cmd);
			data->freqs[data->freq_idx] = (double)f_int;
		}
		double now = get_time_seconds();
		data->start_time = now;
		data->prev_time  = now;
	}

	// 2. read output signal y[n] before modifying state
	double raw = get_element(state, data->element);
	double y = (data->output_mode == BODE_OUTPUT_DIFF)
		? raw - data->prev_value
		: raw;
	data->prev_value = raw;

	// 3. Measure time delta
	double now = get_time_seconds();
	double dt = now - data->prev_time;
	data->prev_time = now;
	double elapsed = now - data->start_time;

	// 4. Compute phase and injection target
	double f = data->freqs[data->freq_idx];
	double A = data->pass ? data->amplitude2 : data->amplitude;
	double phi = 2.0 * M_PI * f * elapsed;
	double u   = A * sin(phi);

	// 5. Correlate against prev_phi: y[n] is the response to u[n-1] which
	// was injected at prev_phi, not the current phi.
	double corr_phi = (data->step_count > 0) ? data->prev_phi : phi;
	double curr_cos_term = y * cos(corr_phi);
	double curr_sin_term = y * sin(corr_phi);

	// 6. Integrate in integer-cycle segments once settled. Settle is both
	// sample- and cycle-based: 50 samples is far less than a period at
	// low frequency, and integrating the switch transient biases the point.
	bool point_done = false;
	bool settled = data->step_count >= data->n_settle
		&& elapsed >= data->n_settle_cycles / f;
	if (settled) {
		if (!data->integrating) {
			data->integrating = true;
			data->seg_idx = 0;
			data->seg_Yc = data->seg_Ys = data->seg_T = 0.0;
			data->seg_start_phi = corr_phi;
			data->prev_cos_term = curr_cos_term;
			data->prev_sin_term = curr_sin_term;
		} else {
			data->seg_Yc += 0.5 * (data->prev_cos_term + curr_cos_term) * dt;
			data->seg_Ys += 0.5 * (data->prev_sin_term + curr_sin_term) * dt;
			data->seg_T  += dt;
			data->dt_sum += dt;
			data->dt_n   += 1;
			data->prev_cos_term = curr_cos_term;
			data->prev_sin_term = curr_sin_term;
			// segment boundary on the integer-cycle phase crossing
			if (corr_phi - data->seg_start_phi
					>= data->cycles_per_seg * 2.0*M_PI
					&& data->seg_T > 0.0) {
				size_t i = data->seg_idx;
				data->seg_re[i] = 2.0 * data->seg_Ys
					/ (A * data->seg_T);
				data->seg_im[i] = 2.0 * data->seg_Yc
					/ (A * data->seg_T);
				data->seg_dur[i] = data->seg_T;
				data->seg_idx++;
				data->seg_start_phi += data->cycles_per_seg
					* 2.0*M_PI;
				data->seg_Yc = data->seg_Ys = data->seg_T = 0.0;
				if (data->seg_idx >= data->n_segments)
					point_done = true;
			}
		}
	}

	// 7. write injection to pipeline state (software mode only; FG drives hardware directly)
	if (data->fg_fd < 0)
		write_element(state, data->element, u, data->inject_mode);
	data->phase    = phi;
	data->prev_phi = phi;
	data->step_count++;

	if (!point_done) return 0;

	// 8. frequency complete: time-weighted coherent mean of the segment
	// phasors is the point estimate; segment scatter gives the error bars
	size_t ns = data->n_segments;
	double Ttot = 0.0, mre = 0.0, mim = 0.0;
	for (size_t i = 0; i < ns; i++) Ttot += data->seg_dur[i];
	for (size_t i = 0; i < ns; i++) {
		double w = data->seg_dur[i] / Ttot;
		mre += w * data->seg_re[i];
		mim += w * data->seg_im[i];
	}
	double mag = hypot(mre, mim);
	double phase = atan2(mim, mre);
	double mag_var = 0.0, ph_var = 0.0, mag_sum = 0.0;
	for (size_t i = 0; i < ns; i++) {
		double mi = hypot(data->seg_re[i], data->seg_im[i]);
		double pi = atan2(data->seg_im[i], data->seg_re[i]);
		double dp = pi - phase;
		while (dp >  M_PI) dp -= 2.0*M_PI;
		while (dp < -M_PI) dp += 2.0*M_PI;
		mag_var += (mi - mag) * (mi - mag);
		ph_var  += dp * dp;
		mag_sum += mi;
	}
	mag_var /= (double)(ns - 1);
	ph_var  /= (double)(ns - 1);
	double mag_sem = sqrt(mag_var / (double)ns);
	double ph_sem  = sqrt(ph_var / (double)ns);
	double qual = mag_sum > 0.0 ? mag / (mag_sum / (double)ns) : 0.0;

	size_t k = data->freq_idx;
	if (data->pass == 0) {
		data->mags[k] = mag;
		data->phases[k] = phase;
		data->mag_sem[k] = mag_sem;
		data->ph_sem[k] = ph_sem;
		data->quality[k] = qual;
	} else {
		data->mags2[k] = mag;
		data->phases2[k] = phase;
	}
	log_info("%sf=%7.2f Hz  |H|=%8.4f +/- %.4f (%6.1f dB)  ph=%6.1f "
		"+/- %.1f deg  q=%.3f",
		data->pass ? "[amp2] " : "", f, mag, mag_sem,
		20.0*log10(mag + 1e-300), phase * 180.0 / M_PI,
		ph_sem * 180.0 / M_PI, qual);
	if (qual < 0.9)
		log_warn("bode_plot: LOW QUALITY point at %.2f Hz (q=%.3f) "
			"-- disturbance line or drift inside the window; "
			"treat with suspicion", f, qual);

	// 9. advance to next frequency / pass
	data->freq_idx++;
	data->step_count = 0;
	data->phase = 0.0;
	data->prev_phi = 0.0;
	if (data->freq_idx < data->n_freqs) {
		advance_frequency(data);
		return 0;
	}
	if (data->pass == 0 && data->amplitude2 > 0.0) {
		data->pass = 1;
		data->freq_idx = 0;
		advance_frequency(data);
		log_info("bode_plot: main sweep done; starting linearity "
			"sweep at amplitude %G", data->amplitude2);
		return 0;
	}

	// 10. End of sweep: fit, save, plot
	bode_fit(data);
	int ret = generate_plot(data);
	aylp_set_done(&state->header);
	return ret;
}


int bode_plot_fini(struct aylp_device *self)
{
	struct aylp_bode_plot_data *data = self->device_data;
	if (data->fg_fd >= 0) {
		char cmd[64];
		snprintf(cmd, sizeof cmd, "DAQC2.fgOFF(%d, %d)\n",
			data->fg_board, data->fg_channel);
		fg_send(data->fg_fd, cmd);
		close(data->fg_fd);
	}
	xfree(data->fg_port);
	xfree(data->config_note);
	xfree(data->freqs);
	xfree(data->mags);
	xfree(data->phases);
	xfree(data->mag_sem);
	xfree(data->ph_sem);
	xfree(data->quality);
	xfree(data->mags2);
	xfree(data->phases2);
	xfree(data->output_file);
	xfree(data);
	return 0;
}
