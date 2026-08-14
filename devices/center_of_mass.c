#include <math.h>
#include <pthread.h>
#include <string.h>
#include <time.h>

#include "anyloop.h"
#include "logging.h"
#include "thread_pool.h"
#include "center_of_mass.h"
#include "xalloc.h"


// TODO: optimize? for example, if the matrix is contiguous and 8x8 in size, the
// compiler can do some really nice optimizations here. but for variable sizes,
// and non-contiguous matrices ... it gets a little harder.
// With -O3, on my old laptop, I timed that this could do 10 million 8x8
// matrices in 4.302 s, whereas the naive solution with contiguous memory took
// 2.447 s. So even with this naive non-contiguous solution, we still get at
// worst 0.4 μs per matrix, which is ... not the worst. That's 43 μs for a 10x10
// grid of submatrices, with no multithreading. GCC is quite intelligent.
/** Center of mass task for gsl_matrix_uchar.
* Will write (y,x) coords (in that order) of the center of mass of the uchar
* matrix at src->mat to dst[0] and dst[1], subtracting src->threshold from
* each pixel (clamping to zero) before accumulating. */
void com_mat_uchar(struct aylp_com_src *src, double *dst)
{
	double y = 0.0, x = 0.0, s = 0.0;
	// take weighted average
	for (size_t i=0; i < src->mat.size1; i++) {
		for (size_t j=0; j < src->mat.size2; j++) {
			unsigned char raw = src->mat.data[i*src->mat.tda + j];
			unsigned char el = raw > src->threshold
				? raw - src->threshold : 0;
			y += i*el;
			x += j*el;
			s += el;
		}
	}
	// if all pixels are at or below threshold, output centre
	if (!s) { dst[0] = 0.0; dst[1] = 0.0; return; }
	// set final values, scaling from 0:size-1 (given by y/s or x/s) to -1:1
	dst[0] = -1.0 + 2*y/(s*(src->mat.size1-1));
	dst[1] = -1.0 + 2*x/(s*(src->mat.size2-1));
}


// This is nice and fast, but doesn't work for us, because we want to find the
// com of a slice of non-contiguous memory....
// TODO: actually, it seems like contiguous matrices are one of the most likely
// scenarios in practice, e.g. if we're reading from a camera. Consider
// switching to this at runtime!
// /** Hardcoded 8x8 center of mass task so the compiler can optimize it.
// * Will look at 64 uchars at src, and write the output as two doubles in dst
// * in (y,x) order. */
// void com_8x8(unsigned char *src, double *dst)
// {
// 	double y = 0.0, x = 0.0, s = 0.0;
// 	for (char i = 0; i < 8; i++) {
// 		for (char j = 0; j < 8; j++) {
// 			y += i * src->data[8*i+j];
// 			x += j * src->data[8*i+j];
// 			s += src->data[8*i+j];
// 		}
// 	}
// 	dst[0] = y / s;
// 	dst[1] = x / s;
// }


int center_of_mass_init(struct aylp_device *self)
{
	int err;
	self->device_data = xcalloc(1, sizeof(struct aylp_center_of_mass_data));
	struct aylp_center_of_mass_data *data = self->device_data;

	// default params
	data->region_height = 0;
	data->region_width = 0;
	data->thread_count = 1;
	data->threshold = 0;
	data->min_peak = 0;	// 0 => any nonzero sum counts as a beam
	data->ref_cut = 0.0;		// 0 => no partial-beam test
	data->ref_warmup = 200;
	data->ref_rate = 0.01;
	data->ref_floor = 0.25;
	data->had_beam = true;	// so the first loss logs
	data->track = false;
	data->registration = false;
	data->init_y = -1;	// <0 => acquire from the brightest pixel
	data->init_x = -1;
	data->reacquire_after = 30;
	data->acquire_seconds = 0.0;	// no acquisition phase by default
	data->acquire_height = 0;	// 0 => the whole image
	data->acquire_width = 0;
	// parse parameters
	if (!self->params) {
		log_error("No params object found.");
		return -1;
	}
	json_object_object_foreach(self->params, key, val) {
		if (key[0] == '_') {
			// keys starting with _ are comments
		} else if (!strcmp(key, "region_height")) {
			data->auto_region_height=json_object_is_type(val,json_type_string)
				&& !strcmp(json_object_get_string(val),"auto");
			data->region_height=data->auto_region_height?0:
				json_object_get_uint64(val);
			log_trace("region_height = %zu", data->region_height);
		} else if (!strcmp(key, "region_width")) {
			data->auto_region_width=json_object_is_type(val,json_type_string)
				&& !strcmp(json_object_get_string(val),"auto");
			data->region_width=data->auto_region_width?0:
				json_object_get_uint64(val);
			log_trace("region_width = %zu", data->region_width);
		} else if (!strcmp(key, "thread_count")) {
			data->thread_count = json_object_get_uint64(val);
			if (data->thread_count == 0) {
				log_error("Correcting 0 threads to 1 thread");
				data->thread_count = 1;
			}
			log_trace("thread_count = %zu", data->thread_count);
		} else if (!strcmp(key, "threshold")) {
			data->auto_threshold=json_object_is_type(val,json_type_string)
				&& !strcmp(json_object_get_string(val),"auto");
			data->threshold=data->auto_threshold?0:
				(unsigned char)json_object_get_int(val);
			log_trace("threshold = %u", data->threshold);
		} else if (!strcmp(key, "min_peak")) {
			data->auto_min_peak=json_object_is_type(val,json_type_string)
				&& !strcmp(json_object_get_string(val),"auto");
			data->min_peak=data->auto_min_peak?0:
				(unsigned char)json_object_get_int(val);
			log_trace("min_peak = %u", data->min_peak);
		} else if (!strcmp(key, "ref_cut")) {
			data->ref_cut = json_object_get_double(val);
			if (data->ref_cut < 0.0 || data->ref_cut > 1.0) {
				log_error("ref_cut must be in [0, 1]; it is a "
					"fraction of the frame's typical row"
				);
				return -1;
			}
			log_trace("ref_cut = %G", data->ref_cut);
		} else if (!strcmp(key, "ref_warmup")) {
			data->ref_warmup = json_object_get_uint64(val);
			if (!data->ref_warmup) {
				log_error("ref_warmup must be nonzero");
				return -1;
			}
			log_trace("ref_warmup = %zu", data->ref_warmup);
		} else if (!strcmp(key, "ref_rate")) {
			data->ref_rate = json_object_get_double(val);
			if (data->ref_rate <= 0.0 || data->ref_rate > 1.0) {
				log_error("ref_rate must be in (0, 1]");
				return -1;
			}
			log_trace("ref_rate = %G", data->ref_rate);
		} else if (!strcmp(key, "ref_floor")) {
			data->ref_floor = json_object_get_double(val);
			if (data->ref_floor <= 0.0 || data->ref_floor >= 1.0) {
				log_error("ref_floor must be in (0, 1)");
				return -1;
			}
			log_trace("ref_floor = %G", data->ref_floor);
		} else if (!strcmp(key, "track")) {
			data->track = json_object_get_boolean(val);
			log_trace("track = %d", data->track);
		} else if (!strcmp(key, "registration")) {
			data->registration = json_object_get_boolean(val);
		} else if (!strncmp(key, "registration_", 13)) {
			log_warn("%s is obsolete: registration now derives sampling, "
				"noise, conditioning and robust scale automatically", key);
		} else if (!strcmp(key, "init_y")) {
			data->init_y=json_object_is_type(val,json_type_string)
				&& !strcmp(json_object_get_string(val),"auto") ? -1
				: json_object_get_int64(val);
			log_trace("init_y = %ld", data->init_y);
		} else if (!strcmp(key, "init_x")) {
			data->init_x=json_object_is_type(val,json_type_string)
				&& !strcmp(json_object_get_string(val),"auto") ? -1
				: json_object_get_int64(val);
			log_trace("init_x = %ld", data->init_x);
		} else if (!strcmp(key, "reacquire_after")) {
			data->auto_reacquire=json_object_is_type(val,json_type_string)
				&& !strcmp(json_object_get_string(val),"auto");
			data->reacquire_after=data->auto_reacquire?30:
				json_object_get_uint64(val);
			if (!data->reacquire_after) {
				log_error("reacquire_after must be nonzero");
				return -1;
			}
			log_trace("reacquire_after = %zu", data->reacquire_after);
		} else if (!strcmp(key, "acquire_seconds")) {
			data->acquire_seconds = json_object_get_double(val);
			if (data->acquire_seconds < 0.0) {
				log_error("acquire_seconds must be >= 0");
				return -1;
			}
			log_trace("acquire_seconds = %G", data->acquire_seconds);
		} else if (!strcmp(key, "acquire_height")) {
			data->acquire_height = json_object_get_uint64(val);
			log_trace("acquire_height = %zu", data->acquire_height);
		} else if (!strcmp(key, "acquire_width")) {
			data->acquire_width = json_object_get_uint64(val);
			log_trace("acquire_width = %zu", data->acquire_width);
		} else {
			log_warn("Unknown parameter \"%s\"", key);
		}
	}
	if ((!data->region_height && !data->auto_region_height)
	|| (!data->region_width && !data->auto_region_width)) {
		log_error("You must provide nonzero region_height and "
			"region_width params"
		);
		return -1;
	}
	if ((data->init_y < 0) != (data->init_x < 0)) {
		log_error("Provide both init_y and init_x, or neither");
		return -1;
	}
	if (data->min_peak && !data->track) {
		log_warn("min_peak only applies in track mode; ignoring it");
		data->min_peak = 0;
	}
	if (data->min_peak && data->min_peak <= data->threshold) {
		// a pixel has to clear `threshold` to weigh anything at all, so a
		// min_peak at or below it can never reject a frame the old code
		// would have accepted -- it is silently doing nothing
		log_warn("min_peak %u is at or below threshold %u, so it can "
			"never reject a frame; raise it above the background "
			"peak of a beamless frame", data->min_peak,
			data->threshold
		);
	}
	if (data->ref_cut > 0.0 && !data->track) {
		log_warn("ref_cut only applies in track mode; ignoring it");
		data->ref_cut = 0.0;
	}
	if (data->registration && !data->track) {
		log_error("registration requires track mode");
		return -1;
	}
	if (data->registration && !data->auto_region_height
			&& !data->auto_region_width && (data->region_height < 7
			|| data->region_width < 7)) {
		log_error("registration needs a tracking window at least 7x7");
		return -1;
	}
	if (data->ref_cut > 0.0 && !data->auto_region_height
			&& !data->auto_region_width) {
		if (data->region_height < 4) {
			// with three rows or fewer there is no meaningful
			// profile to compare against, and after ref_floor
			// trims the weak ones there would be nothing left
			log_error("ref_cut needs region_height >= 4 to have a "
				"row profile to compare; it is %zu",
				data->region_height
			);
			return -1;
		}
		data->ref = xcalloc(data->region_height, sizeof(double));
		data->rows = xcalloc(data->region_height, sizeof(double));
		data->rho = xcalloc(data->region_height, sizeof(double));
	}

	if (data->track) {
		if (data->registration && !data->auto_region_height
				&& !data->auto_region_width) {
			data->registration_ref = xcalloc(data->region_height
				* data->region_width, sizeof(double));
			// Every pixel with a two-pixel derivative margin is a valid
			// candidate.  Size the work arrays from the configured window:
			// a fixed cap makes feature selection depend on ROI size.
			data->registration_sample_capacity =
				(data->region_height - 4) * (data->region_width - 4);
			data->registration_sample_y = xcalloc(
				data->registration_sample_capacity,
				sizeof(size_t));
			data->registration_sample_x = xcalloc(
				data->registration_sample_capacity,
				sizeof(size_t));
			data->registration_sample_ref = xcalloc(
				data->registration_sample_capacity, sizeof(double));
			data->registration_sample_gy = xcalloc(
				data->registration_sample_capacity, sizeof(double));
			data->registration_sample_gx = xcalloc(
				data->registration_sample_capacity, sizeof(double));
			data->registration_residuals = xcalloc(
				data->registration_sample_capacity, sizeof(double));
			data->registration_errors = xcalloc(
				data->registration_sample_capacity, sizeof(double));
		}
		// the tracking window is a single region by definition, so there
		// is nothing to hand out to a thread pool
		if (data->thread_count > 1) {
			log_warn("track mode uses one region; "
				"ignoring thread_count = %zu", data->thread_count
			);
			data->thread_count = 1;
		}
		self->proc = &center_of_mass_proc_track;
		self->fini = &center_of_mass_fini;
	} else if (data->thread_count > 1) {
		// start threads
		data->threads = xmalloc(data->thread_count * sizeof(pthread_t));
		for (size_t t = 0; t < data->thread_count; t++) {
			err = pthread_create(&data->threads[t],
				0, task_runner, &data->queue
			);
			if (err) {
				log_error("Couldn't create pthread: %s",
					strerror(err)
				);
				return -1;
			}
		}
		log_info("Started %zu threads", data->thread_count);
		self->proc = &center_of_mass_proc_threaded;
		self->fini = &center_of_mass_fini_threaded;
	} else {
		// no threading
		self->proc = &center_of_mass_proc;
		self->fini = &center_of_mass_fini;
	}

	// set types and units
	self->type_in = AYLP_T_MATRIX_UCHAR;
	self->units_in = AYLP_U_ANY;
	self->type_out = AYLP_T_VECTOR;
	self->units_out = AYLP_U_MINMAX;
	return 0;
}


int center_of_mass_proc(struct aylp_device *self, struct aylp_state *state)
{
	struct aylp_center_of_mass_data *data = self->device_data;
	size_t max_y = state->matrix_uchar->size1;
	size_t max_x = state->matrix_uchar->size2;
	size_t y_subap_count = max_y / data->region_height;
	size_t x_subap_count = max_x / data->region_width;
	size_t subap_count = y_subap_count * x_subap_count;
	if (UNLIKELY(!subap_count)) {
		log_error("Refusing to process zero subapertures; "
			"region size is %zu by %zu but image is %zu by %zu",
			data->region_height, data->region_width, max_y, max_x
		);
		return -1;
	}
	// allocate the com vector if needed
	if (UNLIKELY(!data->com || data->com->size < subap_count*2)) {
		xfree_type(gsl_vector, data->com);
		data->com = xmalloc_type(gsl_vector, subap_count*2);
	}

	size_t n = 0;
	for (size_t i=0; i < y_subap_count; i++) {
		for (size_t j=0; j < x_subap_count; j++) {
			double y = 0.0, x = 0.0, s = 0.0;
			// inlined version of com_mat_uchar basically
			for (size_t l=0; l < data->region_height; l++) {
				for (size_t m=0; m < data->region_width; m++) {
					unsigned char raw;
					raw = state->matrix_uchar->data[
						(i*data->region_height + l)
						* state->matrix_uchar->tda
						+ j*data->region_width + m
					];
					unsigned char el = raw > data->threshold
						? raw - data->threshold : 0;
					y += l*el;
					x += m*el;
					s += el;
				}
			}
			if (!s) {
				data->com->data[2*n] = 0.0;
				data->com->data[2*n+1] = 0.0;
			} else {
				data->com->data[2*n] = -1.0
					+ 2*y/(s*(data->region_height-1));
				data->com->data[2*n+1] = -1.0
					+ 2*x/(s*(data->region_width-1));
			}
			n += 1;
		}
	}

	// zero-copy update of pipeline state
	state->vector = data->com;
	// housekeeping on the header
	state->header.type = self->type_out;
	state->header.units = self->units_out;
	state->header.log_dim.y = data->com->size;
	state->header.log_dim.x = 1;
	return 0;
}


int center_of_mass_fini(struct aylp_device *self)
{
	struct aylp_center_of_mass_data *data = self->device_data;
	// The reject rate is the number that decides whether gating chopped
	// frames is workable at all, and it is the one thing the rate-limited
	// per-episode warnings cannot tell you. Every held frame hands the loop
	// the previous frame's error again, so the integrator keeps winding on
	// stale data -- a few percent is free, a third is added loop delay.
	if (data && data->track && data->n_held) {
		log_info("center_of_mass: held %zu of %zu frames (%.2f%%) over "
			"%zu episodes", data->n_held, data->n_frames,
			100.0 * data->n_held / data->n_frames, data->n_episodes
		);
	}
	if (data) {
		xfree(data->ref);
		xfree(data->rows);
		xfree(data->rho);
		xfree(data->registration_ref);
		xfree(data->registration_sample_y);
		xfree(data->registration_sample_x);
		xfree(data->registration_sample_ref);
		xfree(data->registration_sample_gy);
		xfree(data->registration_sample_gx);
		xfree(data->registration_residuals);
		xfree(data->registration_errors);
	}
	xfree(self->device_data);
	return 0;
}


/** Place the tracking window on the brightest pixel of the whole image.
* Used to acquire on the first frame, and to recover after the beam has been
* lost for reacquire_after frames. Note that this locks onto whatever is
* brightest — if a stray reflection outpeaks the beam, pass init_y/init_x. */
static bool acquire_window(
	struct aylp_center_of_mass_data *data, gsl_matrix_uchar *img
) {
	unsigned char best = 0;
	unsigned char auto_gate=0;
	if(data->auto_threshold||data->auto_min_peak) {
		size_t bins[256]={0},n=0,acc=0,q=0;
		for(size_t y=0;y<img->size1;y++)for(size_t x=0;x<img->size2;x++)
			if(!y||!x||y+1==img->size1||x+1==img->size2) {
				bins[img->data[y*img->tda+x]]++;n++;
			}
		while(q<255&&(acc+=bins[q])<(95*n+99)/100)q++;
		if(data->auto_threshold)data->threshold=(unsigned char)q;
		auto_gate=(unsigned char)fmin(255,q+3);
	}
	size_t by = img->size1 / 2;
	size_t bx = img->size2 / 2;
	for (size_t i = 0; i < img->size1; i++) {
		for (size_t j = 0; j < img->size2; j++) {
			unsigned char v = img->data[i*img->tda + j];
			if (v > best) { best = v; by = i; bx = j; }
		}
	}
	// The brightest pixel of a beamless frame is just the loudest speck of
	// read noise, and it lands somewhere new every frame. Moving the window
	// there is worse than not moving at all: if the beam returns where it
	// left, an unmoved window is already on it. Park at the image centre
	// only when we have never had a placement to keep.
	if ((data->min_peak && best < data->min_peak)
			|| (data->auto_min_peak && best < auto_gate)) {
		if (!data->acquired) {
			data->win_y = img->size1 / 2;
			data->win_x = img->size2 / 2;
		}
		return false;
	}
	data->win_y = by;
	data->win_x = bx;
	return true;
}


/** Slide the window centre so that a win_h by win_w window lies fully inside the
* image. Callers must already have checked that the window is no bigger than the
* image. */
static void clamp_window(
	struct aylp_center_of_mass_data *data,
	size_t win_h, size_t win_w, size_t max_y, size_t max_x
) {
	size_t half_y = win_h / 2;
	size_t half_x = win_w / 2;
	if (data->win_y < half_y) data->win_y = half_y;
	if (data->win_x < half_x) data->win_x = half_x;
	if (data->win_y > max_y - win_h + half_y)
		data->win_y = max_y - win_h + half_y;
	if (data->win_x > max_x - win_w + half_x)
		data->win_x = max_x - win_w + half_x;
}


static double com_monotonic_s(void)
{
	struct timespec t;
	clock_gettime(CLOCK_MONOTONIC, &t);
	return t.tv_sec + 1e-9 * t.tv_nsec;
}

/** Calibrate the two brightness gates from the window border. The border is
 * deliberately used instead of the whole image so the beam cannot raise its
 * own noise floor. The 95th percentile tolerates isolated hot pixels while
 * following ordinary sensor/background changes. */
static void com_auto_gates(
	struct aylp_center_of_mass_data *data, const gsl_matrix_uchar *img,
	size_t org_y, size_t org_x, size_t h, size_t w
)
{
	if(!data->auto_threshold&&!data->auto_min_peak)return;
	size_t bins[256]={0},n=0,acc=0,q=0;
	for(size_t y=0;y<h;y++)for(size_t x=0;x<w;x++)
		if(!y||!x||y+1==h||x+1==w) {
			bins[img->data[(org_y+y)*img->tda+org_x+x]]++;n++;
		}
	if(!n)return;
	while(q<255&&(acc+=bins[q])<(95*n+99)/100)q++;
	if(data->auto_threshold)data->threshold=(unsigned char)q;
	if(data->auto_min_peak) {
		double signal=data->auto_peak_samples
			? data->auto_mean_peak-q : 0.0;
		double gate=q+fmax(3.0,ceil(0.10*fmax(0.0,signal)));
		data->min_peak=(unsigned char)fmin(255.0,gate);
	}
}

static void com_auto_accept_peak(
	struct aylp_center_of_mass_data *data, unsigned char peak
)
{
	if(!data->auto_min_peak)return;
	if(!data->auto_peak_samples)data->auto_mean_peak=peak;
	else data->auto_mean_peak+=0.05*(peak-data->auto_mean_peak);
	data->auto_peak_samples++;
}


/** Result of the partial-beam test, for logging. */
struct com_cut {
	bool cut;	// some rows are lit and others are not
	double lo;	// dimmest significant row, relative to the reference
	double mid;	// the frame's typical row, relative to the reference
	size_t n_sig;	// rows that carried enough reference to be judged
};

/** Is this frame's beam cut across rows?
*
* Divides the frame's row profile by the learned reference profile. On a whole
* beam every row comes back at the same ratio -- whatever the frame's brightness
* happens to be relative to the reference, uniformly across rows -- so the
* comparison below is between the dimmest row and the frame's own typical row,
* and the frame's overall level cancels. On a cut beam the ratios split into two
* groups, near the frame's level where the shutter let light through and near
* zero where it did not, and the dimmest row falls far under the typical one.
*
* The median is the typical row rather than the mean or the max: on a frame cut
* clean in half the mean sits halfway between the two groups and the max is one
* noisy row, while the median stays with whichever group holds more rows, which
* is all this needs to see a gap. */
static struct com_cut com_check_cut(
	struct aylp_center_of_mass_data *data, size_t win_h
) {
	struct com_cut r = {0};
	double ref_max = 0.0;
	for (size_t i = 0; i < win_h; i++)
		if (data->ref[i] > ref_max) ref_max = data->ref[i];
	if (ref_max <= 0.0)
		return r;	// no reference yet; nothing to say
	double floor = data->ref_floor * ref_max;

	// ratios of the rows the reference says should carry real signal,
	// insertion-sorted as they are collected. win_h is a window height --
	// tens of rows at most -- so this is a handful of compares, and it keeps
	// the median exact rather than approximating it
	for (size_t i = 0; i < win_h; i++) {
		if (data->ref[i] < floor)
			continue;
		double v = data->rows[i] / data->ref[i];
		size_t k = r.n_sig++;
		while (k > 0 && data->rho[k-1] > v) {
			data->rho[k] = data->rho[k-1];
			k--;
		}
		data->rho[k] = v;
	}
	// three rows is the fewest that can show a dim one against a typical
	// one; below that the test would be comparing a row to itself
	if (r.n_sig < 3)
		return r;

	r.lo = data->rho[0];
	r.mid = data->rho[r.n_sig / 2];
	// mid <= 0 means every significant row lost its light while the window
	// as a whole still has some -- the beam is not where the reference says
	// it is, which is not a frame to steer on either
	r.cut = r.mid <= 0.0 || r.lo < data->ref_cut * r.mid;
	return r;
}


/** Fold this frame's row profile into the reference.
*
* Before the bootstrap finishes this takes the row-wise maximum, because the
* source may already be chopping while it runs: a cut only ever removes light,
* so each row's largest value over enough frames is its uncut value, while a
* mean would learn a blend of whole and cut beams and sit too low. Afterwards it
* is an EMA over accepted frames only, which follows slow drift in power, focus
* and alignment without letting cut frames pull it down toward themselves. */
static void com_update_ref(struct aylp_center_of_mass_data *data, size_t win_h)
{
	if (!data->ref_ready) {
		for (size_t i = 0; i < win_h; i++)
			if (data->rows[i] > data->ref[i])
				data->ref[i] = data->rows[i];
		if (++data->ref_seen >= data->ref_warmup) {
			data->ref_ready = true;
			log_info("center_of_mass: reference profile learned "
				"from %zu frames; partial-beam gate is now "
				"active", data->ref_seen
			);
		}
		return;
	}
	for (size_t i = 0; i < win_h; i++)
		data->ref[i] += data->ref_rate * (data->rows[i] - data->ref[i]);
}


/** Throw the learned reference away.
*
* Called when the window re-acquires: the beam has been gone long enough that
* the window is being replaced from scratch, and a profile learned at the old
* position describes a beam that may not be there any more. Relearning costs
* ref_warmup frames; keeping a stale reference costs every frame after. */
static void com_reset_ref(struct aylp_center_of_mass_data *data)
{
	if (!data->ref)
		return;
	for (size_t i = 0; i < data->region_height; i++)
		data->ref[i] = 0.0;
	data->ref_seen = 0;
	data->ref_ready = false;
}


/** Seconds between beam-lost/beam-back log lines. */
#define COM_LOG_INTERVAL 1.0

static double com_bilinear(
	const gsl_matrix_uchar *img, double y, double x
) {
	size_t y0 = (size_t)y, x0 = (size_t)x;
	double fy = y-y0, fx = x-x0;
	const unsigned char *a = img->data + y0*img->tda + x0;
	return (1-fy)*((1-fx)*a[0] + fx*a[1])
		+ fy*((1-fx)*a[img->tda] + fx*a[img->tda+1]);
}

/** Solve a dense 4x4 system with partial pivoting. */
static bool com_solve4(double a[4][4], double b[4], double x[4])
{
	for (size_t k=0; k<4; k++) {
		size_t p=k;
		for (size_t i=k+1; i<4; i++)
			if (fabs(a[i][k]) > fabs(a[p][k])) p=i;
		if (fabs(a[p][k]) < 1e-9) return false;
		if (p != k) {
			for (size_t j=k; j<4; j++) {
				double t=a[k][j]; a[k][j]=a[p][j]; a[p][j]=t;
			}
			double t=b[k]; b[k]=b[p]; b[p]=t;
		}
		for (size_t i=k+1; i<4; i++) {
			double q=a[i][k]/a[k][k];
			for (size_t j=k; j<4; j++) a[i][j]-=q*a[k][j];
			b[i]-=q*b[k];
		}
	}
	for (int i=3; i>=0; i--) {
		double s=b[i];
		for (size_t j=i+1; j<4; j++) s-=a[i][j]*x[j];
		x[i]=s/a[i][i];
	}
	return true;
}

static inline void com_normal_add(
	double n[4][4], double r[4], double gy, double gx, double ref,
	double e, double wt
)
{
	double wy=wt*gy,wxx=wt*gx,wr=wt*ref;
	n[0][0]+=wy*gy; n[0][1]+=wy*gx; n[0][2]+=wy*ref; n[0][3]+=wy;
	n[1][1]+=wxx*gx; n[1][2]+=wxx*ref; n[1][3]+=wxx;
	n[2][2]+=wr*ref; n[2][3]+=wr; n[3][3]+=wt;
	double we=wt*e;
	r[0]+=gy*we; r[1]+=gx*we; r[2]+=ref*we; r[3]+=we;
}

static inline void com_normal_mirror(double n[4][4])
{
	for(size_t i=1;i<4;i++)for(size_t j=0;j<i;j++)n[i][j]=n[j][i];
}

struct com_feature { double g, gy, gx; size_t y, x; };
static int com_feature_desc(const void *a, const void *b)
{
	double x=((const struct com_feature*)a)->g;
	double y=((const struct com_feature*)b)->g;
	return x>y?-1:x<y;
}
static int com_double_asc(const void *a, const void *b)
{
	double x=*(const double*)a,y=*(const double*)b;
	return x<y?-1:x>y;
}

/** In-place selection of the kth value. Residual scale only needs a median;
 * sorting every residual made the robust pass O(n log n) on every frame. */
static double com_select_kth(double *a, size_t n, size_t k)
{
	size_t lo=0, hi=n-1;
	while (lo < hi) {
		double pivot=a[lo+(hi-lo)/2];
		size_t i=lo, j=hi;
		for (;;) {
			while (a[i] < pivot) i++;
			while (a[j] > pivot) j--;
			if (i >= j) break;
			double t=a[i]; a[i++]=a[j]; a[j--]=t;
		}
		if (k <= j) hi=j;
		else lo=j+1;
	}
	return a[k];
}

/** Choose as many keyframe gradients as the measured noise and two-axis
* information require. The stop rule is expressed in predicted translation
* uncertainty, not pixels-per-beam, brightness, diameter or pattern class. */
static bool com_select_registration_samples(
	struct aylp_center_of_mass_data *data
) {
	size_t h=data->region_height,w=data->region_width,n=0;
	size_t cap=data->registration_sample_capacity;
	struct com_feature *f=xmalloc(cap*sizeof *f);
	double bg=0,n_bg=0;
	for(size_t y=0;y<h;y++)for(size_t x=0;x<w;x++)
		if(!y||!x||y+1==h||x+1==w) {
			bg+=data->registration_ref[y*w+x];n_bg++;
		}
	if(n_bg)bg/=n_bg;
	// Second differences vanish on real linear structure but retain sensor
	// noise. Their median supplies a pattern-independent noise estimate.
	size_t noise[1021]={0},nn=0;
	for(size_t y=2;y+2<h;y++) for(size_t x=2;x+2<w;x++) {
			double gy=.5*(data->registration_ref[(y+1)*w+x]
				-data->registration_ref[(y-1)*w+x]);
			double gx=.5*(data->registration_ref[y*w+x+1]
				-data->registration_ref[y*w+x-1]);
		struct com_feature z={gy*gy+gx*gx,gy,gx,y,x};
		if(n<cap) {
			size_t k=n++;f[k]=z;
			while(k){size_t p=(k-1)/2;if(f[p].g<=f[k].g)break;struct com_feature q=f[p];f[p]=f[k];f[k]=q;k=p;}
		} else if(z.g>f[0].g) {
			f[0]=z;
			for(size_t k=0;;){size_t a=2*k+1,b=a+1;if(a>=n)break;size_t m=b<n&&f[b].g<f[a].g?b:a;if(f[k].g<=f[m].g)break;struct com_feature q=f[k];f[k]=f[m];f[m]=q;k=m;}
		}
		double ly=data->registration_ref[(y+1)*w+x]
			-2*data->registration_ref[y*w+x]
			+data->registration_ref[(y-1)*w+x];
		double lx=data->registration_ref[y*w+x+1]
			-2*data->registration_ref[y*w+x]
			+data->registration_ref[y*w+x-1];
		size_t q=(size_t)fmin(1020.0,fabs(ly)+fabs(lx));noise[q]++;nn++;
	}
	size_t acc=0,med=0;while(med<1020&&(acc+=noise[med])<nn/2)med++;
	double sigma=fmax(0.5,med/1.349);
	data->registration_noise=sigma;
	qsort(f,n,sizeof *f,com_feature_desc);
	double yy=0,xx=0,yx=0; size_t keep=0;
	// A larger window contains more independent places where structure can
	// deform or be occluded. Require four square-root sample budgets so robust
	// fitting does not collapse to the same 16 brightest pixels at every ROI.
	// This grows sublinearly (112 -> 432, 384 -> 1520) and is derived entirely
	// from the available candidate population.
	size_t min_keep=4*(size_t)ceil(sqrt((double)cap));
	if(min_keep<16)min_keep=16;
	// 0.05 px is an output-precision contract, not a beam property. Stop as
	// soon as both eigen-directions of translation meet it.  The candidate
	// count comes from the window dimensions, not a preferred ROI size.
	double need=sigma*sigma/(0.05*0.05);
	data->registration_info_need=need;
	for(size_t k=0;k<n&&keep<cap;k++) {
		if(f[k].g<=sigma*sigma) break;
		if(fabs(data->registration_ref[f[k].y*w+f[k].x]-bg)
				<=fmax(2.0,3.0*sigma))continue;
		data->registration_sample_y[keep]=f[k].y;
		data->registration_sample_x[keep]=f[k].x;
		data->registration_sample_ref[keep]=
			data->registration_ref[f[k].y*w+f[k].x];
		data->registration_sample_gy[keep]=f[k].gy;
		data->registration_sample_gx[keep]=f[k].gx;
		keep++;
		yy+=f[k].gy*f[k].gy;xx+=f[k].gx*f[k].gx;yx+=f[k].gy*f[k].gx;
		double tr=yy+xx,disc=sqrt(fmax(0.0,(yy-xx)*(yy-xx)+4*yx*yx));
		if(keep>=min_keep && .5*(tr-disc)>=need) break;
	}
	xfree(f);data->registration_n_samples=keep;
	if(keep<16) return false;
	double tr=yy+xx,disc=sqrt(fmax(0.0,(yy-xx)*(yy-xx)+4*yx*yx));
	data->registration_condition=(tr+disc)/fmax(1e-12,tr-disc);
	return .5*(tr-disc)>=need;
}

static bool com_registration_overlap_ok(
	struct aylp_center_of_mass_data *data, const gsl_matrix_uchar *img,
	size_t org_y,size_t org_x
) {
	double yy=0,xx=0,yx=0;
	for(size_t k=0;k<data->registration_n_samples;k++) {
		size_t y=data->registration_sample_y[k],x=data->registration_sample_x[k];
		double cy=org_y+y+data->registration_y,cx=org_x+x+data->registration_x;
		// A handoff is cheap compared with losing correspondence. Keep every
		// selected informative feature inside the interpolation domain; this is
		// determined by the selected pattern itself, not a displacement limit.
		if(cy<1||cx<1||cy+1>=img->size1||cx+1>=img->size2)continue;
		double gy=data->registration_sample_gy[k];
		double gx=data->registration_sample_gx[k];
		yy+=gy*gy;xx+=gx*gx;yx+=gy*gx;
	}
	double tr=yy+xx,disc=sqrt(fmax(0.0,(yy-xx)*(yy-xx)+4*yx*yx));
	return .5*(tr-disc)>=data->registration_info_need;
}

static bool com_install_registration_keyframe(
	struct aylp_center_of_mass_data *data, const gsl_matrix_uchar *img,
	size_t org_y, size_t org_x, double abs_y, double abs_x
) {
	size_t h=data->region_height,w=data->region_width;
	for(size_t y=0;y<h;y++)for(size_t x=0;x<w;x++)
		data->registration_ref[y*w+x]=img->data[(org_y+y)*img->tda+org_x+x];
	data->registration_ref_y=abs_y;data->registration_ref_x=abs_x;
	data->registration_y=0.0;data->registration_x=0.0;
	return com_select_registration_samples(data);
}

/** Background-invariant absolute centre for the slow registration tether.
 * The registration fit treats offset as a nuisance parameter, so its absolute
 * reference must do the same: an unsubtracted frame pedestal pulls an ordinary
 * centroid toward the ROI centre and attenuates real motion. */
static bool com_registration_anchor_centroid(
	const gsl_matrix_uchar *img, size_t org_y, size_t org_x,
	size_t h, size_t w, double *cy, double *cx
)
{
	size_t n_border=2*h+2*w-4,nb=0,n_hot=0;
	double *border=xmalloc(n_border*sizeof *border);
	for(size_t y=0;y<h;y++)for(size_t x=0;x<w;x++)
		if(!y||!x||y+1==h||x+1==w) {
			border[nb++]=img->data[(org_y+y)*img->tda+org_x+x];
		}
	if(!nb){xfree(border);return false;}
	qsort(border,nb,sizeof *border,com_double_asc);
	double bg=border[nb/2];
	// A centroid is not an absolute observation when the target intersects the
	// ROI boundary: translating it changes how much flux is clipped. Detect
	// sustained edge signal against the robust border background and withhold
	// the tether; local registration remains valid on the visible structure.
	for(size_t k=0;k<nb;k++)if(border[k]>bg+5.0)n_hot++;
	xfree(border);
	if(n_hot>=2)return false;
	double s=0,sy=0,sx=0;
	for(size_t y=0;y<h;y++)for(size_t x=0;x<w;x++) {
		double v=img->data[(org_y+y)*img->tda+org_x+x]-bg;
		if(v<=2.0)continue;
		s+=v;sy+=v*(org_y+y);sx+=v*(org_x+x);
	}
	if(s<=0)return false;
	*cy=sy/s;*cx=sx/s;return true;
}

/** Brightness-affine Lucas--Kanade translation against a fixed keyframe.
* Gain and offset are fitted alongside y/x, so uniform illumination changes
* cannot become position. The sparse grid bounds cost independently of beam
* shape; the only requirement is nonzero spatial structure in both axes. */
static bool com_register(
	struct aylp_center_of_mass_data *data, const gsl_matrix_uchar *img,
	size_t org_y, size_t org_x
) {
	double d[4]={0};
	double trial_y=data->registration_y,trial_x=data->registration_x;
	size_t used=0;
	double span=fmin(data->region_height,data->region_width);
	// Re-linearize several times. Two updates were sufficient for broad smooth
	// beams but attenuated known ±4-pixel motion on compact real patterns.
	for(size_t iter=0;iter<1;iter++) {
		double n[4][4]={{0}},r[4]={0};used=0;
		for (size_t k=0;k<data->registration_n_samples;k++) {
			size_t y=data->registration_sample_y[k],x=data->registration_sample_x[k];
			double cy=org_y+y+trial_y,cx=org_x+x+trial_x;
			if (cy<1 || cx<1 || cy+1>=img->size1 || cx+1>=img->size2) continue;
			double cur=com_bilinear(img,cy,cx);
			double ref=data->registration_sample_ref[k];
			double rgy=data->registration_sample_gy[k];
			double rgx=data->registration_sample_gx[k];
			// Use keyframe gradients: an exposure boundary present only in the
			// destination is an outlier, not evidence of rigid image motion.
			double gy=rgy,gx=rgx,e=ref-cur;
			com_normal_add(n,r,gy,gx,ref,e,1.0);
			used++;
		}
		com_normal_mirror(n);
		if(used<16||!com_solve4(n,r,d)||!isfinite(d[0])||!isfinite(d[1])
				||hypot(d[0],d[1])>span/8.0)return false;
		trial_y+=d[0];trial_x+=d[1];
		if(hypot(d[0],d[1])<0.01)break;
	}

	// Median absolute residual estimates this frame's scale without a noise or
	// deformation knob. Tukey's 4.685-sigma consistency constant then removes
	// samples that do not share the dominant rigid motion.
	used=0;
	for(size_t k=0;k<data->registration_n_samples;k++) {
		size_t y=data->registration_sample_y[k],x=data->registration_sample_x[k];
		double cy=org_y+y+trial_y,cx=org_x+x+trial_x;
		if(cy<1||cx<1||cy+1>=img->size1||cx+1>=img->size2)continue;
		double ref=data->registration_sample_ref[k];
		double e=ref-com_bilinear(img,cy,cx)-d[2]*ref-d[3];
		data->registration_errors[k]=e;
		data->registration_residuals[used++]=fabs(e);
	}
	if(used<16)return false;
	data->registration_residual_scale=
		1.4826*com_select_kth(data->registration_residuals,used,used/2);
	double cut=4.685*fmax(0.25,data->registration_residual_scale);
	double nn[4][4]={{0}},rr[4]={0},dd[4]={0};size_t inliers=0;
	for(size_t k=0;k<data->registration_n_samples;k++) {
		size_t y=data->registration_sample_y[k],x=data->registration_sample_x[k];
		double cy=org_y+y+trial_y,cx=org_x+x+trial_x;
		if(cy<1||cx<1||cy+1>=img->size1||cx+1>=img->size2)continue;
		double ref=data->registration_sample_ref[k];
		double rgy=data->registration_sample_gy[k];
		double rgx=data->registration_sample_gx[k];
		double gy=rgy,gx=rgx;
		double e=data->registration_errors[k],q=fabs(e)/cut;
		double one_minus_q2=1-q*q;
		double wt=q<1?one_minus_q2*one_minus_q2:0;
		if(wt<=0)continue;
		inliers++;
		com_normal_add(nn,rr,gy,gx,ref,e,wt);
	}
	com_normal_mirror(nn);
	if(inliers<16||!com_solve4(nn,rr,dd)||!isfinite(dd[0])||!isfinite(dd[1])
			||hypot(dd[0],dd[1])>span/16.0)return false;
	data->registration_quality=(double)inliers/used;
	// Commit atomically.  A rejected frame must not poison the starting point
	// for the next one; this matters especially at the faster 112-pixel ROI,
	// where one failed solve is followed by appreciably more real motion.
	data->registration_y=trial_y+dd[0];
	data->registration_x=trial_x+dd[1];
	return true;
}

/** True at most once per COM_LOG_INTERVAL, so the beam-lost and beam-back
* transitions can be logged without flooding.
*
* These fire once per loss episode, which is the right cadence for a beam that
* occasionally drops out and completely the wrong one for a chopped source: a
* flux gate rejecting alternate frames produces an episode every few frames, and
* at loop rate that is hundreds of log lines a second going to a terminal in the
* middle of the control loop. The suppressed episodes are still counted, and
* fini reports the totals, so nothing is lost but the noise. */
static bool com_log_ok(struct aylp_center_of_mass_data *data)
{
	double now = com_monotonic_s();
	if (now - data->last_log_t < COM_LOG_INTERVAL)
		return false;
	data->last_log_t = now;
	return true;
}


/** Center of mass over a single window that follows the beam.
*
* Each frame the sum is taken over a region_height by region_width box centred
* on the previous frame's center of mass, so anything outside that box — a stray
* reflection elsewhere on the sensor, say — never enters the sum. The image
* itself is untouched, so a udp_sink placed ahead of this device still shows the
* whole frame.
*
* The output is normalized across the *whole image*, not the window. This is the
* important part: the window chases the beam, so a window-relative coordinate
* would sit near zero no matter where the beam actually was, and the loop would
* have no error signal to act on. Normalizing to the image keeps the setpoint at
* the image centre and keeps the error-per-pixel — hence the loop gain —
* independent of the window size.
*/
int center_of_mass_proc_track(struct aylp_device *self, struct aylp_state *state)
{
	struct aylp_center_of_mass_data *data = self->device_data;
	gsl_matrix_uchar *img = state->matrix_uchar;
	size_t max_y = img->size1;
	size_t max_x = img->size2;
	if (UNLIKELY(data->auto_region_height || data->auto_region_width)) {
		// A 243-pixel window preserves the bounded registration cost used by
		// the steering profiles; smaller camera ROIs use their full extent.
		if (data->auto_region_height)
			data->region_height = max_y < 243 ? max_y : 243;
		if (data->auto_region_width)
			data->region_width = max_x < 243 ? max_x : 243;
		data->auto_region_height = data->auto_region_width = false;
		if (data->registration) {
			data->registration_ref = xcalloc(data->region_height
				* data->region_width, sizeof(double));
			data->registration_sample_capacity =
				(data->region_height - 4) * (data->region_width - 4);
			data->registration_sample_y = xcalloc(
				data->registration_sample_capacity, sizeof(size_t));
			data->registration_sample_x = xcalloc(
				data->registration_sample_capacity, sizeof(size_t));
			data->registration_sample_ref = xcalloc(
				data->registration_sample_capacity, sizeof(double));
			data->registration_sample_gy = xcalloc(
				data->registration_sample_capacity, sizeof(double));
			data->registration_sample_gx = xcalloc(
				data->registration_sample_capacity, sizeof(double));
			data->registration_residuals = xcalloc(
				data->registration_sample_capacity, sizeof(double));
			data->registration_errors = xcalloc(
				data->registration_sample_capacity, sizeof(double));
		}
	}
	if (UNLIKELY(data->region_height > max_y
			|| data->region_width > max_x)) {
		log_error("Tracking window is %zu by %zu but image is only "
			"%zu by %zu", data->region_height, data->region_width,
			max_y, max_x
		);
		return -1;
	}
	// allocate the com vector if needed; track mode emits a single y,x pair
	if (UNLIKELY(!data->com || data->com->size < 2)) {
		xfree_type(gsl_vector, data->com);
		data->com = xmalloc_type(gsl_vector, 2);
	}
	if (UNLIKELY(!data->acquired)) {
		bool placed;
		if (data->init_y >= 0) {
			data->win_y = (size_t)data->init_y;
			data->win_x = (size_t)data->init_x;
			placed = true;
		} else {
			placed = acquire_window(data, img);
		}
		// Staying unacquired until a real beam shows up also keeps the
		// acquisition phase honest: its clock must start when there is
		// something to acquire, or it expires against a dark sensor and
		// narrows the window onto noise.
		if (placed) {
			data->acquired = true;
			data->acquiring = data->acquire_seconds > 0.0;
			data->acquire_t0 = com_monotonic_s();
			log_info("center_of_mass: acquired window at (%zu,%zu)%s",
				data->win_y, data->win_x,
				data->acquiring
					? "; starting wide acquisition phase" : ""
			);
		}
	}
	// During the acquisition phase the sum runs over a wider window (the whole
	// image by default). The point is that a wide window is flux-weighted over
	// everything in it, so the centroid is dragged toward whichever spot carries
	// the most light, rather than toward whichever spot the brightest-pixel scan
	// happened to hit first in raster order. When the phase ends we narrow onto
	// wherever that centroid settled.
	size_t win_h = data->region_height;
	size_t win_w = data->region_width;
	if (UNLIKELY(data->acquiring)) {
		win_h = data->acquire_height ? data->acquire_height : max_y;
		win_w = data->acquire_width ? data->acquire_width : max_x;
		if (win_h > max_y) win_h = max_y;
		if (win_w > max_x) win_w = max_x;
	}
	clamp_window(data, win_h, win_w, max_y, max_x);
	size_t org_y = data->win_y - win_h/2;
	size_t org_x = data->win_x - win_w/2;
	com_auto_gates(data,img,org_y,org_x,win_h,win_w);

	// The row profile is only kept when the partial-beam test is armed and
	// the window is at its configured size. During the acquisition phase the
	// window is a different (wider) shape, so a profile taken then does not
	// describe the same thing and must not be mixed into the reference.
	bool profile = data->ref_cut > 0.0 && !data->acquiring;

	double y = 0.0, x = 0.0, s = 0.0;
	unsigned char peak = 0;
	if (data->registration_ready && !data->acquiring && !profile) {
		// Once registration has an absolute origin, the flux centroid is not
		// an output. The normal path only needs the peak for the signal gate;
		// avoiding three floating-point moments over the whole ROI is material
		// at 384x384 and lets the compiler vectorize this contiguous max scan.
		for (size_t l=0;l<win_h;l++) {
			const unsigned char *p=img->data+(org_y+l)*img->tda+org_x;
			for(size_t m=0;m<win_w;m++) if(p[m]>peak)peak=p[m];
		}
		s=peak>data->threshold;
	} else {
		for (size_t l = 0; l < win_h; l++) {
			double row = 0.0;
			for (size_t m = 0; m < win_w; m++) {
				unsigned char raw = img->data[
					(org_y + l) * img->tda + org_x + m
				];
				if (raw > peak) peak = raw;
				unsigned char el = raw > data->threshold
					? raw - data->threshold : 0;
				// accumulate in image coordinates, not window ones
				y += (org_y + l)*el;
				x += (org_x + m)*el;
				s += el;
				row += el;
			}
			// win_h == region_height whenever `profile` is set, so this
			// cannot run off the end of the array
			if (profile) data->rows[l] = row;
		}
	}

	// "Is the beam in frame?" is a question about brightness, and a nonzero
	// sum does not answer it: sensor read noise running a few counts above
	// `threshold` gives a perfectly nonzero sum whose centroid is noise and
	// wanders the full width of the window. Only a real peak separates a
	// beam from a dark frame -- with min_peak 0 this is exactly the old
	// s != 0 test, so configs that do not set it are unaffected.
	//
	// The three tests fail for different reasons and the log line has to say
	// which, so keep them apart rather than folding them into one bool: an
	// empty sum means nothing cleared `threshold`, a low peak means the
	// window holds no bright pixel, and a split row profile with a healthy
	// peak means the beam is only partly there -- a rolling shutter that
	// caught a chopped source mid-cycle. Reporting the wrong one sends you
	// tuning the wrong parameter.
	//
	// The partial-beam test is skipped until the reference has been learned,
	// so the first ref_warmup frames go through on the brightness tests
	// alone. That is deliberate: a gate is worth less than the reference it
	// would be guessing at, and at loop rate the warmup is milliseconds.
	struct com_cut cut = {0};
	const char *why = 0;
	if (UNLIKELY(s == 0.0))
		why = "sum";
	else if (UNLIKELY(peak < data->min_peak))
		why = "peak";
	else if (profile && data->ref_ready) {
		cut = com_check_cut(data, win_h);
		if (UNLIKELY(cut.cut))
			why = "cut";
	}
	bool beam = !why;
	// only whole frames teach the reference what a whole beam looks like
	if (profile && beam)
		com_update_ref(data, win_h);
	double abs_y = s > 0.0 ? y/s : 0.0, abs_x = s > 0.0 ? x/s : 0.0;
	if (beam && data->registration && !data->acquiring) {
		if (!data->registration_ready) {
			// Start the absolute coordinate from the same background-invariant
			// measurement used by the long-term tether. Otherwise a frame pedestal
			// becomes a permanent initial offset that the tether later appears to
			// "drift" away from.
			double anchor_y,anchor_x;
			if(com_registration_anchor_centroid(img,org_y,org_x,
					win_h,win_w,&anchor_y,&anchor_x)) {
				abs_y=anchor_y;abs_x=anchor_x;
			}
			data->registration_ready=com_install_registration_keyframe(
				data,img,org_y,org_x,abs_y,abs_x);
			if(!data->registration_ready) { why="registration"; beam=false; }
		} else {
			if(!com_register(data,img,org_y,org_x)) {
				double reacquire_y,reacquire_x;
				if(com_registration_anchor_centroid(img,org_y,org_x,
						win_h,win_w,&reacquire_y,&reacquire_x)) {
					abs_y=reacquire_y;abs_x=reacquire_x;
					data->registration_rolls++;
					data->registration_ready=com_install_registration_keyframe(
						data,img,org_y,org_x,abs_y,abs_x);
					if(!data->registration_ready){why="registration";beam=false;}
				} else {
					why="registration";
					beam=false;
				}
			} else {
			abs_y=data->registration_ref_y+data->registration_y;
			abs_x=data->registration_ref_x+data->registration_x;
			// Registration is a locally precise displacement measurement, but a
			// rolling chain is dead reckoning: a tiny signed error in every link
			// accumulates without bound. Keep a robust history of its disagreement
			// with the independent flux centre. The median rejects isolated rolling
			// exposures and asymmetric shape excursions; only a persistent offset
			// becomes a slow absolute correction.
			double anchor_cy=0,anchor_cx=0;
			bool overlap_ok=com_registration_overlap_ok(data,img,org_y,org_x);
			bool have_anchor=false;
			// The absolute tether is irrelevant until a keyframe handoff has
			// occurred. Avoiding this full-window centroid scan is the normal
			// path's largest speedup. We still obtain it before any handoff.
			if(data->registration_rolls || !overlap_ok)
				have_anchor=com_registration_anchor_centroid(
					img,org_y,org_x,win_h,win_w,&anchor_cy,&anchor_cx);
			size_t ap=data->registration_anchor_pos;
			if(have_anchor&&data->registration_rolls) {
				data->registration_anchor_y[ap]=anchor_cy-abs_y;
				data->registration_anchor_x[ap]=anchor_cx-abs_x;
			}
			data->registration_anchor_pos=(ap+1)%31;
			if(have_anchor&&data->registration_rolls
					&&data->registration_anchor_n<31)
				data->registration_anchor_n++;
			if(data->registration_anchor_n==31) {
				double ay[31],ax[31];
				memcpy(ay,data->registration_anchor_y,sizeof ay);
				memcpy(ax,data->registration_anchor_x,sizeof ax);
				qsort(ay,31,sizeof *ay,com_double_asc);
				qsort(ax,31,sizeof *ax,com_double_asc);
				double anchor_y=ay[15]/16.0,anchor_x=ax[15]/16.0;
				double anchor_r=hypot(anchor_y,anchor_x);
				if(anchor_r>0.25) {
					anchor_y*=0.25/anchor_r;
					anchor_x*=0.25/anchor_r;
				}
				data->registration_ref_y+=anchor_y;
				data->registration_ref_x+=anchor_x;
				abs_y+=anchor_y;abs_x+=anchor_x;
				// Stored disagreements are against the pre-correction origin.
				// Translate them with it so the median represents residual bias.
				for(size_t a=0;a<31;a++) {
					data->registration_anchor_y[a]-=anchor_y;
					data->registration_anchor_x[a]-=anchor_x;
				}
			}
			// Keep the original keyframe while it still provides enough visible
			// two-axis information. Chaining otherwise-valid keyframes integrates
			// tiny fit errors. When overlap is genuinely exhausted, re-anchor only
			// from an unclipped independent centroid; a clipped centroid is not an
			// absolute position measurement.
			if(!overlap_ok) {
				if(have_anchor) {
					abs_y=anchor_cy;abs_x=anchor_cx;
					data->registration_rolls++;
					data->registration_ready=com_install_registration_keyframe(
						data,img,org_y,org_x,abs_y,abs_x);
				}
			}
			}
		}
	}

	data->n_frames++;
	if (LIKELY(beam)) {
		com_auto_accept_peak(data,peak);
		data->last_y = -1.0 + 2*abs_y/(max_y - 1);
		data->last_x = -1.0 + 2*abs_x/(max_x - 1);
		// recentre the window for the next frame
		// Registration is expressed against a fixed keyframe/window. Moving
		// that window would change coordinates underneath the fit.
		if (!data->registration || data->acquiring) {
			data->win_y = (size_t)(abs_y + 0.5);
			data->win_x = (size_t)(abs_x + 0.5);
		}
		if (UNLIKELY(!data->had_beam)) {
			if (com_log_ok(data))
				log_info("center_of_mass: beam back after %zu "
					"frames (peak %u, dimmest row %.2f of "
					"typical); resuming updates", data->lost,
					peak, cut.mid > 0.0 ? cut.lo/cut.mid : 1.0
				);
			data->had_beam = true;
		}
		data->lost = 0;
		state->header.status &= (aylp_status)~AYLP_BEAM_LOST;
	} else {
		// No beam: hold the last good output rather than emitting a
		// centroid we do not believe. Reporting the noise centroid would
		// feed the loop a fictitious error; reporting (0,0) reads
		// downstream as "perfectly centred" and would let the integrator
		// park and then lurch when the beam reappears. Holding is the
		// only option that leaves the loop where the beam last was.
		// If the beam stays gone the window is stranded and can never
		// find it again, so fall back to a full-image re-acquire.
		data->n_held++;
		if (UNLIKELY(data->had_beam)) {
			data->n_episodes++;
			if (com_log_ok(data)) {
				if (!strcmp(why, "sum"))
					log_warn("center_of_mass: no beam in frame "
						"(no pixel in the window is above "
						"threshold %u; window peak %u); "
						"holding centroid (%.4f,%.4f)",
						data->threshold, peak,
						data->last_y, data->last_x
					);
				else if (!strcmp(why, "peak"))
					log_warn("center_of_mass: no beam in frame "
						"(window peak %u < min_peak %u); "
						"holding centroid (%.4f,%.4f)",
						peak, data->min_peak,
						data->last_y, data->last_x
					);
				else if (!strcmp(why, "registration"))
					log_warn("center_of_mass: image registration was "
						"ill-conditioned; holding centroid (%.4f,%.4f)",
						data->last_y, data->last_x);
				else
					log_warn("center_of_mass: beam cut across "
						"rows (dimmest of %zu rows is %.2f "
						"of typical, under ref_cut %G; "
						"peak %u is fine); holding centroid "
						"(%.4f,%.4f)", cut.n_sig,
						cut.mid > 0.0 ? cut.lo/cut.mid : 0.0,
						data->ref_cut, peak,
						data->last_y, data->last_x
					);
			}
			data->had_beam = false;
		}
		if (++data->lost >= data->reacquire_after) {
			state->header.status |= AYLP_BEAM_LOST;
			if (data->lost == data->reacquire_after) {
				log_warn("center_of_mass: no signal in window "
					"for %zu frames; re-acquiring",
					data->lost
				);
				if(data->auto_reacquire)
					data->reacquire_after=fmin(100000,
						2*data->reacquire_after);
				// re-enter the wide phase: whatever stranded the
				// window will likely strand it again if we drop
				// straight back to the narrow one
				data->acquiring = data->acquire_seconds > 0.0;
				data->acquire_t0 = com_monotonic_s();
				// the reference describes a beam at the old
				// window position, which is exactly what we
				// have just stopped believing in
				com_reset_ref(data);
				data->registration_ready=false;
			}
			acquire_window(data, img);
		}
	}
	// end of the acquisition phase: narrow onto wherever the centroid settled
	if (UNLIKELY(data->acquiring)) {
		double elapsed = com_monotonic_s() - data->acquire_t0;
		if (elapsed >= data->acquire_seconds) {
			data->acquiring = false;
			log_info("center_of_mass: acquisition done after %.2f s; "
				"narrowing window from %zu by %zu to %zu by %zu, "
				"centred on (%zu,%zu)", elapsed, win_h, win_w,
				data->region_height, data->region_width,
				data->win_y, data->win_x
			);
		}
	}
	data->com->data[0] = data->last_y;
	data->com->data[1] = data->last_x;

	// zero-copy update of pipeline state
	state->vector = data->com;
	// Tell the rest of the pipeline whether this centroid is a fresh
	// measurement or last frame's held over. Everything downstream sees an
	// ordinary vector either way -- holding is invisible in the data, which
	// is the point of holding -- so without this flag an integrator has no
	// way to tell that it is winding on the same error over and over.
	if (LIKELY(beam))
		state->header.status &= ~(aylp_status)AYLP_NO_SIGNAL;
	else
		state->header.status |= AYLP_NO_SIGNAL;
	// housekeeping on the header
	state->header.type = self->type_out;
	state->header.units = self->units_out;
	state->header.log_dim.y = 2;
	state->header.log_dim.x = 1;
	return 0;
}


int center_of_mass_proc_threaded(
	struct aylp_device *self, struct aylp_state *state
)
{
	struct aylp_center_of_mass_data *data = self->device_data;
	size_t max_y = state->matrix_uchar->size1;
	size_t max_x = state->matrix_uchar->size2;
	size_t y_subap_count = max_y / data->region_height;
	size_t x_subap_count = max_x / data->region_width;
	// TODO: we assume one task per subaperture, but this might not be ideal
	// for all subaperture counts and sizes. It may be that the overhead
	// associated with processing one task can be comparable to the cost of
	// calculating on one subaperture, and we should be creating tasks that
	// process more than one subaperture at a time.
	size_t n_tasks = y_subap_count * x_subap_count;
	if (!n_tasks) {
		log_error("Refusing to process zero subapertures; "
			"region size is %zu by %zu but image is %zu by %zu",
			data->region_height, data->region_width, max_y, max_x
		);
		return -1;
	}
	// It's unfortunately quite ugly that we malloc here, but it's the
	// simplest fast solution I can think of to the issue of not knowing the
	// size of state->matrix_uchar when init() is run. Remember that `data`
	// is `calloc`ed so we are guaranteed to malloc on first proc.
	if (data->n_tasks < n_tasks) {
		// we *could* realloc instead of this free/malloc combo, but I
		// expect free/malloc to be faster since we don't care about
		// keeping the old memory (https://stackoverflow.com/a/39562813)
		xfree(data->tasks);
		data->tasks = xmalloc(n_tasks * sizeof(struct aylp_task));
		data->n_tasks = n_tasks;
		// allocate the com_src inputs for the tasks
		xfree(data->com_srcs);
		data->com_srcs = xmalloc(n_tasks * sizeof(struct aylp_com_src));
	}
	// allocate the com vector if needed
	if (!data->com || data->com->size < n_tasks*2) {
		xfree_type(gsl_vector, data->com);
		data->com = xmalloc_type(gsl_vector, n_tasks*2);
	}

	// start assigning tasks
	size_t t = 0;
	for (size_t i=0; i < y_subap_count; i++) {
		for (size_t j=0; j < x_subap_count; j++) {
			// set source data
			data->com_srcs[t].mat = gsl_matrix_uchar_submatrix(
				state->matrix_uchar,
				i * data->region_height,
				j * data->region_width,
				data->region_height,
				data->region_width
			).matrix;
			data->com_srcs[t].threshold = data->threshold;
			// set task
			data->tasks[t] = (struct aylp_task){
				.func = (void(*)(void*,void*))com_mat_uchar,
				.src = &data->com_srcs[t],
				.dst = (void *)(data->com->data+2*t),
				.next_task = 0
			};
			task_enqueue(&data->queue, &data->tasks[t]);
			t += 1;
		}
	}
	// wait for threads to finish
	while (data->queue.tasks_processing) {
		sched_yield();	// (is there a better function to call?)
	}
	// zero-copy update of pipeline state
	state->vector = data->com;
	// housekeeping on the header
	state->header.type = self->type_out;
	state->header.units = self->units_out;
	state->header.log_dim.y = data->com->size;
	state->header.log_dim.x = 1;
	return 0;
}


int center_of_mass_fini_threaded(struct aylp_device *self)
{
	struct aylp_center_of_mass_data *data = self->device_data;
	shut_queue(&data->queue);
	for (size_t t = 0; t < data->thread_count; t++) {
		pthread_join(data->threads[t], 0);
	}
	xfree(data->threads);
	xfree(data->tasks);
	xfree(data->com_srcs);
	xfree(self->device_data);
	return 0;
}
