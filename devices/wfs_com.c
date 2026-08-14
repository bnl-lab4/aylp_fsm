#include <math.h>
#include <string.h>
#include <time.h>

#include "anyloop.h"
#include "logging.h"
#include "wfs_com.h"
#include "xalloc.h"


static double wfs_monotonic_s(void)
{
	struct timespec t;
	clock_gettime(CLOCK_MONOTONIC, &t);
	return t.tv_sec + 1e-9 * t.tv_nsec;
}


// index into data->ref / data->ref_set for subaperture (r,c)
static inline size_t subap_idx(struct aylp_wfs_com_data *data,
	size_t r, size_t c
) {
	return r * data->subap_cols + c;
}


// Place the window centre on the brightest pixel of the whole image. Used on
// first acquisition and after AYLP_BEAM_LOST recovery. Like center_of_mass,
// this locks onto whatever is brightest -- pass init_y/init_x if a stray
// reflection can outpeak the real beam.
static void acquire_window(
	struct aylp_wfs_com_data *data, gsl_matrix_uchar *img
) {
	unsigned char best = 0;
	size_t by = img->size1 / 2;
	size_t bx = img->size2 / 2;
	for (size_t i = 0; i < img->size1; i++) {
		for (size_t j = 0; j < img->size2; j++) {
			unsigned char v = img->data[i*img->tda + j];
			if (v > best) { best = v; by = i; bx = j; }
		}
	}
	data->win_y = by;
	data->win_x = bx;
}


// Slide the window centre so the win_h by win_w extended window (subaperture
// grid plus the search margin on every side) lies fully inside the image.
static void clamp_window(struct aylp_wfs_com_data *data,
	size_t max_y, size_t max_x
) {
	size_t half_y = data->win_h / 2;
	size_t half_x = data->win_w / 2;
	if (data->win_y < half_y) data->win_y = half_y;
	if (data->win_x < half_x) data->win_x = half_x;
	if (data->win_y > max_y - data->win_h + half_y)
		data->win_y = max_y - data->win_h + half_y;
	if (data->win_x > max_x - data->win_w + half_x)
		data->win_x = max_x - data->win_w + half_x;
}


// (Re)acquire: place the window, reset every subaperture's reference
// template (the old ones describe pixel content at the old location and are
// meaningless at the new one), and flag the coming frame as a warm-up frame
// so seeding templates isn't mistaken for a lost signal.
static void reacquire(struct aylp_wfs_com_data *data, gsl_matrix_uchar *img)
{
	if (data->init_y >= 0) {
		data->win_y = (size_t)data->init_y;
		data->win_x = (size_t)data->init_x;
	} else {
		acquire_window(data, img);
	}
	clamp_window(data, img->size1, img->size2);
	// Anchor the correlation tracker's absolute output to a LOCAL centroid
	// around the brightest pixel in the acquisition window. A centroid over
	// the whole window would reintroduce center_of_mass's failure mode: a
	// distant reflection can bias the absolute anchor even though the matched
	// filters subsequently ignore it. Correlation gives displacement relative
	// to the first template, so an accurate, robust initial anchor matters.
	size_t org_y = data->win_y - data->win_h/2;
	size_t org_x = data->win_x - data->win_w/2;
	size_t peak_y = data->win_y, peak_x = data->win_x;
	unsigned char peak = 0;
	for (size_t i = 0; i < data->win_h; i++) {
		for (size_t j = 0; j < data->win_w; j++) {
			unsigned char v = img->data[(org_y+i)*img->tda + org_x+j];
			if (v > peak) { peak = v; peak_y = org_y+i; peak_x = org_x+j; }
		}
	}
	size_t radius_y = data->subap_h/2 ? data->subap_h/2 : 1;
	size_t radius_x = data->subap_w/2 ? data->subap_w/2 : 1;
	size_t y0 = peak_y > radius_y ? peak_y-radius_y : org_y;
	size_t x0 = peak_x > radius_x ? peak_x-radius_x : org_x;
	if (y0 < org_y) y0 = org_y;
	if (x0 < org_x) x0 = org_x;
	size_t y1 = peak_y+radius_y+1;
	size_t x1 = peak_x+radius_x+1;
	if (y1 > org_y+data->win_h) y1 = org_y+data->win_h;
	if (x1 > org_x+data->win_w) x1 = org_x+data->win_w;
	double sy = 0.0, sx = 0.0, flux = 0.0;
	for (size_t i = y0; i < y1; i++) {
		for (size_t j = x0; j < x1; j++) {
			unsigned char raw = img->data[i*img->tda + j];
			double v = raw > data->threshold ? raw-data->threshold : 0.0;
			sy += i*v;
			sx += j*v;
			flux += v;
		}
	}
	memset(data->ref_set, 0, data->n_subaps * sizeof(bool));
	// Explicit init coordinates are a known steering reference and remain the
	// authoritative anchor. When acquisition itself chose the brightest pixel,
	// refine that integer location with the local centroid computed above.
	// Note the deliberate consequence of the fixed-anchor branch: the device
	// then reports the beam as sitting exactly at (init_y, init_x) on the
	// acquisition frame regardless of where inside the window it really is, so
	// a static pointing offset smaller than the search margin is absorbed into
	// the templates rather than reported. That is the point -- the configured
	// pixel, not a one-shot first moment, defines zero error -- but it does
	// mean init_y/init_x must be the intended steering target, not merely a
	// hint about where to look.
	bool fixed_anchor = data->init_y >= 0;
	data->abs_y = !fixed_anchor && flux > 0.0 ? sy/flux : (double)data->win_y;
	data->abs_x = !fixed_anchor && flux > 0.0 ? sx/flux : (double)data->win_x;
	data->ref_off_y = data->abs_y - (double)data->win_y;
	data->ref_off_x = data->abs_x - (double)data->win_x;
	// Publish the anchor as the output only on the FIRST acquisition, when
	// there is no previous measurement to hold. Every later call gets here
	// from the lost-signal branch, which has already asserted
	// AYLP_FRAME_REJECTED and promised downstream that the last valid output
	// is being held; overwriting it with a fresh acquisition guess emits a
	// position step on the very frame the device declared untrustworthy.
	if (!data->acquired) {
		data->last_y = -1.0 + 2.0*data->abs_y/((double)img->size1 - 1.0);
		data->last_x = -1.0 + 2.0*data->abs_x/((double)img->size2 - 1.0);
	}
	data->acquired = true;
	data->just_acquired = true;
	data->t_acquire = wfs_monotonic_s();
	data->lost = 0;
	log_info("wfs_com: acquired window at (%zu,%zu)", data->win_y, data->win_x);
}


int wfs_com_init(struct aylp_device *self)
{
	self->device_data = xcalloc(1, sizeof(struct aylp_wfs_com_data));
	struct aylp_wfs_com_data *data = self->device_data;
	self->fini = &wfs_com_fini;
	bool rolling_shutter_seen = false;
	bool rolling_shutter_value = false;

	// defaults
	data->threshold = 0;
	data->search_radius = 2;
	data->reject_edge_matches = false;
	data->ref_beta = 0.02;
	data->ref_beta_init = 0.0;	// <= 0 disables the schedule
	data->ref_beta_tau = 2.0;
	data->min_confidence = 0.5;
	data->flux_floor = 1.0;
	data->min_valid_subaps = 0;	// 0 => derive from n_subaps at the end
	data->row_time = 0.0;		// disabled by default
	data->rolling_shutter = false;
	data->max_row_shear = 0.0;	// disabled by default
	data->init_y = -1;
	data->init_x = -1;
	data->reacquire_after = 10;

	if (!self->params) {
		log_error("No params object found.");
		return -1;
	}
	json_object_object_foreach(self->params, key, val) {
		if (key[0] == '_') {
		} else if (!strcmp(key, "subap_height")) {
			data->subap_h = json_object_get_uint64(val);
		} else if (!strcmp(key, "subap_width")) {
			data->subap_w = json_object_get_uint64(val);
		} else if (!strcmp(key, "subap_rows")) {
			data->subap_rows = json_object_get_uint64(val);
		} else if (!strcmp(key, "subap_cols")) {
			data->subap_cols = json_object_get_uint64(val);
		} else if (!strcmp(key, "threshold")) {
			int threshold = json_object_get_int(val);
			if (threshold < 0 || threshold > 255) {
				log_error("wfs_com: threshold must be in [0,255]");
				return -1;
			}
			data->threshold = (unsigned char)threshold;
		} else if (!strcmp(key, "search_radius")) {
			data->search_radius = json_object_get_uint64(val);
		} else if (!strcmp(key, "reject_edge_matches")) {
			data->reject_edge_matches = json_object_get_boolean(val);
		} else if (!strcmp(key, "ref_beta")) {
			data->ref_beta = json_object_get_double(val);
		} else if (!strcmp(key, "ref_beta_init")) {
			data->ref_beta_init = json_object_get_double(val);
		} else if (!strcmp(key, "ref_beta_tau")) {
			data->ref_beta_tau = json_object_get_double(val);
		} else if (!strcmp(key, "min_confidence")) {
			data->min_confidence = json_object_get_double(val);
		} else if (!strcmp(key, "flux_floor")) {
			data->flux_floor = json_object_get_double(val);
		} else if (!strcmp(key, "min_valid_subaps")) {
			data->min_valid_subaps = json_object_get_uint64(val);
		} else if (!strcmp(key, "row_time")) {
			data->row_time = json_object_get_double(val);
		} else if (!strcmp(key, "rolling_shutter")) {
			rolling_shutter_seen = true;
			rolling_shutter_value = json_object_get_boolean(val);
		} else if (!strcmp(key, "max_row_shear")) {
			data->max_row_shear = json_object_get_double(val);
		} else if (!strcmp(key, "init_y")) {
			data->init_y = json_object_get_int64(val);
		} else if (!strcmp(key, "init_x")) {
			data->init_x = json_object_get_int64(val);
		} else if (!strcmp(key, "reacquire_after")) {
			data->reacquire_after = json_object_get_uint64(val);
		} else {
			log_warn("Unknown parameter \"%s\"", key);
		}
	}
	data->rolling_shutter = rolling_shutter_seen
		? rolling_shutter_value : data->row_time > 0.0;

	if (!data->subap_h || !data->subap_w
			|| !data->subap_rows || !data->subap_cols) {
		log_error("wfs_com: subap_height, subap_width, subap_rows and "
			"subap_cols must all be nonzero");
		return -1;
	}
	if ((data->init_y < 0) != (data->init_x < 0)) {
		log_error("Provide both init_y and init_x, or neither");
		return -1;
	}
	if (!data->search_radius
			|| data->search_radius > AYLP_WFS_COM_MAX_SEARCH_RADIUS) {
		log_error("wfs_com: search_radius must be in [1, %d]",
			AYLP_WFS_COM_MAX_SEARCH_RADIUS);
		return -1;
	}
	if (!(data->ref_beta > 0.0 && data->ref_beta <= 1.0)) {
		log_error("wfs_com: ref_beta must be in (0, 1]");
		return -1;
	}
	if (!isfinite(data->ref_beta_init)) {
		log_error("wfs_com: ref_beta_init must be finite");
		return -1;
	}
	if (data->ref_beta_init > 0.0
			&& (!isfinite(data->ref_beta_tau)
				|| data->ref_beta_tau <= 0.0)) {
		log_error("wfs_com: ref_beta_tau must be positive when "
			"ref_beta_init > 0");
		return -1;
	}
	if (data->ref_beta_init > 1.0) {
		log_error("wfs_com: ref_beta_init must be <= 1");
		return -1;
	}
	if (!isfinite(data->min_confidence) || data->min_confidence < 0.0
			|| data->min_confidence > 1.0) {
		log_error("wfs_com: min_confidence must be in [0,1]");
		return -1;
	}
	if (!isfinite(data->flux_floor) || data->flux_floor < 0.0) {
		log_error("wfs_com: flux_floor must be finite and non-negative");
		return -1;
	}
	if (!isfinite(data->row_time)) {
		log_error("wfs_com: row_time must be finite");
		return -1;
	}
	if (!isfinite(data->max_row_shear) || data->max_row_shear < 0.0) {
		log_error("wfs_com: max_row_shear must be finite and non-negative");
		return -1;
	}
	if (!data->reacquire_after) {
		log_error("wfs_com: reacquire_after must be nonzero");
		return -1;
	}

	size_t margin = 2*data->search_radius;
	if (data->subap_rows > SIZE_MAX/data->subap_cols
			|| data->subap_rows > (SIZE_MAX-margin)/data->subap_h
			|| data->subap_cols > (SIZE_MAX-margin)/data->subap_w) {
		log_error("wfs_com: subaperture geometry overflows size_t");
		return -1;
	}
	data->n_subaps = data->subap_rows * data->subap_cols;
	if (data->n_subaps > SIZE_MAX/data->subap_h
			|| data->n_subaps*data->subap_h > SIZE_MAX/data->subap_w) {
		log_error("wfs_com: reference allocation size overflows size_t");
		return -1;
	}
	data->ext_h = data->subap_h + 2*data->search_radius;
	data->ext_w = data->subap_w + 2*data->search_radius;
	data->win_h = data->subap_rows*data->subap_h + 2*data->search_radius;
	data->win_w = data->subap_cols*data->subap_w + 2*data->search_radius;
	if (data->ext_h > SIZE_MAX/data->ext_w) {
		log_error("wfs_com: extended scratch size overflows size_t");
		return -1;
	}
	if (!data->min_valid_subaps)
		data->min_valid_subaps = (data->n_subaps + 3) / 4;	// ~25%
	if (data->min_valid_subaps > data->n_subaps) {
		log_error("wfs_com: min_valid_subaps (%zu) exceeds the "
			"subaperture count (%zu)", data->min_valid_subaps,
			data->n_subaps);
		return -1;
	}

	data->ref = xcalloc(data->n_subaps * data->subap_h * data->subap_w,
		sizeof(double));
	data->ref_backup = xcalloc(
		data->n_subaps * data->subap_h * data->subap_w, sizeof(double));
	data->ref_set = xcalloc(data->n_subaps, sizeof(bool));
	data->ref_set_backup = xcalloc(data->n_subaps, sizeof(bool));
	data->matched = xcalloc(data->n_subaps, sizeof(bool));
	data->sub_dy = xcalloc(data->n_subaps, sizeof(double));
	data->sub_dx = xcalloc(data->n_subaps, sizeof(double));
	data->ext = xcalloc(data->ext_h * data->ext_w, sizeof(double));
	data->row_sum_w = xcalloc(data->subap_rows, sizeof(double));
	data->row_sum_wt = xcalloc(data->subap_rows, sizeof(double));
	data->row_sum_wdy = xcalloc(data->subap_rows, sizeof(double));
	data->row_sum_wdx = xcalloc(data->subap_rows, sizeof(double));
	data->com = xmalloc_type(gsl_vector, 2);

	self->type_in = AYLP_T_MATRIX_UCHAR;
	self->units_in = AYLP_U_ANY;
	self->type_out = AYLP_T_VECTOR;
	self->units_out = AYLP_U_MINMAX;
	self->proc = &wfs_com_proc;

	log_info("wfs_com: %zux%zu subapertures of %zux%zu px, search radius "
		"%zu, window %zux%zu", data->subap_rows, data->subap_cols,
		data->subap_h, data->subap_w, data->search_radius,
		data->win_h, data->win_w);
	return 0;
}


int wfs_com_proc(struct aylp_device *self, struct aylp_state *state)
{
	struct aylp_wfs_com_data *data = self->device_data;
	gsl_matrix_uchar *img = state->matrix_uchar;
	if (UNLIKELY(!img || !img->data || img->size1 < 2 || img->size2 < 2
			|| img->tda < img->size2)) {
		log_error("wfs_com: invalid image matrix or dimensions");
		return -1;
	}
	size_t max_y = img->size1, max_x = img->size2;

	if (UNLIKELY(data->win_h > max_y || data->win_w > max_x)) {
		log_error("wfs_com: window is %zux%zu but image is only "
			"%zux%zu", data->win_h, data->win_w, max_y, max_x);
		return -1;
	}
	if (UNLIKELY(!data->acquired)) reacquire(data, img);

	// Template SEEDS written while matching are speculative until the
	// frame-level validity and rolling-shear gates have seen every
	// subaperture; a rejected frame must not leave a subaperture marked
	// initialized from content it just judged untrustworthy. (The EWMA
	// re-registration needs no rollback: it runs in the second pass, which a
	// rejected frame never reaches.) The acquisition warm-up must still seed.
	size_t ref_elems = data->n_subaps * data->subap_h * data->subap_w;
	bool refs_speculative = !data->just_acquired;
	if (refs_speculative) {
		memcpy(data->ref_backup, data->ref, ref_elems * sizeof(double));
		memcpy(data->ref_set_backup, data->ref_set,
			data->n_subaps * sizeof(bool));
	}

	// NLMS-style step schedule (see wfs_com.h): fast right after
	// (re)acquisition while templates are still settling, relaxing to
	// the steady-state ref_beta with time constant ref_beta_tau.
	double ref_beta_cur = data->ref_beta;
	if (data->ref_beta_init > 0.0) {
		double dt = wfs_monotonic_s() - data->t_acquire;
		ref_beta_cur = data->ref_beta
			+ (data->ref_beta_init - data->ref_beta)
				* exp(-dt / data->ref_beta_tau);
	}

	size_t org_y = data->win_y - data->win_h/2;
	size_t org_x = data->win_x - data->win_w/2;
	size_t R = data->search_radius;
	size_t grid = 2*R + 1;
	size_t sh = data->subap_h, sw = data->subap_w;

	double sum_w = 0.0, sum_wdy = 0.0, sum_wdx = 0.0;
	size_t valid_count = 0;
	memset(data->matched, 0, data->n_subaps * sizeof(bool));
	memset(data->row_sum_w, 0, data->subap_rows * sizeof(double));
	memset(data->row_sum_wt, 0, data->subap_rows * sizeof(double));
	memset(data->row_sum_wdy, 0, data->subap_rows * sizeof(double));
	memset(data->row_sum_wdx, 0, data->subap_rows * sizeof(double));

	for (size_t r = 0; r < data->subap_rows; r++) {
		for (size_t c = 0; c < data->subap_cols; c++) {
			size_t idx = subap_idx(data, r, c);
			double *ref = data->ref + idx*sh*sw;
			size_t ext_top = org_y + r*sh;
			size_t ext_left = org_x + c*sw;

			// extract + threshold the extended (margin-padded) block
			for (size_t i = 0; i < data->ext_h; i++) {
				for (size_t j = 0; j < data->ext_w; j++) {
					unsigned char raw = img->data[
						(ext_top+i)*img->tda + ext_left+j
					];
					unsigned char el = raw > data->threshold
						? raw - data->threshold : 0;
					data->ext[i*data->ext_w + j] = (double)el;
				}
			}
			// flux of the nominal (unshifted) core cell, used for
			// the flux gate and as an SNR weight
			double flux = 0.0;
			for (size_t i = 0; i < sh; i++)
				for (size_t j = 0; j < sw; j++)
					flux += data->ext[(R+i)*data->ext_w + (R+j)];

			if (UNLIKELY(!data->ref_set[idx])) {
				// seed the reference template from the core
				// cell; no shift estimate is possible yet
				if (flux >= data->flux_floor) {
					for (size_t i = 0; i < sh; i++)
						for (size_t j = 0; j < sw; j++)
							ref[i*sw+j] = data->ext[
								(R+i)*data->ext_w + (R+j)
							];
					data->ref_set[idx] = true;
				}
				continue;
			}
			// exhaustive integer-pixel normalized cross-correlation
			double ncc[AYLP_WFS_COM_MAX_GRID][AYLP_WFS_COM_MAX_GRID];
			double energy_r = 0.0;
			for (size_t i = 0; i < sh*sw; i++)
				energy_r += ref[i]*ref[i];
			double best = -2.0;
			size_t best_sy = R, best_sx = R;
			for (size_t sy = 0; sy < grid; sy++) {
				for (size_t sx = 0; sx < grid; sx++) {
					double corr = 0.0, energy_c = 0.0;
					for (size_t i = 0; i < sh; i++) {
						for (size_t j = 0; j < sw; j++) {
							double v = data->ext[
								(sy+i)*data->ext_w + sx+j
							];
							corr += v * ref[i*sw+j];
							energy_c += v*v;
						}
					}
					double n = corr
						/ sqrt(energy_c*energy_r + 1e-9);
					ncc[sy][sx] = n;
					if (n > best) {
						best = n; best_sy = sy; best_sx = sx;
					}
				}
			}
			double confidence = best < 0.0 ? 0.0 : best;
			if (confidence < data->min_confidence) continue;
			// Gate and weight the patch that ACTUALLY matched. Using
			// the nominal core's flux here makes a valid translated
			// feature appear to go dark merely because it crossed a
			// cell boundary.
			double matched_flux = 0.0;
			for (size_t i = 0; i < sh; i++)
				for (size_t j = 0; j < sw; j++)
					matched_flux += data->ext[
						(best_sy+i)*data->ext_w + best_sx+j];
			if (matched_flux < data->flux_floor) continue;
			// A boundary peak says only "the displacement is at least
			// this large". Treating it as exact lets a distant
			// artifact walk the window by R pixels on every bad frame.
			// Steering configs can reject these censored estimates
			// until an interior match returns.
			if (data->reject_edge_matches
					&& (best_sy == 0 || best_sy == grid-1
						|| best_sx == 0 || best_sx == grid-1))
				continue;

			// Full 2D (non-separable) quadratic peak fit over the 3x3
			// neighborhood around the integer peak, skipped (left at
			// the integer estimate) if the peak sits on the edge of
			// the search range -- widen search_radius if that happens
			// often. Fitting y and x independently (two separate 1D
			// parabolic fits) ignores the cross (row*col) term of the
			// correlation surface, which is a well-documented source
			// of systematic "peak locking" bias in correlation
			// trackers. This fits f(i,j) = A i^2 + B j^2 + C i j +
			// D i + E j + F by ordinary least squares over the 9
			// samples (closed form, since the 3x3 regular grid's
			// design matrix is fixed and its Gram matrix inverts
			// analytically) and solves for the true vertex of that
			// surface, coupling term included.
			double dy_frac = 0.0, dx_frac = 0.0;
			if (best_sy > 0 && best_sy < grid-1
					&& best_sx > 0 && best_sx < grid-1) {
				double fm1m1 = ncc[best_sy-1][best_sx-1];
				double fm10  = ncc[best_sy-1][best_sx];
				double fm1p1 = ncc[best_sy-1][best_sx+1];
				double f0m1  = ncc[best_sy][best_sx-1];
				double f0p1  = ncc[best_sy][best_sx+1];
				double fp1m1 = ncc[best_sy+1][best_sx-1];
				double fp10  = ncc[best_sy+1][best_sx];
				double fp1p1 = ncc[best_sy+1][best_sx+1];

				double row_m1 = fm1m1+fm10+fm1p1;
				double row_p1 = fp1m1+fp10+fp1p1;
				double col_m1 = fm1m1+f0m1+fp1m1;
				double col_p1 = fm1p1+f0p1+fp1p1;
				double S0 = row_m1+row_p1
					+ (f0m1+ncc[best_sy][best_sx]+f0p1);
				double Sii = row_m1+row_p1;
				double Sjj = col_m1+col_p1;
				double Si = row_p1-row_m1;
				double Sj = col_p1-col_m1;
				double Sij = fp1p1-fp1m1-fm1p1+fm1m1;

				double F = (5*S0 - 3*Sii - 3*Sjj)/9.0;
				double A = (S0 - Sjj - 3*F)/2.0;
				double B = (S0 - Sii - 3*F)/2.0;
				double C = Sij/4.0;
				double D = Si/6.0;
				double E = Sj/6.0;

				// Guard: a genuine correlation peak fits a concave
				// (negative-definite) surface. A < 0 with a
				// positive determinant confirms that (and, since
				// det=4AB-C^2, implies B < 0 too); anything else is
				// a degenerate or saddle-shaped local fit, which
				// happens near very low-confidence/noisy matches --
				// fall back to the integer estimate rather than
				// solving a nonsense vertex.
				double det = 4*A*B - C*C;
				if (A < 0.0 && det > 1e-12) {
					dy_frac = (-2*B*D + C*E) / det;
					dx_frac = (-2*A*E + C*D) / det;
					if (dy_frac < -1.0) dy_frac = -1.0;
					if (dy_frac > 1.0) dy_frac = 1.0;
					if (dx_frac < -1.0) dx_frac = -1.0;
					if (dx_frac > 1.0) dx_frac = 1.0;
				}
			}

			double dy = ((double)best_sy - (double)R) + dy_frac;
			double dx = ((double)best_sx - (double)R) + dx_frac;
			double weight = confidence * matched_flux;
			sum_w += weight;
			sum_wdy += weight * dy;
			sum_wdx += weight * dx;
			data->row_sum_w[r] += weight;
			data->row_sum_wt[r] += weight
				* ((double)r*sh + (double)sh/2.0 + dy);
			data->row_sum_wdy[r] += weight * dy;
			data->row_sum_wdx[r] += weight * dx;
			valid_count++;
			// Flag for the second pass. The template is NOT updated
			// here: re-registering it now, with this subaperture's
			// own sub-pixel estimate, closes a private feedback loop
			// (estimate -> template position -> next estimate) that
			// is unstable wherever the correlation surface is
			// ill-conditioned. A subaperture seeing only a smooth
			// monotone gradient -- the outer cells of any beam that
			// doesn't fill the grid -- has an NCC that is ~1 at
			// EVERY offset, so its peak location is noise, yet it
			// sails through min_confidence. Feeding that noise back
			// into its own template ratchets the template outward a
			// little more every frame until the peak pins to the
			// search boundary. Measured on a byte-identical static
			// frame with contrib/legacy/legacy_steering_profiles/steering_par_fsp_wfs_com.json, this
			// fabricated ~0.013 px/frame of drift, then failed every
			// subsequent frame's validity gate, then asserted
			// AYLP_BEAM_LOST and re-acquired roughly every 34 frames
			// -- on a perfectly good, motionless beam. The whole
			// premise of the device is that real motion is COMMON
			// MODE across subapertures, so the consensus shift is
			// the only registration signal that is averaged enough
			// to be stable; see the second pass below.
			data->matched[idx] = true;
			data->sub_dy[idx] = dy;
			data->sub_dx[idx] = dx;
		}
	}

	// Per-row-group weighted-mean shift, computed regardless of row_time:
	// used both for the shear-rejection gate below and, when enabled,
	// as the rolling-shutter regression's data points. If row-groups
	// disagree by more than max_row_shear pixels in y or x, subapertures
	// on different rows genuinely saw the target at different positions
	// within this one frame's readout -- exactly the rolling-shutter
	// signature -- and the whole frame is rejected below rather than
	// blended into one (distorted) global shift.
	double row_y_min = 0.0, row_y_max = 0.0, row_x_min = 0.0, row_x_max = 0.0;
	bool shear_rejected_frame = false;
	if (data->max_row_shear > 0.0) {
		size_t n_rows_seen = 0;
		for (size_t r = 0; r < data->subap_rows; r++) {
			double w = data->row_sum_w[r];
			if (w <= 0.0) continue;
			double y_r = data->row_sum_wdy[r] / w;
			double x_r = data->row_sum_wdx[r] / w;
			if (!n_rows_seen || y_r < row_y_min) row_y_min = y_r;
			if (!n_rows_seen || y_r > row_y_max) row_y_max = y_r;
			if (!n_rows_seen || x_r < row_x_min) row_x_min = x_r;
			if (!n_rows_seen || x_r > row_x_max) row_x_max = x_r;
			n_rows_seen++;
		}
		shear_rejected_frame = n_rows_seen >= 2
			&& (row_y_max - row_y_min > data->max_row_shear
				|| row_x_max - row_x_min > data->max_row_shear);
	}

	if (UNLIKELY(data->just_acquired)) {
		// warm-up frame: every subaperture just seeded its template
		// and has nothing to report yet. Not a lost signal.
		data->just_acquired = false;
		state->header.status |= AYLP_FRAME_REJECTED;
	} else if (valid_count < data->min_valid_subaps || sum_w <= 0.0) {
		// Too few confident subapertures this frame -- e.g. a chop-off
		// frame, a transient occlusion, or genuine signal loss. Hold
		// the last valid output rather than reporting a fake position,
		// same rationale as anyloop:center_of_mass.
		state->header.status |= AYLP_FRAME_REJECTED;
		if (refs_speculative) {
			memcpy(data->ref, data->ref_backup,
				ref_elems * sizeof(double));
			memcpy(data->ref_set, data->ref_set_backup,
				data->n_subaps * sizeof(bool));
		}
		if (++data->lost >= data->reacquire_after) {
			state->header.status |= AYLP_BEAM_LOST;
			if (data->lost == data->reacquire_after) {
				// Report how many subapertures were even SEEDED,
				// not just how many matched: a min_valid_subaps
				// larger than the number of cells the beam ever
				// illuminates can never be met, and without this
				// the symptom is an unexplained continuous
				// re-acquire loop rather than an obviously
				// mistuned threshold.
				size_t seeded = 0;
				for (size_t i = 0; i < data->n_subaps; i++)
					if (data->ref_set[i]) seeded++;
				log_warn("wfs_com: only %zu of %zu required "
					"subapertures matched (%zu of %zu have a "
					"seeded template; flux_floor %.3g) for "
					"%zu frames; re-acquiring", valid_count,
					data->min_valid_subaps, seeded,
					data->n_subaps, data->flux_floor,
					data->lost);
			}
			reacquire(data, img);
		}
	} else if (shear_rejected_frame) {
		// Rolling-shutter shear rejection: hold the last valid output
		// instead of reporting a shear-distorted position. Deliberately
		// NOT counted toward data->lost/reacquire_after -- see
		// max_row_shear's comment in wfs_com.h: the geometry hasn't
		// changed, only this frame's readout timing, so a sustained
		// rolling-shutter condition (e.g. continuous vibration) must not
		// force a repeated window/template reset.
		//
		// Reaching this branch means the frame PASSED the
		// min_valid_subaps gate above, i.e. the beam was demonstrably
		// present and well-matched -- so the consecutive-loss counter
		// resets here exactly as it does on a normally tracked frame.
		// Merely leaving it frozen would let low-signal frames
		// interleaved with shear-rejected ones accumulate toward
		// reacquire_after across frames where the beam was repeatedly
		// seen just fine, eventually firing a spurious reacquire.
		data->lost = 0;
		if (refs_speculative) {
			memcpy(data->ref, data->ref_backup,
				ref_elems * sizeof(double));
			memcpy(data->ref_set, data->ref_set_backup,
				data->n_subaps * sizeof(bool));
		}
		data->shear_rejected++;
		state->header.status |= AYLP_FRAME_REJECTED;
		log_debug("wfs_com: rejected frame for rolling-shutter shear "
			"(row spread y=%.2f x=%.2f px > max_row_shear %.2f px; "
			"%zu total)", row_y_max - row_y_min,
			row_x_max - row_x_min, data->max_row_shear,
			data->shear_rejected);
	} else {
		state->header.status &= (aylp_status)~AYLP_FRAME_REJECTED;
		double global_dy = sum_wdy / sum_w;
		double global_dx = sum_wdx / sum_w;
		if (data->rolling_shutter) {
			// Rolling-shutter correction: rather than blending every
			// subaperture's shift together regardless of row (which
			// mixes measurements taken at genuinely different
			// instants into one smeared number), first collapse to
			// one (weighted) shift per row group, then fit a line
			// against each row group's own capture time and report
			// the position at a fixed reference time (the window's
			// vertical centre) instead of at the blended,
			// ill-defined "average" instant.
			double sw_t = 0.0, swt = 0.0, swy = 0.0, swx = 0.0;
			size_t n_rows_seen = 0;
			for (size_t r = 0; r < data->subap_rows; r++) {
				double w = data->row_sum_w[r];
				if (w <= 0.0) continue;
				double t = data->row_sum_wt[r] / w;
				sw_t += w; swt += w*t;
				swy += data->row_sum_wdy[r];	// already weight*dy
				swx += data->row_sum_wdx[r];	// already weight*dx
				n_rows_seen++;
			}
			if (n_rows_seen >= 2 && sw_t > 0.0) {
				double t_mean = swt / sw_t;
				double y_mean = swy / sw_t;
				double x_mean = swx / sw_t;
				double s_tt = 0.0, s_ty = 0.0, s_tx = 0.0;
				double s_yy = 0.0, s_xx = 0.0;
				for (size_t r = 0; r < data->subap_rows; r++) {
					double w = data->row_sum_w[r];
					if (w <= 0.0) continue;
					double t = data->row_sum_wt[r] / w;
					double dt = t - t_mean;
					double y_r = data->row_sum_wdy[r] / w;
					double x_r = data->row_sum_wdx[r] / w;
					s_tt += w * dt * dt;
					s_ty += w * dt * (y_r - y_mean);
					s_tx += w * dt * (x_r - x_mean);
					s_yy += w * (y_r-y_mean)*(y_r-y_mean);
					s_xx += w * (x_r-x_mean)*(x_r-x_mean);
				}
				double t_ref = (double)(data->subap_rows*sh)/2.0;
				if (s_tt > 1e-12 && n_rows_seen >= 3) {
					double b_y = s_ty / s_tt;
					double b_x = s_tx / s_tt;
					// Shrink each slope by its own statistical
					// significance before extrapolating with
					// it. The lever arm (t_ref - t_mean)
					// reaches several pixels, so an
					// insignificant slope is multiplied into a
					// large bogus correction: measured on a
					// 10 px/250 Hz scene this MORE THAN
					// DOUBLED the error (0.251 -> 0.519 px
					// debiased) versus not correcting at all.
					// With only subap_rows samples, of which
					// a beam that doesn't fill the grid lights
					// 2-3, the slope is mostly noise.
					// Weighting it by b^2/(b^2 + var(b))
					// leaves a well-determined slope almost
					// untouched and collapses a
					// noise-dominated one to zero, degrading
					// gracefully to the plain weighted mean.
					// n_rows_seen >= 3 is required so there is
					// at least one residual degree of freedom
					// to estimate var(b) from at all; a
					// 2-point line is exact by construction
					// and its slope is pure noise.
					double dof = (double)n_rows_seen - 2.0;
					double sse_y = s_yy - b_y*b_y*s_tt;
					double sse_x = s_xx - b_x*b_x*s_tt;
					if (sse_y < 0.0) sse_y = 0.0;
					if (sse_x < 0.0) sse_x = 0.0;
					double vy = sse_y/(dof*s_tt);
					double vx = sse_x/(dof*s_tt);
					b_y *= b_y*b_y/(b_y*b_y + vy + 1e-18);
					b_x *= b_x*b_x/(b_x*b_x + vx + 1e-18);
					global_dy = y_mean + b_y*(t_ref - t_mean);
					global_dx = x_mean + b_x*(t_ref - t_mean);
				} else {
					global_dy = y_mean;
					global_dx = x_mean;
				}
			}
		}

		// Second pass: re-register every contributing subaperture's
		// reference template against the FRAME CONSENSUS shift, which is
		// only known now. Sampling the observation at the consensus
		// displacement and storing it at the template's nominal index
		// divides out the motion, so the template stays pinned to the
		// window lattice and keeps measuring displacement against the
		// acquisition anchor (the invariant the reconstruction below
		// relies on) while still adapting to slow changes in beam SHAPE.
		// Unlike a per-subaperture registration this cannot run away:
		// the consensus averages every contributing subaperture, so one
		// ill-conditioned cell's noisy estimate is attenuated instead of
		// being fed straight back into its own template. Runs only on
		// accepted frames, so a rejected or shear-distorted frame can
		// never teach the matched filters its own distortion.
		//
		// Register against a TRIMMED consensus, not the reported one.
		// The reported position is a flux-weighted average over every
		// valid subaperture, so a persistent asymmetric distortion (a
		// partial occlusion, a stray reflection sitting in one corner)
		// biases it -- that is unavoidable, and down-weighting is the
		// best that can be done for the measurement. Registering the
		// templates with that biased number is avoidable, though, and
		// would be much worse: it writes the bias into the references,
		// where it persists after the distortion clears. Averaging only
		// the subapertures that agree with each other keeps the
		// registration on the undistorted majority.
		double reg_dy = global_dy, reg_dx = global_dx;
		{
			double n = 0.0, tdy = 0.0, tdx = 0.0;
			for (size_t i = 0; i < data->n_subaps; i++) {
				if (!data->matched[i]) continue;
				if (fabs(data->sub_dy[i] - global_dy)
						> AYLP_WFS_COM_CONSENSUS_TOL
					|| fabs(data->sub_dx[i] - global_dx)
						> AYLP_WFS_COM_CONSENSUS_TOL)
					continue;
				n += 1.0;
				tdy += data->sub_dy[i];
				tdx += data->sub_dx[i];
			}
			// Unweighted on purpose: this is a consistency estimate,
			// and flux weighting would let one bright distorted cell
			// dominate the very quantity meant to exclude it. Falls
			// back to the untrimmed consensus if the subapertures
			// agree with nothing.
			if (n > 0.0) { reg_dy = tdy/n; reg_dx = tdx/n; }
		}
		double sy_f = (double)R + reg_dy;
		double sx_f = (double)R + reg_dx;
		// Keep the bilinear footprint (sy0..sy0+1 + sh-1) inside
		// ext_h/w. A consensus shift can exceed +-search_radius once the
		// per-subaperture sub-pixel refinements are averaged in.
		if (sy_f < 0.0) sy_f = 0.0;
		if (sy_f > (double)(2*R)) sy_f = (double)(2*R);
		if (sx_f < 0.0) sx_f = 0.0;
		if (sx_f > (double)(2*R)) sx_f = (double)(2*R);
		size_t sy0 = (size_t)sy_f, sx0 = (size_t)sx_f;
		double wy = sy_f - (double)sy0, wx = sx_f - (double)sx0;
		if (sy0 >= 2*R) { sy0 = 2*R; wy = 0.0; }
		if (sx0 >= 2*R) { sx0 = 2*R; wx = 0.0; }
		for (size_t r = 0; r < data->subap_rows; r++) {
			for (size_t c = 0; c < data->subap_cols; c++) {
				size_t idx = subap_idx(data, r, c);
				if (!data->matched[idx]) continue;
				// Freeze the template of any subaperture that
				// disagrees with the consensus. This is the gate
				// min_confidence was meant to be: a subaperture
				// that is partially occluded, or that a speckle
				// has landed on, still scores NCC ~1 against its
				// own template (a smooth patch correlates ~1 at
				// every offset), so confidence alone never
				// catches it -- but its disagreement with the
				// rigid-translation consensus does. Re-registering
				// it with a shift it plainly did not undergo is
				// what would corrupt its reference.
				if (fabs(data->sub_dy[idx] - global_dy)
						> AYLP_WFS_COM_CONSENSUS_TOL
					|| fabs(data->sub_dx[idx] - global_dx)
						> AYLP_WFS_COM_CONSENSUS_TOL)
					continue;
				double *ref = data->ref + idx*sh*sw;
				size_t ext_top = org_y + r*sh;
				size_t ext_left = org_x + c*sw;
				// re-extract: ext holds one subaperture at a time
				for (size_t i = 0; i < data->ext_h; i++) {
					for (size_t j = 0; j < data->ext_w; j++) {
						unsigned char raw = img->data[
							(ext_top+i)*img->tda
								+ ext_left+j
						];
						unsigned char el =
							raw > data->threshold
							? raw - data->threshold : 0;
						data->ext[i*data->ext_w + j] =
							(double)el;
					}
				}
				for (size_t i = 0; i < sh; i++) {
					size_t y0 = sy0 + i;
					size_t y1 = wy > 0.0 ? y0 + 1 : y0;
					for (size_t j = 0; j < sw; j++) {
						size_t x0 = sx0 + j;
						size_t x1 = wx > 0.0 ? x0 + 1 : x0;
						double v00 = data->ext[
							y0*data->ext_w + x0];
						double v01 = data->ext[
							y0*data->ext_w + x1];
						double v10 = data->ext[
							y1*data->ext_w + x0];
						double v11 = data->ext[
							y1*data->ext_w + x1];
						double v = v00*(1-wy)*(1-wx)
							+ v01*(1-wy)*wx
							+ v10*wy*(1-wx)
							+ v11*wy*wx;
						ref[i*sw+j] =
							(1.0-ref_beta_cur)*ref[i*sw+j]
							+ ref_beta_cur*v;
					}
				}
			}
		}

		// Each match is a residual relative to the template at the
		// CURRENT integer window placement, not an increment relative to
		// last frame's output. Accumulating it repeatedly creates
		// artificial drift whenever a sub-pixel residual persists. The
		// acquisition centroid supplies the fixed template-to-window
		// offset; reconstruct absolute position fresh.
		data->abs_y = (double)data->win_y + data->ref_off_y + global_dy;
		data->abs_x = (double)data->win_x + data->ref_off_x + global_dx;
		// Derive the next frame's integer window placement from the
		// sub-pixel position, while retaining the measured coordinate
		// even if the window itself must clamp near an image edge.
		// Guard against casting a negative double to size_t (which would
		// wrap around to a huge value and clamp to the WRONG edge below)
		// if a still-noisy weighted average ever pushes the position
		// outside the image before clamp_window() gets a chance to pull
		// it back in.
		if (data->abs_y < 0.0) data->abs_y = 0.0;
		if (data->abs_y > (double)(max_y - 1))
			data->abs_y = (double)(max_y - 1);
		if (data->abs_x < 0.0) data->abs_x = 0.0;
		if (data->abs_x > (double)(max_x - 1))
			data->abs_x = (double)(max_x - 1);
		data->win_y = (size_t)llround(data->abs_y);
		data->win_x = (size_t)llround(data->abs_x);
		clamp_window(data, max_y, max_x);

		data->last_y = -1.0 + 2.0*data->abs_y/((double)max_y - 1.0);
		data->last_x = -1.0 + 2.0*data->abs_x/((double)max_x - 1.0);
		data->lost = 0;
		state->header.status &= (aylp_status)~AYLP_BEAM_LOST;
	}

	data->com->data[0] = data->last_y;
	data->com->data[1] = data->last_x;
	state->vector = data->com;
	state->header.type = self->type_out;
	state->header.units = self->units_out;
	state->header.log_dim.y = 2;
	state->header.log_dim.x = 1;
	return 0;
}


int wfs_com_fini(struct aylp_device *self)
{
	struct aylp_wfs_com_data *data = self->device_data;
	if (data) {
		xfree(data->ref);
		xfree(data->ref_backup);
		xfree(data->ref_set);
		xfree(data->ref_set_backup);
		xfree(data->matched);
		xfree(data->sub_dy);
		xfree(data->sub_dx);
		xfree(data->ext);
		xfree(data->row_sum_w);
		xfree(data->row_sum_wt);
		xfree(data->row_sum_wdy);
		xfree(data->row_sum_wdx);
		if (data->com) xfree_type(gsl_vector, data->com);
	}
	xfree(self->device_data);
	return 0;
}
