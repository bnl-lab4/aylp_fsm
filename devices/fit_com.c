#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "anyloop.h"
#include "logging.h"
#include "fit_com.h"
#include "xalloc.h"

#define NP AYLP_FIT_COM_NP
// The jacobian's seven columns plus the residual, and the packed upper triangle
// of their Gram matrix -- which holds JtJ, Jtr and the cost all at once.
#define NCOL (NP+1)
#define GTRI(a,b) ((a)*NCOL - (a)*((a)+1)/2 + (b))
#define NGRAM (NCOL*(NCOL+1)/2)
// Rows of the core whitened at a time. Small enough that the scratch block
// stays in L1 for the pass that then walks it once per column.
#define JAC_BLOCK 4

// The kernels are written against a 4-wide double vector. GCC lowers that to
// whatever the target actually has, so one source compiles for baseline SSE2
// and for AVX2+FMA; target_clones emits both and an ifunc resolver picks at
// load time, so the wide path is used where it exists without the binary
// refusing to start where it does not.
#define VW 4
typedef double v4d __attribute__((vector_size(VW*sizeof(double))));
typedef double v4du __attribute__((vector_size(VW*sizeof(double)), aligned(8)));
#define VLD(p) ((v4d)*(const v4du *)(p))
#define VST(p, v) (*(v4du *)(p) = (v4du)(v))
#define VBC(s) __extension__({ double vbc_ = (s); (v4d){vbc_,vbc_,vbc_,vbc_}; })
#define VSUM(v) __extension__({ v4d vs_ = (v); vs_[0]+vs_[1]+vs_[2]+vs_[3]; })

// -DAYLP_NO_KERNEL_CLONES forces the baseline path, which is how the fallback
// gets tested on a machine that would otherwise always take the wide one.
#if defined(AYLP_NO_KERNEL_CLONES)
#define AYLP_KERNEL
#elif defined(__x86_64__) && defined(__GLIBC__) && defined(__GNUC__) \
		&& !defined(__clang__)
#define AYLP_KERNEL __attribute__((target_clones("arch=x86-64-v3","default")))
#else
#define AYLP_KERNEL
#endif


static inline double mono_us(void)
{
	struct timespec t;
	clock_gettime(CLOCK_MONOTONIC, &t);
	return (double)t.tv_sec*1e6 + (double)t.tv_nsec*1e-3;
}


// Count, sum and sum-of-squares of the window minus the core, from a sample of
// its rows rather than all of them, scaled back up.
//
// The tail enters the fit only through bg. A background taken from a few
// thousand pixels is good to a twentieth of a count -- finer than the count
// quantisation it is estimated from, and far finer than anything downstream can
// use. Sweeping all 61504 px of a 248x248 frame instead costs several
// microseconds of memory traffic on a frame that arrives cold from the camera:
// most of the loop's latency budget, spent refining a number that was already
// exact enough. Sampling makes that cost independent of sensor size and, what a
// fixed-delay controller cares about far more than the last decimal of bg,
// identical every frame.
//
// The sampled rows are taken in a few CONTIGUOUS bands rather than spread every
// d-th row. Both read the same number of bytes, but a band is a sequential run
// the hardware prefetcher will follow, where one row in eight is a stride of
// two kilobytes that it will not -- and on a cold frame that is the whole cost.
// Bands are spread across the window so a background gradient still shows up.
//
// Row accumulators are 32 bits, which is exact for any row up to 66043 px.
#define TAIL_BANDS 4

static inline size_t tail_nband(size_t win_h, size_t max_rows)
{
	return win_h <= max_rows ? 1 : TAIL_BANDS;
}

// Rows [i0,i1) of band k, as offsets into the window.
static inline void tail_band(size_t win_h, size_t max_rows, size_t k,
	size_t *i0, size_t *i1
) {
	size_t nband = tail_nband(win_h, max_rows);
	size_t band_h = win_h <= max_rows ? win_h : max_rows/nband;
	if (!band_h) band_h = 1;
	*i0 = k*win_h/nband;
	*i1 = *i0 + band_h < win_h ? *i0 + band_h : win_h;
}

AYLP_KERNEL
static void tail_sums(const unsigned char *base, size_t tda,
	size_t org_y, size_t org_x, size_t win_h, size_t win_w,
	size_t core_y, size_t core_x, size_t core_h, size_t core_w,
	size_t max_rows, double *n, double *s1, double *s2
) {
	uint64_t a = 0, b = 0, cnt = 0;
	size_t nband = tail_nband(win_h, max_rows);
	for (size_t k = 0; k < nband; k++) {
		size_t i0, i1;
		tail_band(win_h, max_rows, k, &i0, &i1);
		for (size_t i = i0; i < i1; i++) {
			size_t row = org_y + i;
			const unsigned char *r = base + row*tda;
			// The columns of this row that are not in the core: the
			// whole row, or the two pieces either side of it.
			size_t seg[2][2], nseg = 0;
			if (row >= core_y && row < core_y + core_h) {
				if (core_x > org_x) {
					seg[nseg][0] = org_x;
					seg[nseg][1] = core_x;
					nseg++;
				}
				if (core_x + core_w < org_x + win_w) {
					seg[nseg][0] = core_x + core_w;
					seg[nseg][1] = org_x + win_w;
					nseg++;
				}
			} else {
				seg[nseg][0] = org_x;
				seg[nseg][1] = org_x + win_w;
				nseg++;
			}
			for (size_t m = 0; m < nseg; m++) {
				uint32_t ra = 0, rb = 0;
				for (size_t j = seg[m][0]; j < seg[m][1]; j++) {
					uint32_t v = r[j];
					ra += v;
					rb += v*v;
				}
				a += ra; b += rb;
				cnt += seg[m][1] - seg[m][0];
			}
		}
	}
	// Scale the sample up to the population it stands for, so the cost is on
	// the same footing as when every window pixel was iterated on.
	double want = (double)(win_h*win_w - core_h*core_w);
	double scale = cnt ? want/(double)cnt : 0.0;
	*n = want;
	*s1 = scale*(double)a;
	*s2 = scale*(double)b;
}


// Brightest pixel of a rectangle. Two vectorizable passes -- max, then locate
// the row that holds it -- rather than one scalar pass carrying an index,
// because the index dependency is what stops a scan like this vectorizing at
// all. Only ever called on acquisition.
AYLP_KERNEL
static unsigned char rect_argmax(const unsigned char *base, size_t tda,
	size_t y0, size_t x0, size_t h, size_t w, size_t *py, size_t *px
) {
	unsigned char best = 0;
	for (size_t i = 0; i < h; i++) {
		const unsigned char *row = base + (y0+i)*tda + x0;
		for (size_t j = 0; j < w; j++)
			if (row[j] > best) best = row[j];
	}
	*py = y0; *px = x0;
	for (size_t i = 0; i < h; i++) {
		const unsigned char *row = base + (y0+i)*tda + x0;
		unsigned hit = 0;
		for (size_t j = 0; j < w; j++) hit |= (row[j] == best);
		if (!hit) continue;
		for (size_t j = 0; j < w; j++) {
			if (row[j] != best) continue;
			*py = y0+i; *px = x0+j;
			return best;
		}
	}
	return best;
}


// Everything a kernel needs to know about this frame, so it takes one pointer
// instead of reaching back through the device struct on every row.
struct fit_ctx {
	size_t h, w, stride;
	size_t y0, x0;			// core origin, image coords
	double ref_row;
	const double *pix;		// core pixels, pad columns 0
	const double *wt;		// robust weights, pad columns 0
	double tn, ts, ts2;		// window minus core: n, sum(I), sum(I^2)
};


// Gaussian sampled on the pixel grid, by recurrence rather than by
// transcendental.
//
// exp(-(d+s)^2/2r^2) = exp(-d^2/2r^2) * exp(-(2ds+s^2)/2r^2), and that second
// factor itself advances by a constant ratio, so a gaussian sampled at a
// constant step costs two multiplies per sample. The step is constant along a
// row (one column), and ALSO down the core -- because the centre moves
// linearly with row index, dy and the row's first dx both advance by a
// constant. So the entire 2D core costs nine exp() calls per pass instead of
// one per pixel, and a cost-only trial pass is nearly free.
//
// The original code walked outward from the sample nearest the centre so no
// ratio could exceed 1, which a narrow sigma across a 248-wide window would
// otherwise overflow. The core box is bounded at a few sigma, so the largest
// ratio here is bounded too, and a straight left-to-right sweep -- which is
// what lets a row be contiguous vector loads -- is safe.
struct rowwalk {
	double ey, rey, cey;	// exp(-dy^2/2r^2), down the rows
	double gx, rgx, cgx;	// exp(-dx^2/2r^2) at each row's first column
	double rx, crx;		// that row's first column-to-column ratio
	double cstep, cstep6, cstep16;
	double dy, dy_step;
	double dx, dx_step;
	double t;
};

static inline void rowwalk_init(struct rowwalk *rw, const struct fit_ctx *c,
	const double *p
) {
	double sig = p[AYLP_FIT_P_SIGMA];
	double i2 = 0.5/(sig*sig);
	double a = 1.0 - p[AYLP_FIT_P_SY];	// dy gains this per row
	double b = -p[AYLP_FIT_P_SX];		// dx gains this per row
	double row0 = (double)c->y0;
	double t0 = row0 - c->ref_row;
	double dy0 = row0 - (p[AYLP_FIT_P_Y0] + p[AYLP_FIT_P_SY]*t0);
	double dx0 = (double)c->x0 - (p[AYLP_FIT_P_X0] + p[AYLP_FIT_P_SX]*t0);

	rw->ey = exp(-dy0*dy0*i2);
	rw->rey = exp(-(2.0*a*dy0 + a*a)*i2);
	rw->cey = exp(-2.0*a*a*i2);
	rw->gx = exp(-dx0*dx0*i2);
	rw->rgx = exp(-(2.0*b*dx0 + b*b)*i2);
	rw->cgx = exp(-2.0*b*b*i2);
	rw->rx = exp(-(2.0*dx0 + 1.0)*i2);
	rw->crx = exp(-2.0*b*i2);

	double cs = exp(-2.0*i2);
	double cs2 = cs*cs, cs4 = cs2*cs2, cs8 = cs4*cs4;
	rw->cstep = cs;
	rw->cstep6 = cs4*cs2;
	rw->cstep16 = cs8*cs8;

	rw->dy = dy0; rw->dy_step = a;
	rw->dx = dx0; rw->dx_step = b;
	rw->t = t0;
}

static inline void rowwalk_next(struct rowwalk *rw)
{
	rw->ey *= rw->rey; rw->rey *= rw->cey;
	rw->gx *= rw->rgx; rw->rgx *= rw->cgx;
	rw->rx *= rw->crx;
	rw->dy += rw->dy_step;
	rw->dx += rw->dx_step;
	rw->t += 1.0;
}

// Seed the four lanes of a row, and the factor that advances all four by one
// vector: g(j+VW)/g(j) = r(j)^VW * cstep^(1+2+...+(VW-1)), and that factor is
// itself scaled by cstep^(VW*VW) each step.
static inline void rowwalk_vec(const struct rowwalk *rw, v4d *g, v4d *ratio)
{
	double g0 = rw->gx, r0 = rw->rx, cs = rw->cstep;
	double g1 = g0*r0, r1 = r0*cs;
	double g2 = g1*r1, r2 = r1*cs;
	double g3 = g2*r2, r3 = r2*cs;
	v4d r = {r0, r1, r2, r3};
	v4d r2v = r*r;
	*g = (v4d){g0, g1, g2, g3};
	*ratio = r2v*r2v*VBC(rw->cstep6);
}


// The tail -- window minus core -- is flat at bg, so its whole contribution to
// the cost is a quadratic in bg over three numbers gathered once per frame.
static inline double tail_ss(const struct fit_ctx *c, double bg)
{
	return c->tn*bg*bg - 2.0*bg*c->ts + c->ts2;
}


// Sum of squares over a whole core buffer, pad columns included -- they hold a
// zero residual precisely so this can run whole vectors without masking.
// Spelled out rather than left to the compiler because a floating-point
// reduction cannot be reassociated without -ffast-math, so the obvious loop
// compiles to one serial add per element at four cycles of latency each.
AYLP_KERNEL
static double sum_sq(const double *restrict x, size_t n)
{
	v4d acc = VBC(0.0);
	for (size_t k = 0; k < n; k += VW) {
		v4d v = VLD(x+k);
		acc += v*v;
	}
	return VSUM(acc);
}


// Tukey biweight, as its root (see the note on aylp_fit_com_data::weights).
// Returns how many pixels kept any weight at all.
AYLP_KERNEL
static size_t tukey_root(const double *restrict resid, double *restrict wroot,
	size_t n, double cut
) {
	double inv_cut = 1.0/cut;
	size_t kept = 0;
	for (size_t k = 0; k < n; k++) {
		double u = resid[k]*inv_cut;
		double m = 1.0 - u*u;
		if (m < 0.0) m = 0.0;
		wroot[k] = m;
		kept += m > 0.0;
	}
	return kept;
}


// Weighted cost of parameter vector p, with the per-pixel residuals. This is
// the pass a rejected trial step pays for, so it does nothing it does not have
// to: no jacobian, no transcendental, two multiplies per pixel of gaussian.
AYLP_KERNEL
static double eval_cost(const struct fit_ctx *c, const double *p,
	double *restrict resid
) {
	struct rowwalk rw;
	rowwalk_init(&rw, c, p);
	v4d vamp = VBC(p[AYLP_FIT_P_AMP]);
	v4d vbg = VBC(p[AYLP_FIT_P_BG]);
	v4d vc16 = VBC(rw.cstep16);
	v4d acost = VBC(0.0);
	for (size_t i = 0; i < c->h; i++, rowwalk_next(&rw)) {
		v4d vey = VBC(rw.ey), g, ratio;
		rowwalk_vec(&rw, &g, &ratio);
		const double *pixr = c->pix + i*c->stride;
		const double *wr = c->wt + i*c->stride;
		double *rr = resid + i*c->stride;
		for (size_t j = 0; j < c->stride; j += VW) {
			v4d q = vey*g;
			v4d res = vamp*q + (vbg - VLD(pixr+j));
			v4d sr = VLD(wr+j)*res;
			VST(rr+j, res);
			acost += sr*sr;
			g *= ratio;
			ratio *= vc16;
		}
		// Pad columns carry weight 0 so they cannot reach the cost, but
		// their residual is meaningless; zero it so every later reduction
		// over the buffer can run whole vectors without masking.
		for (size_t j = c->w; j < c->stride; j++) rr[j] = 0.0;
	}
	return VSUM(acost) + tail_ss(c, p[AYLP_FIT_P_BG]);
}


// Cost, residuals, and the Gauss-Newton normal equations, all from one walk of
// the core. The residual is model-minus-observation, so the step solves
// (JtJ + lambda diag) d = -Jtr.
//
// Written as the Gram matrix of the whitened design matrix [J | r], because
// the obvious form -- accumulate JtJ[a][b] and Jtr[a] in the pixel loop -- asks
// for 36 vector accumulators live at once against sixteen registers, and
// spills every one of them to the stack. That turns each of the 36 multiply-adds
// into a load, a multiply, an add and a store, and it measured 5.9 ns/px.
//
// So: weighting a least-squares problem by w is the same as scaling each row of
// the design matrix and the residual by sqrt(w) and dropping the weight, and
// the Tukey weight is a square by construction (w = m^2, m >= 0) so its root is
// already to hand and none is taken here. Pass one writes those eight scaled
// columns to a scratch block small enough to stay in L1; pass two walks it once
// per column, holding that column and the <= 8 accumulators it feeds in
// registers. The residual is carried as the ninth column, so cost and Jtr fall
// out of the same triangle as JtJ instead of needing their own accumulators.
//
// All seven jacobian columns are accumulated unconditionally even when the
// slopes are frozen: solve_step overwrites the frozen rows and columns anyway,
// and a branch here would cost more than the arithmetic it skips.
// One column of pass 2: accumulate column A against every column at or after
// it. A must be a literal, so the accumulator array has a constant extent and
// lives in registers.
#define GRAM_COL(A) do { \
	v4d gacc[NCOL]; \
	for (size_t b = (A); b < NCOL; b++) gacc[b] = VBC(0.0); \
	const double *ca = scr + (size_t)(A)*blk; \
	for (size_t k = 0; k < blk; k += VW) { \
		v4d x = VLD(ca+k); \
		for (size_t b = (A); b < NCOL; b++) \
			gacc[b] += x*VLD(scr + b*blk + k); \
	} \
	for (size_t b = (A); b < NCOL; b++) \
		gram[GTRI((A),b)] += VSUM(gacc[b]); \
} while (0)

AYLP_KERNEL
static double eval_jac(const struct fit_ctx *c, const double *p,
	double *restrict JtJ, double *restrict Jtr, double *restrict resid,
	double *restrict scr
) {
	struct rowwalk rw;
	rowwalk_init(&rw, c, p);
	double sig = p[AYLP_FIT_P_SIGMA];
	double inv_s2 = 1.0/(sig*sig);
	v4d vamp = VBC(p[AYLP_FIT_P_AMP]);
	v4d vbg = VBC(p[AYLP_FIT_P_BG]);
	v4d vc16 = VBC(rw.cstep16);
	v4d vis3 = VBC(inv_s2/sig);
	v4d vstep = VBC((double)VW);
	v4d vstep_s2 = VBC((double)VW*inv_s2);
	double gram[NGRAM];
	for (size_t k = 0; k < NGRAM; k++) gram[k] = 0.0;

	size_t stride = c->stride;
	for (size_t i0 = 0; i0 < c->h; i0 += JAC_BLOCK) {
		size_t nb = c->h - i0 < JAC_BLOCK ? c->h - i0 : JAC_BLOCK;
		size_t blk = nb*stride;

		// ---- pass 1: the whitened design matrix [J | r] ----
		for (size_t i = 0; i < nb; i++, rowwalk_next(&rw)) {
			v4d vey = VBC(rw.ey), g, ratio;
			rowwalk_vec(&rw, &g, &ratio);
			v4d vt = VBC(rw.t);
			v4d vdyk = VBC(rw.dy*inv_s2);
			v4d vdy2 = VBC(rw.dy*rw.dy);
			v4d vdx = {rw.dx, rw.dx+1.0, rw.dx+2.0, rw.dx+3.0};
			v4d vdxs = vdx*VBC(inv_s2);
			const double *pixr = c->pix + (i0+i)*stride;
			const double *wr = c->wt + (i0+i)*stride;
			double *rr = resid + (i0+i)*stride;
			for (size_t j = 0; j < stride; j += VW) {
				v4d q = vey*g;
				v4d res = vamp*q + (vbg - VLD(pixr+j));
				v4d sw = VLD(wr+j);
				VST(rr+j, res);
				v4d aq = vamp*q;
				v4d su = sw*(aq*vdyk);
				v4d sv = sw*(aq*vdxs);
				double *o = scr + i*stride + j;
				VST(o + AYLP_FIT_P_Y0*blk, su);
				VST(o + AYLP_FIT_P_X0*blk, sv);
				VST(o + AYLP_FIT_P_SY*blk, su*vt);
				VST(o + AYLP_FIT_P_SX*blk, sv*vt);
				VST(o + AYLP_FIT_P_SIGMA*blk,
					sw*(aq*(vdx*vdx + vdy2)*vis3));
				VST(o + AYLP_FIT_P_AMP*blk, sw*q);
				VST(o + AYLP_FIT_P_BG*blk, sw);
				VST(o + NP*blk, sw*res);
				g *= ratio;
				ratio *= vc16;
				vdx += vstep;
				vdxs += vstep_s2;
			}
			for (size_t j = c->w; j < stride; j++) rr[j] = 0.0;
		}

		// ---- pass 2: upper triangle of the Gram matrix ----
		// Unrolled over the column index by hand. With it as a loop
		// variable the accumulator count is a runtime value, so GCC
		// keeps the accumulators on the stack and the whole point of
		// this restructuring -- holding them in registers -- is lost.
		GRAM_COL(0); GRAM_COL(1); GRAM_COL(2); GRAM_COL(3);
		GRAM_COL(4); GRAM_COL(5); GRAM_COL(6); GRAM_COL(7);
	}

	for (size_t a = 0; a < NP; a++) {
		Jtr[a] = gram[GTRI(a,NP)];
		for (size_t b = a; b < NP; b++)
			JtJ[a*NP+b] = JtJ[b*NP+a] = gram[GTRI(a,b)];
	}
	// The tail's jacobian is zero in every column but bg, where it is 1.
	double bg = p[AYLP_FIT_P_BG];
	Jtr[AYLP_FIT_P_BG] += c->tn*bg - c->ts;
	JtJ[AYLP_FIT_P_BG*NP + AYLP_FIT_P_BG] += c->tn;
	return gram[GTRI(NP,NP)] + tail_ss(c, bg);
}


// Cholesky solve of the (symmetric positive definite) damped normal
// equations, in place on a copy. Returns false if not positive definite,
// which the caller answers by raising the damping.
static bool solve_step(const double *JtJ, const double *Jtr,
	const bool *active, double lambda, double *step
) {
	double A[NP*NP];
	memcpy(A, JtJ, sizeof(A));
	for (size_t a = 0; a < NP; a++) {
		if (active[a]) {
			// Marquardt scaling: damp by the diagonal, not by the
			// identity, so parameters with wildly different natural
			// scales (a background in counts against a slope in
			// px/row) are damped comparably.
			double d = A[a*NP+a];
			A[a*NP+a] = d + lambda*(d > 0.0 ? d : 1.0);
		} else {
			// Frozen parameter: unit row/column, zero rhs => zero step.
			for (size_t b = 0; b < NP; b++) A[a*NP+b] = A[b*NP+a] = 0.0;
			A[a*NP+a] = 1.0;
		}
	}
	double L[NP*NP];
	memset(L, 0, sizeof(L));
	for (size_t a = 0; a < NP; a++) {
		for (size_t b = 0; b <= a; b++) {
			double s = A[a*NP+b];
			for (size_t k = 0; k < b; k++) s -= L[a*NP+k]*L[b*NP+k];
			if (a == b) {
				if (s <= 1e-300) return false;
				L[a*NP+a] = sqrt(s);
			} else {
				L[a*NP+b] = s/L[b*NP+b];
			}
		}
	}
	double z[NP];
	for (size_t a = 0; a < NP; a++) {
		double s = active[a] ? -Jtr[a] : 0.0;
		for (size_t k = 0; k < a; k++) s -= L[a*NP+k]*z[k];
		z[a] = s/L[a*NP+a];
	}
	for (size_t a = NP; a-- > 0; ) {
		double s = z[a];
		for (size_t k = a+1; k < NP; k++) s -= L[k*NP+a]*step[k];
		step[a] = s/L[a*NP+a];
	}
	return true;
}


// Slide the window so it lies fully inside the image.
static void clamp_window(struct aylp_fit_com_data *data,
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


static void clamp_params(struct aylp_fit_com_data *data, double *p,
	size_t max_y, size_t max_x
) {
	if (!(p[AYLP_FIT_P_SIGMA] >= data->sigma_min))
		p[AYLP_FIT_P_SIGMA] = data->sigma_min;
	if (p[AYLP_FIT_P_SIGMA] > data->sigma_max)
		p[AYLP_FIT_P_SIGMA] = data->sigma_max;
	if (!(p[AYLP_FIT_P_AMP] >= 0.0)) p[AYLP_FIT_P_AMP] = 0.0;
	if (!(p[AYLP_FIT_P_BG] >= 0.0)) p[AYLP_FIT_P_BG] = 0.0;
	if (p[AYLP_FIT_P_BG] > 255.0) p[AYLP_FIT_P_BG] = 255.0;
	if (!(p[AYLP_FIT_P_Y0] >= 0.0)) p[AYLP_FIT_P_Y0] = 0.0;
	if (p[AYLP_FIT_P_Y0] > (double)(max_y-1)) p[AYLP_FIT_P_Y0] = (double)(max_y-1);
	if (!(p[AYLP_FIT_P_X0] >= 0.0)) p[AYLP_FIT_P_X0] = 0.0;
	if (p[AYLP_FIT_P_X0] > (double)(max_x-1)) p[AYLP_FIT_P_X0] = (double)(max_x-1);
	// A slope steep enough to smear the beam across the whole window in one
	// readout is not a measurement, it is a diverging fit.
	double smax = (double)max_y;
	for (int k = AYLP_FIT_P_SY; k <= AYLP_FIT_P_SX; k++) {
		if (!(p[k] > -smax)) p[k] = -smax;
		if (p[k] > smax) p[k] = smax;
	}
}


// Clamp [centre-half, centre+half] into [lo,hi] and cap its length, shrinking
// around the centre so the beam keeps the middle of whatever survives.
static void plan_span(double centre, double half, long lo, long hi, size_t cap,
	size_t *org, size_t *len
) {
	long a = (long)floor(centre - half);
	long b = (long)ceil(centre + half);
	if (a < lo) a = lo;
	if (a > hi) a = hi;
	if (b > hi) b = hi;
	if (b < lo) b = lo;
	if (b < a) b = a;
	if ((size_t)(b - a + 1) > cap) {
		long c = (long)llround(centre);
		if (c < a) c = a;
		if (c > b) c = b;
		a = c - (long)(cap/2);
		if (a < lo) a = lo;
		b = a + (long)cap - 1;
		if (b > hi) {
			b = hi;
			a = b - (long)cap + 1;
			if (a < lo) a = lo;
		}
	}
	*org = (size_t)a;
	*len = (size_t)(b - a + 1);
}


// Choose the box the solver iterates on, from the warm-started parameters.
static void plan_core(struct aylp_fit_com_data *data,
	size_t org_y, size_t org_x
) {
	const double *p = data->p;
	double sig = p[AYLP_FIT_P_SIGMA];
	if (!(sig >= data->sigma_min) || sig > data->sigma_max)
		sig = data->sigma_init;
	// One extra pixel of slack, so a centre that has moved a fraction of a
	// pixel since the previous frame cannot shave the edge off the box.
	double rad = data->fit_radius*sig + 1.0;
	if (rad < 3.0) rad = 3.0;

	// Rows carrying signal are those with |dy| <= rad. dy is linear in the
	// row index with slope (1 - slope_y), so as the beam tracks the readout
	// the band widens; if it tracks it exactly the band is the whole frame,
	// hence the cap. A band that wide is a diverging slope, not a
	// measurement, and the gate downstream will say so.
	double a = 1.0 - p[AYLP_FIT_P_SY];
	if (fabs(a) < 0.25) a = a < 0.0 ? -0.25 : 0.25;
	double row_c = (p[AYLP_FIT_P_Y0] - p[AYLP_FIT_P_SY]*data->ref_row)/a;
	double half_y = rad/fabs(a);
	if (half_y > 4.0*rad) half_y = 4.0*rad;

	// The x centre sweeps by slope_x across that band of rows.
	double col_c = p[AYLP_FIT_P_X0]
		+ p[AYLP_FIT_P_SX]*(row_c - data->ref_row);
	double half_x = rad + fabs(p[AYLP_FIT_P_SX])*half_y;
	if (half_x > 4.0*rad) half_x = 4.0*rad;

	plan_span(row_c, half_y, (long)org_y, (long)(org_y + data->win_h - 1),
		data->max_core, &data->core_y, &data->core_h);
	plan_span(col_c, half_x, (long)org_x, (long)(org_x + data->win_w - 1),
		data->max_core, &data->core_x, &data->core_w);
	// Widen to a whole vector where the window and the cap allow. The
	// kernels run whole vectors either way and pad the remainder with zero
	// weights, so those lanes are paid for regardless -- better that they
	// carry real pixels than nothing.
	size_t want = (data->core_w + VW - 1) / VW * VW;
	if (want <= data->win_w && want <= data->max_core
			&& want > data->core_w) {
		size_t room_l = data->core_x - org_x;
		size_t room_r = (org_x + data->win_w)
			- (data->core_x + data->core_w);
		size_t grow = want - data->core_w;
		size_t left = grow/2;
		if (left > room_l) left = room_l;
		if (grow - left > room_r) left = grow - room_r;
		data->core_x -= left;
		data->core_w = want;
	}
	data->core_stride = (data->core_w + VW - 1) / VW * VW;

	size_t need = data->core_h * data->core_stride;
	if (UNLIKELY(need > data->core_cap
			|| data->core_stride > data->core_cap_stride)) {
		data->core_cap_stride = data->core_stride;
		xfree(data->pix);
		xfree(data->weights);
		xfree(data->resid);
		xfree(data->resid_alt);
		xfree(data->scratch);
		data->core_cap = need;
		size_t bytes = (need*sizeof(double) + 31) / 32 * 32;
		data->pix = alloc_check(aligned_alloc(32, bytes));
		data->weights = alloc_check(aligned_alloc(32, bytes));
		data->resid = alloc_check(aligned_alloc(32, bytes));
		data->resid_alt = alloc_check(aligned_alloc(32, bytes));
		size_t sbytes = (NCOL*JAC_BLOCK*data->core_stride*sizeof(double)
			+ 31) / 32 * 32;
		data->scratch = alloc_check(aligned_alloc(32, sbytes));
	}
}


// Ask for every byte of the frame this device will touch, all at once.
//
// The frame has just been DMAed in by the camera, so it is cold in every cache,
// and at 248x248 the misses -- not the arithmetic -- were the largest single
// item in the frame's latency. But which bytes are wanted is known before any
// of them is needed: the core box is already planned, and the tail's bands are
// fixed by the window. Issuing the loads up front lets the misses overlap each
// other and the work that follows, instead of each one stalling the loop that
// asks for it.
static void prefetch_frame(const struct aylp_fit_com_data *data,
	const unsigned char *base, size_t tda, size_t org_y, size_t org_x
) {
	for (size_t k = 0; k < tail_nband(data->win_h, data->tail_rows); k++) {
		size_t i0, i1;
		tail_band(data->win_h, data->tail_rows, k, &i0, &i1);
		for (size_t i = i0; i < i1; i++) {
			const unsigned char *r = base + (org_y+i)*tda + org_x;
			for (size_t j = 0; j < data->win_w; j += 64)
				__builtin_prefetch(r + j);
		}
	}
	for (size_t i = 0; i < data->core_h; i++)
		for (size_t j = 0; j < data->core_w; j += 64)
			__builtin_prefetch(base + (data->core_y+i)*tda
				+ data->core_x + j);
}


// Copy the core out of the frame as doubles, once, so the solver's passes are
// contiguous aligned loads instead of a widening conversion per pixel per
// pass.
AYLP_KERNEL
static void gather_core(struct aylp_fit_com_data *data,
	const unsigned char *base, size_t tda
) {
	for (size_t i = 0; i < data->core_h; i++) {
		const unsigned char *row = base + (data->core_y+i)*tda
			+ data->core_x;
		double *dst = data->pix + i*data->core_stride;
		double *w = data->weights + i*data->core_stride;
		for (size_t j = 0; j < data->core_w; j++) {
			dst[j] = (double)row[j];
			w[j] = 1.0;
		}
		for (size_t j = data->core_w; j < data->core_stride; j++) {
			dst[j] = 0.0;
			w[j] = 0.0;
		}
	}
}


// Seed the parameter vector from the image itself: brightest pixel (or the
// configured init point) for position, window border for background, and the
// flux-weighted second moment for sigma. Slopes start at zero.
static void acquire(struct aylp_fit_com_data *data, gsl_matrix_uchar *img)
{
	size_t max_y = img->size1, max_x = img->size2;
	if (data->init_y >= 0) {
		data->win_y = (size_t)data->init_y;
		data->win_x = (size_t)data->init_x;
	} else {
		size_t by, bx;
		rect_argmax(img->data, img->tda, 0, 0, max_y, max_x, &by, &bx);
		data->win_y = by; data->win_x = bx;
	}
	clamp_window(data, max_y, max_x);
	size_t org_y = data->win_y - data->win_h/2;
	size_t org_x = data->win_x - data->win_w/2;

	// background from the window border, peak from the interior
	double bsum = 0.0; size_t bn = 0;
	for (size_t j = 0; j < data->win_w; j++) {
		bsum += img->data[org_y*img->tda + org_x+j];
		bsum += img->data[(org_y+data->win_h-1)*img->tda + org_x+j];
		bn += 2;
	}
	for (size_t i = 1; i + 1 < data->win_h; i++) {
		bsum += img->data[(org_y+i)*img->tda + org_x];
		bsum += img->data[(org_y+i)*img->tda + org_x+data->win_w-1];
		bn += 2;
	}
	double bg = bn ? bsum/(double)bn : 0.0;
	size_t py, px;
	unsigned char peak = rect_argmax(img->data, img->tda, org_y, org_x,
		data->win_h, data->win_w, &py, &px);
	if (data->init_y >= 0) { py = (size_t)data->init_y; px = (size_t)data->init_x; }

	// sigma from the radial second moment: for a 2D gaussian,
	// <r^2> = 2 sigma^2 over the flux above background.
	// Taken LOCALLY, over a box a few sigma_init wide around the peak. Over
	// the whole window a stray reflection inflates it badly (its distance
	// enters as r^2), the initial model is then broad enough to span both
	// blobs, and the fit settles between them -- at which point residuals
	// are large at BOTH and no amount of robust reweighting can tell which
	// one is the beam. Seeding narrow and on the peak keeps the stray
	// firmly in the tail, where the reweighting can discard it.
	long rad = (long)ceil(3.0*data->sigma_init);
	long lo_y = (long)py - rad, hi_y = (long)py + rad;
	long lo_x = (long)px - rad, hi_x = (long)px + rad;
	if (lo_y < (long)org_y) lo_y = (long)org_y;
	if (lo_x < (long)org_x) lo_x = (long)org_x;
	if (hi_y > (long)(org_y+data->win_h-1)) hi_y = (long)(org_y+data->win_h-1);
	if (hi_x > (long)(org_x+data->win_w-1)) hi_x = (long)(org_x+data->win_w-1);
	double sw = 0.0, sr2 = 0.0;
	for (long i = lo_y; i <= hi_y; i++) {
		double dy = (double)i - (double)py;
		for (long j = lo_x; j <= hi_x; j++) {
			double dx = (double)j - (double)px;
			double v = (double)img->data[i*(long)img->tda + j] - bg;
			if (v <= 0.0) continue;
			sw += v; sr2 += v*(dy*dy + dx*dx);
		}
	}
	double sigma = sw > 0.0 ? sqrt(sr2/sw/2.0) : data->sigma_init;
	if (!(sigma >= data->sigma_min)) sigma = data->sigma_init;
	if (sigma > data->sigma_max) sigma = data->sigma_max;

	data->p[AYLP_FIT_P_Y0] = (double)py;
	data->p[AYLP_FIT_P_X0] = (double)px;
	data->p[AYLP_FIT_P_SY] = 0.0;
	data->p[AYLP_FIT_P_SX] = 0.0;
	data->p[AYLP_FIT_P_SIGMA] = sigma;
	data->p[AYLP_FIT_P_AMP] = (double)peak - bg > 0.0 ? (double)peak - bg : 1.0;
	data->p[AYLP_FIT_P_BG] = bg;
	data->lambda = 1e-3;
	data->acquired = true;
	data->lost = 0;
	log_info("fit_com: acquired at (%.2f,%.2f), sigma %.2f, amp %.1f, bg %.1f",
		data->p[AYLP_FIT_P_Y0], data->p[AYLP_FIT_P_X0], sigma,
		data->p[AYLP_FIT_P_AMP], bg);
}


int fit_com_init(struct aylp_device *self)
{
	self->device_data = xcalloc(1, sizeof(struct aylp_fit_com_data));
	struct aylp_fit_com_data *data = self->device_data;
	self->fini = &fit_com_fini;

	// defaults
	data->window_h = 0;		// 0 => whole image
	data->window_w = 0;
	data->init_y = -1;
	data->init_x = -1;
	data->max_iter = 10;
	data->tol = 1e-4;
	data->max_us = 10.0;
	data->robust_k = 2.5;
	data->robust_iter = 1;
	data->sigma_init = 2.0;
	data->sigma_min = 0.5;
	data->sigma_max = 20.0;
	data->fit_radius = 3.5;
	data->max_core = 64;
	data->tail_rows = 32;
	data->min_amplitude = 5.0;
	data->max_residual = 0.0;	// 0 => disabled
	data->max_step = 0.0;		// 0 => disabled
	data->moment_output = false;
	data->fit_gaussian = true;
	data->moment_cut = 5.0;
	data->moment_col_stride = 1;
	data->moment_max_y_skew = 0.0;
	data->moment_min_y_width = 0.0;
	data->pwm_period_frames = 0.0;
	data->pwm_dark_flux = 100.0;
	data->pwm_bright_flux = 1000.0;
	data->pwm_full_start = 4;
	data->pwm_full_end_margin = 2;
	data->pwm_qy = 0.03;
	data->pwm_qx = 0.03;
	data->pwm_ry_full = 1.0;
	data->pwm_rx_full = 1.0;
	data->reacquire_after = 10;
	data->row_time = 0.0;
	bool fit_slope = true;

	if (!self->params) {
		log_error("No params object found.");
		return -1;
	}
	json_object_object_foreach(self->params, key, val) {
		if (key[0] == '_') {
		} else if (!strcmp(key, "window_height")) {
			data->window_h = json_object_get_uint64(val);
		} else if (!strcmp(key, "window_width")) {
			data->window_w = json_object_get_uint64(val);
		} else if (!strcmp(key, "init_y")) {
			data->init_y = json_object_get_int64(val);
		} else if (!strcmp(key, "init_x")) {
			data->init_x = json_object_get_int64(val);
		} else if (!strcmp(key, "max_iter")) {
			data->max_iter = json_object_get_uint64(val);
		} else if (!strcmp(key, "max_us")) {
			data->max_us = json_object_get_double(val);
		} else if (!strcmp(key, "tol")) {
			data->tol = json_object_get_double(val);
		} else if (!strcmp(key, "robust_k")) {
			data->robust_k = json_object_get_double(val);
		} else if (!strcmp(key, "robust_iter")) {
			data->robust_iter = json_object_get_uint64(val);
		} else if (!strcmp(key, "sigma_init")) {
			data->sigma_init = json_object_get_double(val);
		} else if (!strcmp(key, "sigma_min")) {
			data->sigma_min = json_object_get_double(val);
		} else if (!strcmp(key, "sigma_max")) {
			data->sigma_max = json_object_get_double(val);
		} else if (!strcmp(key, "fit_radius")) {
			data->fit_radius = json_object_get_double(val);
		} else if (!strcmp(key, "max_core")) {
			data->max_core = json_object_get_uint64(val);
		} else if (!strcmp(key, "tail_rows")) {
			data->tail_rows = json_object_get_uint64(val);
		} else if (!strcmp(key, "min_amplitude")) {
			data->min_amplitude = json_object_get_double(val);
		} else if (!strcmp(key, "max_residual")) {
			data->max_residual = json_object_get_double(val);
		} else if (!strcmp(key, "max_step")) {
			data->max_step = json_object_get_double(val);
		} else if (!strcmp(key, "moment_output")) {
			data->moment_output = json_object_get_boolean(val);
		} else if (!strcmp(key, "fit_gaussian")) {
			data->fit_gaussian = json_object_get_boolean(val);
		} else if (!strcmp(key, "moment_cut")) {
			data->moment_cut = json_object_get_double(val);
		} else if (!strcmp(key, "moment_col_stride")) {
			data->moment_col_stride = json_object_get_uint64(val);
		} else if (!strcmp(key, "moment_max_y_skew")) {
			data->moment_max_y_skew = json_object_get_double(val);
		} else if (!strcmp(key, "moment_min_y_width")) {
			data->moment_min_y_width = json_object_get_double(val);
		} else if (!strcmp(key, "moment_reject_inferred")) {
			data->moment_reject_inferred = json_object_get_boolean(val);
		} else if (!strcmp(key, "pwm_period_frames")) {
			data->pwm_period_frames = json_object_get_double(val);
		} else if (!strcmp(key, "pwm_dark_flux")) {
			data->pwm_dark_flux = json_object_get_double(val);
		} else if (!strcmp(key, "pwm_bright_flux")) {
			data->pwm_bright_flux = json_object_get_double(val);
		} else if (!strcmp(key, "pwm_full_start")) {
			data->pwm_full_start = json_object_get_uint64(val);
		} else if (!strcmp(key, "pwm_full_end_margin")) {
			data->pwm_full_end_margin = json_object_get_uint64(val);
		} else if (!strcmp(key, "reacquire_after")) {
			data->reacquire_after = json_object_get_uint64(val);
		} else if (!strcmp(key, "row_time")) {
			data->row_time = json_object_get_double(val);
		} else if (!strcmp(key, "fit_slope")) {
			fit_slope = json_object_get_boolean(val);
		} else {
			log_warn("Unknown parameter \"%s\"", key);
		}
	}

	if ((data->init_y < 0) != (data->init_x < 0)) {
		log_error("Provide both init_y and init_x, or neither");
		return -1;
	}
	if (!data->fit_gaussian && !data->moment_output) {
		log_error("fit_com: fit_gaussian=false requires moment_output=true");
		return -1;
	}
	if (data->fit_gaussian && !data->max_iter) {
		log_error("fit_com: max_iter must be nonzero");
		return -1;
	}
	if (!isfinite(data->tol) || data->tol < 0.0) {
		log_error("fit_com: tol must be finite and non-negative");
		return -1;
	}
	if (!isfinite(data->max_us) || data->max_us < 0.0) {
		log_error("fit_com: max_us must be finite and non-negative");
		return -1;
	}
	if (!isfinite(data->robust_k) || data->robust_k < 0.0) {
		log_error("fit_com: robust_k must be finite and non-negative");
		return -1;
	}
	if (!isfinite(data->sigma_min) || !isfinite(data->sigma_max)
			|| data->sigma_min < 0.05
			|| data->sigma_max <= data->sigma_min) {
		// The gaussian recurrence spans the core box in steps of one
		// pixel; below about a twentieth of a pixel of width that span
		// is enough decades to overflow, and a sub-pixel beam is not
		// something a centroid on this sensor can resolve anyway.
		log_error("fit_com: need 0.05 <= sigma_min < sigma_max");
		return -1;
	}
	if (!isfinite(data->sigma_init) || data->sigma_init <= 0.0) {
		log_error("fit_com: sigma_init must be positive");
		return -1;
	}
	if (!isfinite(data->fit_radius) || data->fit_radius < 2.0) {
		log_error("fit_com: fit_radius must be at least 2 sigmas");
		return -1;
	}
	if (data->max_core < 8) {
		log_error("fit_com: max_core must be at least 8");
		return -1;
	}
	if (!data->tail_rows) {
		log_error("fit_com: tail_rows must be nonzero");
		return -1;
	}
	if (!isfinite(data->min_amplitude) || data->min_amplitude < 0.0) {
		log_error("fit_com: min_amplitude must be finite, non-negative");
		return -1;
	}
	if (!isfinite(data->max_residual) || data->max_residual < 0.0) {
		log_error("fit_com: max_residual must be finite, non-negative");
		return -1;
	}
	if (!isfinite(data->max_step) || data->max_step < 0.0) {
		log_error("fit_com: max_step must be finite, non-negative");
		return -1;
	}
	if (!isfinite(data->moment_cut) || data->moment_cut < 0.0) {
		log_error("fit_com: moment_cut must be finite and non-negative");
		return -1;
	}
	if (!data->moment_col_stride) {
		log_error("fit_com: moment_col_stride must be nonzero");
		return -1;
	}
	if (!isfinite(data->moment_max_y_skew)
			|| data->moment_max_y_skew < 0.0) {
		log_error("fit_com: moment_max_y_skew must be finite and non-negative");
		return -1;
	}
	if (!isfinite(data->moment_min_y_width)
			|| data->moment_min_y_width < 0.0
			|| data->moment_min_y_width > 1.0) {
		log_error("fit_com: moment_min_y_width must be in [0, 1]");
		return -1;
	}
	if (!isfinite(data->pwm_period_frames) || data->pwm_period_frames < 0.0
			|| !isfinite(data->pwm_dark_flux) || data->pwm_dark_flux < 0.0
			|| !isfinite(data->pwm_bright_flux)
			|| data->pwm_bright_flux <= data->pwm_dark_flux) {
		log_error("fit_com: invalid PWM period/flux thresholds");
		return -1;
	}
	if (!isfinite(data->row_time) || data->row_time < 0.0) {
		log_error("fit_com: row_time must be finite and non-negative");
		return -1;
	}
	if (!data->reacquire_after) {
		log_error("fit_com: reacquire_after must be nonzero");
		return -1;
	}
	if (data->window_h && data->window_h < 3) {
		log_error("fit_com: window_height must be at least 3");
		return -1;
	}
	if (data->window_w && data->window_w < 3) {
		log_error("fit_com: window_width must be at least 3");
		return -1;
	}

	for (size_t a = 0; a < NP; a++) data->active[a] = true;
	data->active[AYLP_FIT_P_SY] = fit_slope;
	data->active[AYLP_FIT_P_SX] = fit_slope;
	data->lambda = 1e-3;

	self->type_in = AYLP_T_MATRIX_UCHAR;
	self->units_in = AYLP_U_ANY;
	self->type_out = AYLP_T_VECTOR;
	self->units_out = AYLP_U_MINMAX;
	self->proc = &fit_com_proc;
	data->com = xmalloc_type(gsl_vector, 2);

	log_info("fit_com: %s model, window %zux%zu (0 = whole frame), "
		"core <= %zu px at %.1f sigma, max_iter %zu, max_us %.1f",
		fit_slope ? "sheared (row-dependent centre)" : "static",
		data->window_h, data->window_w, data->max_core,
		data->fit_radius, data->max_iter, data->max_us);
	return 0;
}


int fit_com_proc(struct aylp_device *self, struct aylp_state *state)
{
	struct aylp_fit_com_data *data = self->device_data;
	data->inferred_last = false;
	gsl_matrix_uchar *img = state->matrix_uchar;
	if (UNLIKELY(!img || !img->data || img->size1 < 3 || img->size2 < 3
			|| img->tda < img->size2)) {
		log_error("fit_com: invalid image matrix or dimensions");
		return -1;
	}
	double t_start = data->max_us > 0.0 ? mono_us() : 0.0;
	size_t max_y = img->size1, max_x = img->size2;

	// (re)size the window to this image
	size_t want_h = data->window_h ? data->window_h : max_y;
	size_t want_w = data->window_w ? data->window_w : max_x;
	if (want_h > max_y) want_h = max_y;
	if (want_w > max_x) want_w = max_x;
	if (UNLIKELY(want_h != data->win_h || want_w != data->win_w)) {
		data->win_h = want_h;
		data->win_w = want_w;
		data->acquired = false;
	}
	// A fixed epoch in IMAGE coordinates, so the instant the reported
	// position belongs to does not wander as the window tracks the beam --
	// a moving epoch reads downstream as a wobbling loop delay.
	data->ref_row = 0.5*((double)max_y - 1.0);

	bool fresh = !data->acquired;
	if (UNLIKELY(fresh)) acquire(data, img);
	clamp_window(data, max_y, max_x);
	size_t org_y = data->win_y - data->win_h/2;
	size_t org_x = data->win_x - data->win_w/2;
	size_t npix = data->win_h * data->win_w;

	// Fix the evaluated box for the whole frame, from the warm-started
	// parameters. Fixed, because a box that moved between a trial step and
	// the point it is compared against would be comparing two different
	// cost functions.
	plan_core(data, org_y, org_x);
	prefetch_frame(data, img->data, img->tda, org_y, org_x);
	gather_core(data, img->data, img->tda);
	// Detect a rolling-shutter illumination edge against the previous complete
	// Gaussian. Uniform brightness changes scale every row equally; a shutter
	// edge makes the observed/expected ratio collapse on only part of the spot.
	// Mask those rows before fitting so their zeros cannot pull the inferred
	// full-beam centre toward the visible fragment.
	unsigned char shutter_row[data->core_h];
	memset(shutter_row, 0, sizeof shutter_row);
	double row_ratio[data->core_h], ratio_max = 0.0;
	for (size_t i = 0; i < data->core_h; i++) {
		double obs = 0.0;
		double y = data->core_y+i;
		for (size_t j = 0; j < data->core_w; j++) {
			double v = data->pix[i*data->core_stride+j]
				- data->p[AYLP_FIT_P_BG];
			if (v > 0.0) obs += v;
		}
		double dy = y-data->p[AYLP_FIT_P_Y0];
		double row_shape = exp(-dy*dy/(2.0*data->p[AYLP_FIT_P_SIGMA]
			*data->p[AYLP_FIT_P_SIGMA]));
		row_ratio[i] = row_shape > 0.01 ? obs/row_shape : -1.0;
		if (row_ratio[i] > ratio_max) ratio_max = row_ratio[i];
	}
	size_t shutter_low = 0, shutter_high = 0;
	if (!fresh && ratio_max > 0.0)
		for (size_t i = 0; i < data->core_h; i++) {
			if (row_ratio[i] < 0.0) continue;
			if (row_ratio[i] < 0.10*ratio_max) shutter_low++;
			else shutter_high++;
		}
	bool shutter_masked = false;
	if (shutter_masked)
		for (size_t i = 0; i < data->core_h; i++)
			if (row_ratio[i] >= 0.0 && row_ratio[i] < 0.10*ratio_max) {
				shutter_row[i] = 1;
				memset(data->weights+i*data->core_stride, 0,
					data->core_w*sizeof *data->weights);
			}
	tail_sums(img->data, img->tda, org_y, org_x, data->win_h, data->win_w,
		data->core_y, data->core_x, data->core_h, data->core_w,
		data->tail_rows, &data->tail_n, &data->tail_s, &data->tail_s2);

	struct fit_ctx ctx = {
		.h = data->core_h, .w = data->core_w,
		.stride = data->core_stride,
		.y0 = data->core_y, .x0 = data->core_x,
		.ref_row = data->ref_row,
		.pix = data->pix, .wt = data->weights,
		.tn = shutter_masked ? 0 : data->tail_n,
		.ts = data->tail_s,
		.ts2 = data->tail_s2,
	};
	size_t nbuf = data->core_h * data->core_stride;
	size_t npad = data->core_h * (data->core_stride - data->core_w);

	// Warm start from the previous frame: between frames the beam moves far
	// less than its own width, so the fit begins inside the basin and
	// converges in a couple of iterations.
	double p[NP];
	memcpy(p, data->p, sizeof(p));
	bool frame_active[NP];
	memcpy(frame_active, data->active, sizeof frame_active);
	if (shutter_masked) {
		frame_active[AYLP_FIT_P_SIGMA] = false;
		frame_active[AYLP_FIT_P_BG] = false;
		frame_active[AYLP_FIT_P_SY] = false;
		frame_active[AYLP_FIT_P_SX] = false;
	}

	double JtJ[NP*NP], Jtr[NP], step[NP], trial[NP];
	double *resid = data->resid, *resid_alt = data->resid_alt;
	double cost = eval_cost(&ctx, p, resid);

	// The robust weighting needs residuals from a model that is already
	// close, or a redescending weight rejects the beam instead of the
	// artifact. robust_iter buys that with a plain first iteration.
	//
	// Starting it at iteration 0 instead is tempting -- the initial pass
	// above already evaluated the WARM-STARTED parameters, and on a clean
	// synthetic stray that rejects it in one iteration rather than four.
	// But those warm-start residuals also carry the frame's noise, and
	// cutting on them before the fit has seated down-weights one flank of
	// the beam more than the other: measured over a moving, noisy scene it
	// costs a factor of five in position rms. So the plain iteration stays.
	//
	// What the opening residual IS good for is deciding whether the warm
	// start can be trusted at all. On the frame that just acquired there is
	// none, and if the beam has jumped the previous solution no longer
	// describes this frame; either way the plain iteration is not optional,
	// whatever robust_iter was configured to.
	double rms0 = sqrt((cost > 0.0 ? cost : 0.0)/(double)npix);
	size_t rob_start = data->robust_iter;
	if ((fresh || !(rms0 <= 2.0*data->last_rms)) && rob_start < 1)
		rob_start = 1;

	double lambda = data->lambda;
	size_t it = 0;
	bool converged = false, budget_hit = false;
	double iter_us = 0.0;
	for (; data->fit_gaussian && it < data->max_iter; it++) {
		// A latency guard, not an accuracy one: the loop is warm-started
		// and an ordinary frame converges in two or three iterations
		// well inside the budget, so what this bounds is the frame that
		// does not -- see doc/devices/fit_com.md for what that costs on
		// a scene that genuinely wants more. Jitter in measurement latency
		// feeds straight into fsp's fixed-delay model, and a late
		// measurement is worth less to the controller than a slightly
		// less converged one delivered on time.
		if (data->max_us > 0.0 && it) {
			// Predict the next iteration from the last, with margin:
			// an iteration that needs extra damping attempts costs
			// more than the one before it, and the epilogue -- the
			// residual reduction and the gates -- still has to run.
			double used = mono_us() - t_start;
			if (used + 1.5*iter_us > data->max_us) {
				budget_hit = true;
				break;
			}
		}
		double t_iter = data->max_us > 0.0 ? mono_us() : 0.0;
		// IRLS: once the fit has seated, downweight pixels the model
		// cannot explain -- a stray reflection, a hot pixel, clipping.
		// Crucially the model already contains the motion, so a sheared
		// beam fits and is NOT treated as an outlier.
		if (data->robust_k > 0.0 && it >= rob_start) {
			// The scale is the rms over the whole WINDOW, as it was
			// when every window pixel was iterated on: the tail's
			// contribution is the same closed form in bg. Keeping
			// that unchanged matters, because the scale is what sets
			// how aggressive the cutoff is.
			double ss = sum_sq(resid, nbuf)
				+ tail_ss(&ctx, p[AYLP_FIT_P_BG]);
			double scale = sqrt(ss/(double)npix);
			if (scale > 1e-9) {
				// Tukey biweight, not Huber. Huber only tapers as
				// cut/|r| and never reaches zero, so a bright
				// stray reflection -- which is not a sparse
				// outlier but hundreds of pixels of coherent
				// support -- keeps enough leverage to drag the
				// centre several px. A redescending weight
				// excludes it outright once the fit is seated on
				// the real beam, which acquisition guarantees by
				// starting from the brightest pixel.
				double cut = data->robust_k*scale;
				size_t kept = tukey_root(resid, data->weights,
					nbuf, cut);
				// Pad columns have a zero residual, so the loop
				// above just gave them weight 1. Never let
				// reweighting starve the fit of data either.
				bool starved = kept - npad < 4*NP;
				for (size_t i = 0; i < data->core_h; i++) {
					double *w = data->weights
						+ i*data->core_stride;
					if (shutter_row[i])
						for (size_t j = 0; j < data->core_w; j++)
							w[j] = 0.0;
					else if (starved)
						for (size_t j = 0;
								j < data->core_w; j++)
							w[j] = 1.0;
					for (size_t j = data->core_w;
							j < data->core_stride; j++)
						w[j] = 0.0;
				}
			}
		}
		// One pass: cost at the current point under the current weights,
		// the normal equations, and the residuals -- all from the same
		// walk over the core rather than three separate ones.
		cost = eval_jac(&ctx, p, JtJ, Jtr, resid, data->scratch);
		bool stepped = false;
		for (int attempt = 0; attempt < 8; attempt++) {
			if (!solve_step(JtJ, Jtr, frame_active, lambda, step)) {
				lambda *= 10.0;
				continue;
			}
			for (size_t a = 0; a < NP; a++)
				trial[a] = p[a] + (frame_active[a] ? step[a] : 0.0);
			clamp_params(data, trial, max_y, max_x);
			// Trial residuals go to the spare buffer, so an accepted
			// step just swaps pointers instead of costing another
			// full pass to re-derive what it already computed.
			double tc = eval_cost(&ctx, trial, resid_alt);
			if (tc < cost) {
				double rel = (cost - tc)/(cost > 0.0 ? cost : 1.0);
				memcpy(p, trial, sizeof(p));
				double *sw = resid; resid = resid_alt; resid_alt = sw;
				cost = tc;
				lambda = lambda > 1e-8 ? lambda*0.3 : lambda;
				stepped = true;
				// A fit that has never been reweighted has not
				// converged -- it has converged to the UNWEIGHTED
				// problem, which is a different one, and on a
				// scene with an artifact a materially wrong one.
				// Without this the loop can stop at the first
				// iteration and the robust weighting never runs
				// at all: the cost is dominated by background
				// pixels the fit cannot improve, so the relative
				// improvement is tiny however far the centre
				// still has to move.
				if (rel < data->tol && it >= rob_start)
					converged = true;
				break;
			}
			lambda *= 10.0;
			if (lambda > 1e12) break;
		}
		if (data->max_us > 0.0) iter_us = mono_us() - t_iter;
		if (!stepped || converged) break;
	}
	data->resid = resid;
	data->resid_alt = resid_alt;
	data->lambda = lambda < 1e-8 ? 1e-8 : (lambda > 1.0 ? 1.0 : lambda);
	// iterations actually run: the loop leaves `it` at the one it broke out
	// of, or at max_iter if it ran them all
	data->n_iter_last = it < data->max_iter ? it + 1 : it;
	data->budget_hit = budget_hit;

	// unweighted rms residual over the window, the honest goodness-of-fit
	double ss = sum_sq(resid, nbuf) + tail_ss(&ctx, p[AYLP_FIT_P_BG]);
	double rms = sqrt((ss > 0.0 ? ss : 0.0)/(double)npix);

	// A centre pinned to an image boundary is a constrained optimizer
	// result, not a measured centroid.  Require one fitted sigma of support
	// on every side so a clipped spot cannot be reported as valid motion.
	double centre_margin = p[AYLP_FIT_P_SIGMA];
	double step_y = p[AYLP_FIT_P_Y0] - data->p[AYLP_FIT_P_Y0];
	double step_x = p[AYLP_FIT_P_X0] - data->p[AYLP_FIT_P_X0];
	bool fit_ok = data->fit_gaussian
		&& isfinite(p[AYLP_FIT_P_Y0]) && isfinite(p[AYLP_FIT_P_X0])
		&& p[AYLP_FIT_P_Y0] >= centre_margin
		&& p[AYLP_FIT_P_Y0] <= (double)(max_y-1) - centre_margin
		&& p[AYLP_FIT_P_X0] >= centre_margin
		&& p[AYLP_FIT_P_X0] <= (double)(max_x-1) - centre_margin
		&& (fresh || data->max_step <= 0.0
			|| hypot(step_y, step_x) <= data->max_step)
		&& p[AYLP_FIT_P_AMP] >= data->min_amplitude
		&& p[AYLP_FIT_P_SIGMA] > data->sigma_min
		&& p[AYLP_FIT_P_SIGMA] < data->sigma_max
		&& (data->max_residual <= 0.0 || rms <= data->max_residual);

	// PWM illumination plus a rolling shutter can cut a bright spot at a row
	// boundary. The raw intensity centroid of that fragment is not the beam
	// centre. Use moments on complete frames, but when row skew/width identifies
	// a shutter fragment publish the robust Gaussian's inferred full-spot centre.
	double out_y = p[AYLP_FIT_P_Y0], out_x = p[AYLP_FIT_P_X0];
	double moment_sum = 0.0, moment_peak = 0.0;
	bool moment_ok = false, fragment = false;
	double prev_y = (data->last_y+1.0)*0.5*(max_y-1.0);
	double prev_x = (data->last_x+1.0)*0.5*(max_x-1.0);
	if (data->moment_have_full) data->moment_full_age++;
	if (data->moment_output) {
		double row_flux[data->win_h];
		memset(row_flux, 0, sizeof row_flux);
		double bs = 0.0; size_t bn = 0;
		for (size_t i = 0; i < data->win_h; i++)
			for (size_t j = 0; j < data->win_w; j++)
				if (i == 0 || j == 0 || i == data->win_h-1
						|| j == data->win_w-1) {
					bs += img->data[(org_y+i)*img->tda + org_x+j]; bn++;
				}
		double bg = bn ? bs/bn : 0.0, sy = 0.0, sx = 0.0;
		double syy = 0.0, syyy = 0.0;
		double sw2 = 0.0, sy2 = 0.0, sx2 = 0.0;
		for (size_t i = 0; i < data->win_h; i++) {
			double row_sum = 0.0;
			for (size_t j = 0; j < data->win_w;
					j += data->moment_col_stride) {
				double v = (double)img->data[(org_y+i)*img->tda + org_x+j] - bg;
				if (v > moment_peak) moment_peak = v;
				if (v > 2.0) {
					sw2 += v; sy2 += v*(org_y+i); sx2 += v*(org_x+j);
					row_sum += v;
				}
				if (v <= data->moment_cut) continue;
				moment_sum += v;
				double ay = (double)(org_y+i);
				sy += v*ay; sx += v*(org_x+j);
				syy += v*ay*ay; syyy += v*ay*ay*ay;
			}
			row_flux[i] = row_sum;
		}
		if (moment_sum > 100.0) { out_y = sy/moment_sum; out_x = sx/moment_sum; }
		double y2 = 0.0, y3 = 0.0;
		if (moment_sum > 100.0) {
			// Derive central moments from the raw moments accumulated above.
			// This removes a full-window pass from the sensor-to-command path.
			y2 = syy - 2.0*out_y*sy + out_y*out_y*moment_sum;
			y3 = syyy - 3.0*out_y*syy + 3.0*out_y*out_y*sy
				- out_y*out_y*out_y*moment_sum;
			if (y2 < 0.0 && y2 > -1e-9*syy) y2 = 0.0;
		}
		double moment_y_width = moment_sum > 0.0 ? sqrt(y2/moment_sum) : 0.0;
		double moment_y_skew = y2 > 0.0
			? fabs((y3/moment_sum) / pow(y2/moment_sum, 1.5)) : INFINITY;
		double out_y2 = sw2 > 100.0 ? sy2/sw2 : out_y;
		double out_x2 = sw2 > 100.0 ? sx2/sw2 : out_x;
		double margin = data->sigma_init;
		double threshold_shift = hypot(out_y-out_y2, out_x-out_x2);
		if (data->pwm_period_frames > 0.0) {
			/* A dark onset is a much more reliable shutter-phase marker than
			 * spot skew: this beam is intrinsically skewed even when fully on.
			 * Frames safely inside the on plateau are calibration measurements;
			 * transition frames get flux-dependent measurement noise and dark
			 * frames are pure causal predictions. */
			bool dark = sw2 < data->pwm_dark_flux;
			if (dark && !data->pwm_prev_dark) {
				data->pwm_phase = 0;
				data->pwm_phase_locked = true;
			} else if (data->pwm_phase_locked) {
				data->pwm_phase++;
			}
			data->pwm_prev_dark = dark;
			size_t last_full = (size_t)floor(data->pwm_period_frames)
				- data->pwm_full_end_margin;
			bool full = data->pwm_phase_locked
				&& data->pwm_phase >= data->pwm_full_start
				&& data->pwm_phase <= last_full;
			bool bright = sw2 >= data->pwm_bright_flux;
			bool measured = full && bright && sw2 > 100.0;
			if (measured) {
				data->pwm_y=out_y2; data->pwm_x=out_x2;
				data->pwm_filter_have=true;
			}
			fragment = !measured;
			moment_ok = data->pwm_filter_have;
			if (moment_ok) { out_y=data->pwm_y; out_x=data->pwm_x; }
		} else {
		fragment = shutter_masked || (moment_sum > 100.0
			&& moment_peak >= data->min_amplitude
			&& (threshold_shift > 1.5
				|| (data->moment_max_y_skew > 0.0
					&& moment_y_skew > data->moment_max_y_skew)
				|| (data->moment_min_y_width > 0.0
					&& moment_y_width < data->moment_min_y_width
						* data->sigma_init)));
		moment_ok = moment_sum > 100.0 && moment_peak >= data->min_amplitude
			&& sw2 > 100.0 && threshold_shift <= 1.5
			&& out_y >= margin && out_y <= max_y-1.0-margin
			&& out_x >= margin && out_x <= max_x-1.0-margin
			&& (fresh || data->max_step <= 0.0
				|| hypot(out_y-prev_y, out_x-prev_x) <= data->max_step);
		if (fragment) {
			double rfmax = 0.0;
			for (size_t i = 0; i < data->win_h; i++)
				if (row_flux[i] > rfmax) rfmax = row_flux[i];
			double loy = fmax(margin, prev_y-4.0);
			double hiy = fmin(max_y-1.0-margin, prev_y+4.0);
			double measured_x = sx/moment_sum;
			double best_y = prev_y, best_x = measured_x, best_ss = INFINITY,
				best_ff = 0.0;
			size_t nr = 0;
			for (size_t i = 0; i < data->win_h; i++)
				if (row_flux[i] >= 0.10*rfmax) nr++;
			// Register the observed row profile against the latest complete
			// frame. Brightness is a fitted scale and therefore cannot move it.
			for (double yc=loy; yc<=hiy+0.25 && nr>=3; yc+=0.5) {
				double tf=0.0, tt=0.0, ff=0.0;
				for (size_t i=0; i<data->win_h; i++) {
					if (row_flux[i] < 0.10*rfmax) continue;
					double q=(org_y+i)-(yc-data->moment_template_y)
						-data->moment_template_org_y;
					long q0=(long)floor(q); double t=0.0, f=q-q0;
					if (q0>=0 && (size_t)(q0+1)<data->moment_template_n)
						t=(1-f)*data->moment_row_template[q0]
							+f*data->moment_row_template[q0+1];
					tf+=t*row_flux[i]; tt+=t*t; ff+=row_flux[i]*row_flux[i];
				}
				double ssr=tt>0.0 ? ff-tf*tf/tt : INFINITY;
				if (ssr<best_ss) { best_ss=ssr; best_ff=ff; best_y=yc; }
			}
			double profile_rel = best_ff > 0.0 ? best_ss/best_ff : INFINITY;
			if (data->moment_have_full) {
				out_y = data->moment_full_y + 0.25*(best_y-data->moment_full_y);
				out_x = data->moment_full_x + 0.25*(best_x-data->moment_full_x);
			} else { out_y = best_y; out_x = best_x; }
			moment_ok = nr >= 3 && isfinite(best_ss)
				&& profile_rel < 0.50
				&& out_y >= margin && out_y <= max_y-1.0-margin
				&& out_x >= margin && out_x <= max_x-1.0-margin
				&& (fresh || data->max_step <= 0.0
					|| hypot(out_y-prev_y,out_x-prev_x)<=data->max_step);
		}
		if (!fragment && moment_ok) {
			if (data->moment_template_cap < data->win_h) {
				data->moment_row_template = xrealloc(data->moment_row_template,
					data->win_h*sizeof *data->moment_row_template);
				data->moment_template_cap = data->win_h;
			}
			memcpy(data->moment_row_template, row_flux,
				data->win_h*sizeof *row_flux);
			data->moment_template_n = data->win_h;
			data->moment_template_org_y = org_y;
			data->moment_template_y = out_y;
		}
		}
	}
	bool ok = data->moment_output ? moment_ok : fit_ok;

	if (LIKELY(ok)) {
		data->inferred_last = fragment;
		if (data->moment_output && !fragment) {
			// Full frames alone update motion; a fragment prediction must never
			// train on itself. Divide displacement by elapsed frames: complete
			// frames can be separated by several PWM-sliced readouts.
			if (data->moment_have_full && data->moment_full_age) {
				double dy = (out_y-data->moment_full_y)/data->moment_full_age;
				double dx = (out_x-data->moment_full_x)/data->moment_full_age;
				data->moment_dy = 0.8*data->moment_dy + 0.2*dy;
				data->moment_dx = 0.8*data->moment_dx + 0.2*dx;
			}
			data->moment_full_y = out_y; data->moment_full_x = out_x;
			data->moment_full_age = 0; data->moment_have_full = true;
		}
		if (fit_ok && !fragment) memcpy(data->p, p, sizeof(p));
		data->last_rms = rms;
		data->vel_y = data->row_time > 0.0
			? p[AYLP_FIT_P_SY]/data->row_time : 0.0;
		data->vel_x = data->row_time > 0.0
			? p[AYLP_FIT_P_SX]/data->row_time : 0.0;
		data->last_y = -1.0 + 2.0*out_y/((double)max_y - 1.0);
		data->last_x = -1.0 + 2.0*out_x/((double)max_x - 1.0);
		data->win_y = (size_t)llround(out_y);
		data->win_x = (size_t)llround(out_x);
		clamp_window(data, max_y, max_x);
		data->lost = 0;
		state->header.status &= (aylp_status)~AYLP_FRAME_REJECTED;
		state->header.status &= (aylp_status)~AYLP_BEAM_LOST;
		if (fragment && data->moment_reject_inferred)
			state->header.status |= AYLP_FRAME_REJECTED;
	} else {
		// The model could not describe this frame. Hold the last valid
		// output rather than publish a diverged fit, and say so.
		state->header.status |= AYLP_FRAME_REJECTED;
		if (++data->lost >= data->reacquire_after) {
			state->header.status |= AYLP_BEAM_LOST;
			if (data->lost == data->reacquire_after) {
				log_warn("fit_com: fit failed for %zu frames "
					"(amp %.1f < %.1f? sigma %.2f, rms %.1f); "
					"re-acquiring", data->lost,
					p[AYLP_FIT_P_AMP], data->min_amplitude,
					p[AYLP_FIT_P_SIGMA], rms);
			}
			data->acquired = false;
		}
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


int fit_com_fini(struct aylp_device *self)
{
	struct aylp_fit_com_data *data = self->device_data;
	if (data) {
		xfree(data->pix);
		xfree(data->weights);
		xfree(data->resid);
		xfree(data->resid_alt);
		xfree(data->scratch);
		xfree(data->moment_row_template);
		if (data->com) xfree_type(gsl_vector, data->com);
	}
	xfree(self->device_data);
	return 0;
}
