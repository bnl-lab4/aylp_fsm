anyloop:fit_com
================

Types and units: `[T_MATRIX_UCHAR, U_ANY] -> [T_VECTOR, U_MINMAX]`.

A **joint spatio-temporal beam fit**. Output convention (normalized to the
whole image, `[-1, 1]`, 0 = image centre) matches `anyloop:center_of_mass` and
`anyloop:wfs_com` exactly, so it is a drop-in replacement upstream of
`anyloop:fsp`/`anyloop:pid`.

Why
---

A rolling-shutter frame does not have *a* beam position. Each row is exposed at
a different instant, so if the beam moves during readout the frame is sheared,
and any single position extracted from it is ill-defined — there are as many
answers as there are rows.

`center_of_mass` ignores this and reports the first moment of the smear.
`wfs_com` measures a position and then tries to undo the shear afterwards, by
fitting a line through 2–3 row-group means — two parameters from three points,
with no redundancy, so the fitted slope is mostly noise (see
`doc/devices/wfs_com.md`, which documents that correction measuring net-harmful
and being disabled by default).

This device puts the motion **inside the model** and fits it directly:

```
I(row,col) = bg + amp * exp(-[(row - y(row))^2 + (col - x(row))^2] / 2 sigma^2)

  y(row) = y0 + slope_y * (row - ref_row)
  x(row) = x0 + slope_x * (row - ref_row)
```

Seven parameters — `y0, x0, slope_y, slope_x, sigma, amp, bg` — against every
pixel in the window. On a 20x20 window that is 400 observations for 7 unknowns,
against the row regression's 3 for 2. The same physical quantity, estimated
with two orders of magnitude more data.

The consequence that matters: **shear stops being a defect to detect and reject,
and becomes the signal that identifies the slope.** Nothing is ever discarded
for moving. That removes an entire failure mode, because a shear threshold is
necessarily a *velocity* threshold — shear is proportional to speed — so it
preferentially drops frames at the fast part of every oscillation and keeps the
turning points. Those dropouts are correlated with the disturbance phase, which
means the loop runs open (holding its last command) exactly when the beam is
moving fastest, and any adaptive predictor downstream identifies its model from
a phase-biased sample.

Measured against `wfs_com` on the same synthetic sheared scene (3 px @ 60 Hz,
32x32 ROI, fs 3788):

| | raw RMS | debiased RMS | bias | us/frame |
|---|---|---|---|---|
| `fit_com` | 0.0043 px | 0.0043 px | (-0.000, -0.000) | 7.7 |
| `wfs_com` | 0.6602 px | 0.1099 px | (+0.456, +0.465) | 13 |

The accuracy columns are from the original survey. `fit_com`'s timing was 62
us/frame when that survey was run and is re-stated here for the current
implementation, measured as described under *Performance*; `wfs_com`'s is
unchanged code and unchanged number.

The bias column is the second structural difference. `wfs_com` measures
displacement relative to a reference captured at acquisition, so whatever
offset the beam had at that instant is absorbed and reported as zero. This
device solves for absolute position every frame, so there is no anchor to
absorb anything and no acquisition-time offset to carry.

Algorithm
---------

Per frame:

1. Place the window (whole image by default) around the previous position.
2. Warm-start the parameter vector from the previous frame. Between frames the
   beam moves far less than its own width, so the solver begins inside the
   basin and typically converges in one or two iterations.
3. Levenberg–Marquardt with Marquardt (diagonal) scaling, so parameters with
   very different natural scales — a background in counts against a slope in
   px/row — are damped comparably. The damping is carried between frames.
4. From the second iteration, reweight pixels by a **Tukey biweight** on the
   residuals: a stray reflection or hot pixel is excluded outright. Because the
   model already contains the motion, a *sheared* beam fits and is never
   mistaken for an outlier.
5. Gate the result. A fit whose amplitude has collapsed, whose sigma has hit a
   bound, or whose residual exceeds `max_residual` is rejected: the last valid
   output is held and `AYLP_FRAME_REJECTED` asserted. After `reacquire_after`
   consecutive failures, `AYLP_BEAM_LOST` is asserted and the parameters are
   re-seeded from the image.

Acquisition seeds position from the brightest pixel (or `init_y`/`init_x`),
background from the window border, and sigma from a **local** radial second
moment taken over a box a few `sigma_init` wide. That locality matters: over
the whole window a stray reflection inflates the sigma seed (its distance
enters as r²), the initial model is then broad enough to span both blobs, and
the fit settles between them — at which point residuals are large at *both* and
no amount of robust reweighting can tell which one is the beam. Seeding narrow
and on the peak keeps the stray in the tail where the reweighting can discard
it. Measured: a stray 7 px away pulls the centre 2.41 px without this, 0.27 px
with it.

Performance
-----------

At 248x248, whole frame, on a Coffee Lake i9 (measured over a ring of
pre-rendered frames large enough to fall out of last-level cache, which is what
the device actually sees -- a frame just DMAed in by the camera, cold, while its
own scratch stays warm):

| | p50 | p90 | position rms |
|---|---|---|---|
| before | 7401 us | 7819 us | 0.0216 px |
| now | 6.3 us | 7.2 us | 0.0212 px |

Same answer, about 1200x less time. Four things get it there.

**The solver iterates on a core box, not on the window.** Past a few sigma the
model is flat to well below one count, so those pixels say nothing about
position, width or amplitude -- only about the background. And there the model
*is* a constant, so `sum((bg - I)^2)` over them is `n*bg^2 - 2*bg*sum(I) +
sum(I^2)`: three numbers, gathered once per frame, that stand in for every one
of those pixels at every iteration, exactly. This is what decouples cost from
sensor size -- a 248x248 frame and a 32x32 one now do the same solver work --
and it is the single biggest factor. The box is re-planned every frame from the
warm-started sigma and slopes, including the widening a sheared beam needs, so
it always contains the beam.

**Those three background numbers are sampled, not swept.** A background taken
from a few thousand pixels is good to a twentieth of a count, finer than the
quantisation it is estimated from. Sweeping all 61504 instead cost several
microseconds of memory traffic on a cold frame -- most of the budget, spent
refining a number that was already exact enough. The sampled rows are taken in
a few contiguous bands rather than every d-th row: same bytes, but a band is a
sequential run the prefetcher will follow where a 2 kB stride is not.

**The frame is prefetched.** Everything the device will read is known before any
of it is needed, so it is all asked for at once and the misses overlap each
other instead of each stalling the loop that wants it.

**The kernels are vectorised, and the normal equations are built as a Gram
matrix.** The gaussian is sampled by recurrence in both axes -- the centre moves
linearly with row index, so `dy` and each row's first `dx` advance by a
constant, and the whole 2D core costs nine `exp()` calls per pass instead of one
per pixel. The jacobian pass is the expensive one, and the obvious way to write
it wants 36 vector accumulators live against sixteen registers; every one spills,
and each multiply-add becomes a load, a multiply, an add and a store. Instead,
each row of the design matrix and its residual are scaled by `sqrt(w)` -- which
a Tukey weight already has, being a square by construction -- and then `JtJ`,
`Jtr` and the cost are all just the upper triangle of the Gram matrix of
`[J | r]`, accumulated one column at a time with its accumulators in registers.
That alone was 1.8x.

The wide kernels are selected at load time by an ifunc, so the binary still runs
where AVX2 and FMA are absent -- about 2.5x slower there, which the latency
guard below then bounds instead.

Bounding the latency
--------------------

`max_us` (default 10) is a hard cap on the whole call. Before each iteration the
solver predicts the next one from the last, with margin, and stops if it would
overrun; `budget_hit` in the device data records that it did.

**This trades convergence for punctuality, and on a hard scene the trade is
real.** Measured on a stray reflection 7 px from the beam -- the case that most
wants iterations -- the robust weighting needs about four LM iterations to
redescend it and nine to settle:

| iterations | centre error |
|---|---|
| 3 | 2.34 px |
| 4 | 0.39 px |
| 6 | 0.28 px |
| 9 | 0.27 px |

At roughly 2 us an iteration those nine do not fit in 10 us, so with the guard
on that scene lands nearer 3 iterations than 9. On an ordinary scene it costs
nothing -- the warm-started fit converges in two iterations and the guard never
fires. Set `max_us` to 0 to let the solver run to `max_iter` instead, which is
the right choice if a late measurement is cheaper for you than a less converged
one. It is the wrong choice for a fixed-delay model downstream, which is why the
default is the cap.

One invariant is not negotiable to the guard: **a fit that has never been
reweighted has not converged**, it has converged to the unweighted problem. The
solver will not declare convergence before the robust weighting has run at
least once, whatever the cost improvement says. Without that it could stop at
the first iteration and never reweight at all -- the cost is dominated by
background pixels the fit cannot improve, so the relative improvement is tiny
however far the centre still has to move.

What this does **not** do
--------------------------

- **Fit a beam that isn't gaussian.** This is the real exposure, and it is
  different in kind from a correlation tracker's. `wfs_com` correlates against
  a *learned* template, so it does not care what shape the beam is, only that
  the shape is stable. This device assumes the shape. Where the beam is
  speckled, clipped, or strongly astigmatic, the model is wrong everywhere
  rather than at identifiable outlier pixels, and the result is a systematic
  bias in `y0`/`x0` that robust weighting cannot remove. Check the residual map
  on real frames before trusting it: if residuals are noise-like the model
  fits; if they show structure that moves, they do not.
- **Survive two comparable blobs in the window.** Robust weighting has a
  breakdown point. The defence is the window — keep it tight enough to exclude
  known stray reflections rather than relying on the weighting.
- **Guarantee a fixed iteration count.** LM is iterative. `max_iter` and
  `max_us` bound it, but the count genuinely varies frame to frame -- what is
  bounded is the time, not the number of steps.
- **Meet 10 us without AVX2 and FMA.** The figures above are with the wide
  kernels; the baseline path an older CPU takes is about 2.5x slower, and there
  `max_us` will be cutting the solver off at one or two iterations rather than
  bounding something that already fits.
- **Resolve shear finer than the model.** The motion model is linear in row
  index — constant velocity within one readout. Genuinely curved intra-frame
  motion is not captured, and neither is a bottom-up or column-wise rolling
  shutter without changing the axis convention.

Parameters
----------

- `window_height`, `window_width` (default 0 = whole image): fit window,
  tracked on the previous position. Smaller is faster and excludes distant
  artifacts; must comfortably contain the beam (±4 sigma is a good target).
- `init_y`, `init_x` (optional): initial position and window placement. Unlike
  `wfs_com`, these only seed the search — the reported position is absolute, so
  no acquisition-time offset is absorbed.
- `sigma_init` (default 2.0): starting beam width, and the scale of the local
  box used to seed sigma at acquisition. Set it near the true value.
- `sigma_min`, `sigma_max` (defaults 0.5, 20.0): bounds. Hitting a bound
  rejects the frame, so they double as a sanity gate.
- `min_amplitude` (default 5.0): fitted amplitude below this means no beam.
- `max_residual` (default 0, disabled): rms residual in counts above which the
  frame is rejected. Set it from the residual you actually measure, plus
  headroom, or leave it off.
- `max_iter` (default 10): LM iteration cap.
- `max_us` (default 10): latency cap for the whole call, in microseconds. 0
  disables it and leaves `max_iter` as the only bound. See *Bounding the
  latency* above for what it costs on a scene that wants more iterations.
- `fit_gaussian` (default true): run the Gaussian LM solver. Set false only
  with `moment_output: true` to use the intensity-moment/PWM tracker without
  allowing a Gaussian fit to affect validity or output. Acquisition, moment
  gates, held-frame flags and beam-loss/reacquisition remain active.
- `moment_col_stride` (default 1): sample every Nth image column in the
  intensity-moment pass. Every sensor row remains represented for shutter
  detection, while broad-beam latency falls approximately with N. Validate
  position noise on raw frames before raising it; `1` uses every pixel.
- `tol` (default 1e-4): relative cost improvement at which to stop. Note this is
  relative to the *whole window's* cost, which on a large frame is dominated by
  background pixels the fit cannot improve, so it fires earlier than its face
  value suggests.
- `fit_radius` (default 3.5): half-width, in sigmas, of the box the solver
  actually iterates on. Everything outside it is treated as flat background,
  exactly (see *Performance*). At 3.5 sigma the truncated model is 0.4 counts
  for a 200-count beam -- below the quantisation. Cost scales as its square, so
  this is the knob to turn if you need latency and can afford a coarser tail;
  measured position rms was unchanged from 2.5 to 4.0 sigma on a clean beam.
- `max_core` (default 64): hard cap on either side of that box, so a runaway
  sigma cannot buy an unbounded pass. A beam wider than about 9 sigma-pixels is
  truncated by this rather than by `fit_radius`.
- `tail_rows` (default 32): how many window rows the background sample is taken
  from. Only affects the precision of `bg` and of the reported rms, both of
  which are far finer than needed at the default; lower it if the frame is very
  large and you want the gather cheaper still.
- `robust_k` (default 2.5): Tukey cutoff in residual sigmas. 0 disables
  reweighting. Smaller rejects more aggressively.
- `robust_iter` (default 1): plain iterations before reweighting begins. The
  weighting needs residuals from a model that is already close, or a
  redescending weight rejects the beam rather than the artifact. Setting it to 0
  reweights straight from the warm-started residuals, which on a clean synthetic
  stray rejects it in one iteration rather than four -- but those residuals also
  carry the frame's noise, and cutting on them before the fit has seated
  down-weights one flank of the beam more than the other: measured over a
  moving, noisy scene it cost a factor of five in position rms. Leave it at 1.
  Regardless of the setting, the first iteration is forced to run unweighted on
  the frame that just acquired, and on any frame whose opening residual says the
  previous solution no longer describes it -- in both cases there is no warm
  start to trust.
- `fit_slope` (default true): fit the intra-frame motion. False freezes both
  slopes at zero, reducing this to a static 5-parameter beam fit.
- `row_time` (default 0): seconds per image row of readout. **The position fit
  does not need it** — the slope is fitted in px/row and the reported position
  is the model evaluated at the ROI centre row, both independent of it. It
  converts the fitted slope to a physical velocity and fixes the epoch the
  reported position belongs to, which is what a delay-modelling controller
  downstream needs. A wrong value cannot corrupt the position.
- `reacquire_after` (default 10): consecutive failed fits before
  `AYLP_BEAM_LOST` and a re-seed.

The reference epoch is the ROI centre row, fixed in **image** coordinates, so
the instant the reported position belongs to does not wander as the window
tracks the beam. A moving epoch reads downstream as a wobbling loop delay.

Diagnostics kept in the device data, useful when logging: `vel_y`/`vel_x`
(px/s, zero unless `row_time` is set), `last_rms`, and `n_iter_last`.
