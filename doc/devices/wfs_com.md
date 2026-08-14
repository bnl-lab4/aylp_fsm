anyloop:wfs_com
================

Types and units: `[T_MATRIX_UCHAR, U_ANY] -> [T_VECTOR, U_MINMAX]`.

A Shack-Hartmann-style **matched-filter** replacement for `anyloop:center_of_mass`'s
`track` mode. Where `center_of_mass` sums the whole window into a single
flux-weighted first moment, `wfs_com` splits the window into a grid of
subapertures, estimates each subaperture's own shift by normalized
cross-correlation against a self-calibrating reference template, and combines
the subapertures into one global `[y, x]` position by a confidence-weighted
average. Output convention (normalized to the whole image, `[-1, 1]`, 0 =
image centre) matches `center_of_mass` exactly, so it is a drop-in
replacement upstream of `anyloop:fsp`/`anyloop:pid`.

Why: a single first moment cannot distinguish "the beam moved" from "the beam
got asymmetrically distorted" -- both shift the flux-weighted average
identically. An array of local, independently-confidence-scored shift
estimates can: a genuine translation moves every subaperture the same way,
while distortion (partial scintillation, a developing speckle, a stray
reflection entering one corner) shows up as subapertures disagreeing, and
disagreeing/low-confidence subapertures are automatically down-weighted
rather than silently biasing the result the way they would in a plain sum.

Algorithm
---------

Per subaperture, per frame:

1. Extract an `(subap_height + 2*search_radius)` by
   `(subap_width + 2*search_radius)` block (the subaperture plus a search
   margin on every side), threshold it the same way `center_of_mass` does.
2. If this subaperture has no reference template yet (first frame after
   acquisition, or after a beam-loss re-acquire), seed the template from this
   frame's core cell and move on -- no shift estimate is possible yet.
3. Otherwise, exhaustively search every integer offset within
   `+-search_radius` pixels, scoring each by normalized cross-correlation
   (NCC) against the reference template. The offset with the highest NCC is
   the integer-pixel shift estimate; a full two-dimensional quadratic fit to
   the surrounding 3x3 NCC surface, including its cross-axis term, refines it
   to sub-pixel accuracy. NCC doubles as the confidence score (0..1): a
   well-matched, undistorted subaperture scores near 1, a poorly-matched one
   (partial occlusion, a speckle that wasn't there when the template was
   learned, mismatched shape) scores low.
4. Subapertures below `min_confidence` are excluded from this frame's
   reconstruction entirely, and their template is left untouched -- a bad
   match cannot corrupt its own reference. Subapertures at or above threshold
   contribute a weight of `confidence * matched-candidate flux` to a weighted
   average of the
   per-subaperture shift estimates, which becomes the frame's global `(dy,
   dx)`.
5. Once that consensus shift is known -- and only if the frame passes the
   validity and shear gates below -- a second pass re-registers each
   contributing subaperture's template against it: the observation is sampled
   at the consensus displacement and stored at the template's nominal index,
   which divides the motion back out, so the template stays pinned to the
   window lattice (the invariant the absolute reconstruction depends on) while
   still EWMA-adapting (`ref_beta`) to slow changes in beam *shape* such as
   focus drift.

   The registration shift is a *trimmed* consensus: subapertures whose own
   estimate sits more than half a pixel from the frame's weighted mean are
   dropped and the rest re-averaged, and those same outliers keep their
   templates frozen for that frame. This is a consistency tolerance on the
   rigid-translation model, not a motion limit -- a real translation of any
   size moves every subaperture alike and produces no disagreement. Trimming
   matters because the reported position, being a flux-weighted average over
   everything valid, is unavoidably pulled by a persistent asymmetric
   distortion; writing that pull into the references would make it outlast
   the distortion that caused it.

   The consensus, rather than each subaperture's own estimate, is what
   re-registers templates because per-subaperture registration is unstable.
   Registering a template with the same noisy estimate that template produced
   closes a private positive-feedback loop, and the correlation surface of any
   subaperture that sees only a smooth monotone gradient -- the outer cells of
   any beam that does not fill the grid -- is nearly degenerate: its NCC is
   ~1 at *every* offset, so its peak position is noise while its
   `min_confidence` score stays near-perfect. Feeding that noise back ratchets
   the template a little further out every frame until the peak pins to the
   search boundary. On a byte-identical static frame this fabricated roughly
   0.013 px/frame of drift, then failed the validity gate on every subsequent
   frame, then asserted `AYLP_BEAM_LOST` and re-acquired about every 34
   frames -- on a perfectly good, motionless beam. Since the device's whole
   model is that real motion is *common mode* across subapertures, the
   consensus is the only registration signal averaged over enough
   subapertures to be stable.
6. The global shift is a residual relative to the current integer window and
   acquisition template—not an increment from the previous output. Absolute
   position is reconstructed afresh as `window centre + acquisition anchor
   offset + residual`; accumulating the residual would repeatedly count a
   persistent sub-pixel offset and create artificial drift. The reconstructed
   position becomes the normalized output and, rounded/clamped, the next
   frame's window placement.

If fewer than `min_valid_subaps` subapertures clear `min_confidence` in a
frame (an occlusion, a very dim frame, or a chopper's "off" phase), the last
valid output is held and `AYLP_FRAME_REJECTED` is asserted rather than
reporting a spurious position, exactly as
`center_of_mass` holds through a zero-signal frame. After `reacquire_after`
such frames in a row, `AYLP_BEAM_LOST` is asserted, every subaperture's
reference template is reset, and the window re-acquires from the brightest
pixel of the whole image (or `init_y`/`init_x`, if given) -- old templates
describe pixel content at the old window location and are meaningless at a
new one.

Rolling-shutter correction (`rolling_shutter`)
----------------------------------------

With `rolling_shutter: true`, subaperture rows are no longer blended together
regardless of when they were actually captured. Instead: each row group's
confidence-weighted shift is computed on its own, a weighted linear fit
of (row-group shift) against the row group's matched-candidate mean row index
is taken independently for y and
x, and the reported position is that fit evaluated at a fixed reference time
(the window's vertical centre) rather than at the ill-defined "average"
instant a naive blend implicitly reports. This assumes the standard
top-row-first rolling shutter (row 0 captured first, later rows captured
progressively later) and constant velocity across one frame's readout, which
is the same small-motion assumption the rest of the tracker already makes.
`rolling_shutter: false` (the default) disables this and reproduces the plain
confidence-weighted average across every subaperture regardless of row.

**Measured efficacy: do not enable this without checking it helps.** Scored
against a shear-free control run of the same scene, and *debiased* (see below),
the correction is neutral at 2-10 px / 30-60 Hz and actively harmful at
10 px / 250 Hz. The cause is structural: the fit has only `subap_rows` data
points, of which a beam that does not fill the window lights 2-3, so the
intra-frame slope carries little signal and a lot of noise -- and the
extrapolation to the centre row multiplies that noise by a lever arm of several
pixels. Before the significance shrinkage described below, this more than
doubled the error (0.251 -> 0.519 px) on the fast scene. It is now clamped to
roughly break-even, but "harmless" is not "useful": enable it only with enough
subaperture rows that the beam lights at least three of them, and verify
against a control.

Two guards keep it from doing damage. The slope is only fitted when at least
three row groups have data -- a two-point line is exact by construction, so its
slope is pure noise -- and each fitted slope is then shrunk by
`b²/(b² + var(b))`, which leaves a well-determined slope essentially untouched
and collapses a noise-dominated one to zero, degrading gracefully to the plain
weighted mean.

A measurement caveat that cost real confusion: the raw RMS of this device on a
moving target is dominated by a *constant* offset whenever `init_y`/`init_x` are
set, because the beam's displacement at the acquisition instant is absorbed and
reappears as a fixed bias on every later frame. That bias has nothing to do with
rolling shutter and swamps the shear term -- measured against the raw total, the
correction appeared to deliver ~13% when the honest figure was zero. Always
debias, and always compare against a shear-free control, before concluding
anything about this correction.

The legacy `row_time` parameter's magnitude does not matter -- only whether it
is positive.
This was assumed otherwise in an earlier draft of this doc and disproven by
actually running the calibration sweep in
`contrib/diagnostics/test_wfs_com_rolling_shutter.c`: every positive value tested (0.25x
through 2x the true per-row skew) produced *identical* tracking error, and
only `row_time <= 0` (correction disabled) differed. The reason is algebraic,
not empirical noise: the regression fits `shift` against `row_time * row
index`, and evaluates the fit at `row_time * (window-centre row index)`.
Both the fitted independent variable and the evaluation point scale by the
same `row_time` factor, which cancels exactly in the extrapolation --
`row_time` currently only selects *which reconstruction* runs (naive blend
vs. per-row-group linear extrapolation to the window's vertical centre, in
row-index units), not a genuine physical time axis. Any positive value (e.g.
`1.0`) enables the correction identically to a datasheet-accurate one. This
is a real limitation to be aware of, not just a convenience: because no
actual physical time enters the calculation, this correction cannot (yet)
be combined with anything that needs a true timestamp -- e.g. reconciling
row-group measurements against the loop's own external clock. Extending it
to do that would break the current scale-invariance and reintroduce a real
calibration requirement.

Rolling-shutter rejection (`max_row_shear`)
--------------------------------------------

Independent of, and usable without, the rolling-shutter correction above:
`max_row_shear > 0` makes the device **reject** shear-distorted frames
outright rather than trying to correct them. Every frame, each row group's
confidence-weighted shift is computed (the same per-row-group numbers the
`row_time` regression uses, computed regardless of `row_time`). If the
spread between row groups — `max - min` across all row groups with data —
exceeds `max_row_shear` pixels in **either** y or x, the frame is rejected:
the last valid output is held, and `AYLP_FRAME_REJECTED` is set in the state
header for downstream logic. A running count is kept in `shear_rejected` and
each rejection is logged at debug level with the measured spread, so a
threshold can be tuned by watching the log.

Rejection is deliberately **not** counted toward `reacquire_after` and never
sets `AYLP_BEAM_LOST`: the beam is still there and the window geometry is
still valid, only this frame's readout timing is untrustworthy. Counting it
as signal loss would make a *sustained* rolling-shutter condition (continuous
vibration, exactly the case this targets) reset the window and every
reference template over and over.

More than that, a rejected frame **resets** the consecutive-loss counter,
exactly as a normally tracked frame does. A frame can only be rejected for
shear after passing the `min_valid_subaps` gate, so reaching that branch is
positive evidence the beam is present and well-matched. Merely freezing the
counter instead is not enough: genuinely low-signal frames interleaved with
shear-rejected ones would then accumulate toward `reacquire_after` across
frames where the beam was repeatedly seen just fine, and eventually fire a
spurious re-acquire. This is not hypothetical — before the counter reset was
added, the regression suite's plain rolling-shutter scenario (no occlusions,
no chopping, just shear) triggered two spurious re-acquires, and the resulting
template resets cut its shear-detection rate from 44/72 frames to 28/72.

Reference-template updates are transactional at frame scope. The EWMA
re-registration runs in a second pass that is only reached once a frame has
passed both the minimum-valid-subaperture gate and the shear gate, so a
rejected frame never reaches it at all; template-initialized flags and the
seeds written with them, which necessarily happen in the first pass, are
rolled back if either gate fails. Thus a rejected distorted frame cannot
slowly teach that distortion to the matched filters or leave a late-seeded
subaperture marked initialized with a rolled-back template. The initial
warm-up frame is the exception because its purpose is to seed the templates.

Downstream FSP recognizes `AYLP_FRAME_REJECTED`: it holds the previous actuator
command and performs prediction-only advancement of the command delay, shaped
plant, modal state, and broadband time history. It does not correct or train an
observer from the duplicated held coordinate. This differs from beam loss,
which deliberately drives command authority to zero during reacquisition.

Requires at least 2 row groups with confident subapertures to evaluate at
all; below that, no rejection can happen (there is nothing to compare). This
means `subap_rows` must be >= 2 for `max_row_shear` to ever do anything.

`max_row_shear` and `row_time` are orthogonal and can be combined: rejection
is evaluated first, so a frame whose shear exceeds the threshold is dropped
rather than passed to the linear correction. Setting a threshold too low
rejects most frames and starves the loop of measurements; too high and
nothing is rejected. Since the shear scales with target velocity, expect the
rejected fraction to rise with vibration amplitude — that is the intended
behavior, but check that enough frames survive to keep the loop closed.

What this does **not** do
--------------------------

- **Rolling shutter beyond the `row_time` correction and `max_row_shear`
  rejection above.** The correction assumes a simple linear
  (constant-velocity-within-frame) skew model; it does not calibrate itself,
  and it does nothing for shutter effects that aren't well described by a
  single per-row time constant (e.g. an electronic rolling shutter with a
  non-uniform per-row exposure schedule). Both mechanisms also assume the
  standard top-row-first readout — a bottom-up or column-wise rolling shutter
  is not handled and would need the row/column axis convention changed.
  Neither can see shear finer than one subaperture row group, since shift
  estimates only exist at row-group granularity.
- **Beam chopping.** This device only prevents a chop-off frame from being
  misread as a valid (near-zero) position -- it falls into the same
  hold-last-value / eventual-reacquire path as any other low-signal frame.
  It does **not** perform synchronous/lock-in demodulation against a chop
  reference; that would need an external chop-state input this device
  doesn't have.
- **Higher-order wavefront modes.** The confidence-weighted average recovers
  common-mode (rigid) tip/tilt robustly, which is what a mirror actuator can
  correct. It does not decompose disagreement between subapertures into
  higher-order Zernike modes -- it only uses that disagreement to down-weight
  the subapertures producing it.

Parameters
----------

- `subap_height`, `subap_width` (required): pixel size of one subaperture.
- `subap_rows`, `subap_cols` (required): subaperture grid shape. The tracked
  window is `subap_rows*subap_height + 2*search_radius` by
  `subap_cols*subap_width + 2*search_radius` pixels.
- `threshold` (default 0): subtracted from each pixel (clamped to zero)
  before use, same as `center_of_mass`.
- `search_radius` (default 2, max 6): integer-pixel correlation search range
  per subaperture per frame. If a subaperture's true motion regularly exceeds
  this, the sub-pixel refinement is skipped for that frame (the peak sits on
  the search boundary) -- widen this if that shows up often.
- `reject_edge_matches` (default false): discard a subaperture match when its
  NCC peak lies on the search boundary. Such a peak is a censored estimate—it
  says the displacement is at least the search radius, not that it equals the
  radius. Enable this in steering loops to prevent a distant artifact from
  walking the window by one search-radius step per bad frame.
- `ref_beta` (default 0.02): EWMA rate for the reference-template update.
  Smaller = templates adapt to real slow changes (focus drift) more
  cautiously; larger = faster adaptation but more sensitivity to letting a
  borderline-confident frame drift the template.
- `min_confidence` (default 0.5): normalized-correlation floor for a
  subaperture to be included in the frame's reconstruction and to update its
  template. Note what this can and cannot catch: it rejects a subaperture
  whose content no longer *resembles* its template, but not one whose
  correlation surface is *flat*. A smooth monotone gradient correlates at
  ~1 against a shifted copy of itself at every offset, so an outer, weakly
  illuminated cell scores near-perfect confidence while its peak position is
  noise. Those subapertures still enter the weighted average (down-weighted
  only by their low flux); what protects the tracker from them is the
  consensus registration described in the algorithm section, not this
  threshold. Raising `flux_floor` so that barely-lit cells never seed at all
  is the effective lever if they are biasing the reconstruction.
- `flux_floor` (default 1.0): minimum thresholded flux in a subaperture's core
  cell when seeding, or in the actual best-matched candidate once initialized.
- `min_valid_subaps` (default ~25% of the grid): confident subapertures
  required for a frame to be treated as valid signal, rather than held/lost.
- `rolling_shutter` (default false): enables centre-row rolling-shutter
  reconstruction in row-index space. This is the preferred interface.
- `row_time` (legacy alias, default 0): a positive value also enables the
  correction, but its magnitude does not affect the result. It remains for
  existing configurations and the diagnostic sweep; new configurations
  should use `rolling_shutter` so they do not imply a false time calibration.
- `max_row_shear` (default 0, disabled): maximum tolerated spread (in pixels)
  between row groups' shift estimates before the whole frame is rejected as
  rolling-shutter distorted and the last valid output is held. Independent of
  `row_time`; see "Rolling-shutter rejection" above. Needs `subap_rows >= 2`.
- `init_y`, `init_x` (optional): fixed acquisition reference and window centre
  in image pixels, used on every reacquisition. If omitted, the window is
  acquired from the brightest pixel and its integer location is refined by a
  local centroid. Supplying the coordinates avoids a brighter reflection
  becoming the absolute steering reference.

  These coordinates define **zero error**, not merely where to look, and that
  distinction matters when comparing against `center_of_mass`. When they are
  given, the acquisition frame reports the beam as sitting exactly at
  (`init_y`, `init_x`) whatever its true position inside the window, so a
  static pointing offset up to about the search margin is absorbed into the
  reference templates instead of being reported; a loop closed on that output
  parks the beam at the configured pixel. `center_of_mass`'s track mode, and
  this device when `init_y`/`init_x` are omitted, instead anchor to a measured
  centroid and so report a real static offset. Two configs that differ only in
  whether these are set therefore do not emit the same absolute coordinate,
  even though their gain and phase response is the same.
- `reacquire_after` (default 10): consecutive low-confidence frames before
  `AYLP_BEAM_LOST` is asserted and the window/templates reset.
