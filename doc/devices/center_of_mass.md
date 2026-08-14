anyloop:center_of_mass
======================

Types and units: `[T_MATRIX_UCHAR, U_ANY] -> [T_VECTOR, U_MINMAX]`.

This device breaks up an image into one or more regions, and calculates the
center-of-mass coordinate of that image. For example, this device might be used
with only one region to determine the center-of-mass coordinates of a beam on a
camera, which can then be used to control a tip-tilt mirror to recenter said
beam. This device is also used with many regions for getting error signals from
a wavefront sensor.

An example configuration for a wavefront sensor with 8x8-pixel subapertures:

```json
{
  "uri": "anyloop:center_of_mass",
    "params": {
      "region_height": 8,
      "region_width": 8,
      "thread_count": 1
    }
}
```

Pipeline data is replaced with a vector of interleaved center-of-mass y and x
coordinates (a vector of length 2N, where N is the number of regions of
interest). For example, if the input has four regions of interest, the output
will be [y1,x1,y2,x2,y3,x3,y4,x4] where each y,x is from -1 to 1, where 0 means
perfectly centered in the region of interest. It is assumed that the input is
written in order of increasing x coordinate, then increasing y coordinate.

Parameters
----------

- `region_height` (integer or `"auto"`) (required)
  - Height of each region to find the center of mass of. The image will be split
    up into regions of this height, from the top going down. Excess data will be
    ignored. Set this to 0 to set the region height to the logical height of the
    whole image.
  - In tracked mode, `"auto"` selects the first frame's height, capped at 243
    pixels to retain the registration timing budget.
- `region_width` (integer or `"auto"`) (required)
  - Width of each region to find the center of mass of. The image will be split
    up into regions of this width, from left to right. Excess data will be
    ignored. Set this to 0 to set the region width to the logical height of the
    whole image.
  - In tracked mode, `"auto"` selects the first frame's width, capped at 243
    pixels to retain the registration timing budget.
- `thread_count` (integer) (optional)
  - Number of threads to use for the calculation. Set this to 1 (default) for no
    multithreading. Ignored when `track` is set.
- `threshold` (integer 0-255, or `"auto"`) (optional)
  - Subtracted from each pixel before the centroid is formed. `"auto"` takes
    the 95th percentile of the tracking-window border on each frame, avoiding a
    fixed camera-noise setting. Numeric values retain the original fixed-gate
    behavior.
- `track` (boolean) (optional)
  - Confine the sum to a single `region_height` by `region_width` window centred
    on the previous frame's center of mass. Defaults to false.
- `registration` (boolean) (optional)
  - In tracking mode, report rigid image translation against a valid keyframe
    instead of its flux centroid. The original keyframe is retained while enough
    of its two-axis structure remains visible, avoiding dead-reckoning error from
    unnecessary keyframe chains. If overlap is exhausted, a new keyframe is
    installed only from an independently measured, unclipped flux centre.
    Border-clipped centroids are rejected because lost flux makes them move even
    under a perfectly rigid translation. A rolling median of valid centroid
    disagreement supplies a bounded slow absolute correction. The fit includes per-frame intensity
    gain and background offset, so brightness changes and evolving asymmetric
    illumination do not become fictitious motion. It assumes only that the
    tracked pattern has spatial structure in both axes; it does not assume a
    Gaussian beam, a particular size, or a single spot. Noise is estimated from
    the keyframe, samples are added until both translation axes reach the same
    uncertainty target, and each frame derives its robust residual scale from
    its own median residual. Candidate and scratch storage are derived from the
    configured window dimensions; there is no preferred ROI size or fixed
    sample ceiling. There are deliberately no beam-specific sampling, gradient,
    or outlier parameters. Defaults to false.
- `min_peak` (integer 0-255, or `"auto"`) (optional)
  - The brightest pixel inside the window must reach this for the frame to count
    as holding the beam; on frames that fall short, the last good center of mass
    is held instead of being updated. Defaults to 0, which disables the test.
    Only meaningful with `track`. Distinct from `threshold`, which merely shapes
    the weighting: this decides whether the frame is used at all. Must be set
    above `threshold` to reject anything. `"auto"` starts three counts above
    the border background, then uses 10% of the accepted-frame mean peak above
    that background.
- `ref_cut` (float 0-1) (optional)
  - Reject frames holding only part of the beam. Each frame's row profile is
    divided by a learned reference profile; if the dimmest significant row falls
    below this fraction of the frame's typical row, the last good center of mass
    is held instead. Defaults to 0, which disables the test. Only meaningful with
    `track`, and needs `region_height` of at least 4. See "Partial beams" below.
- `ref_warmup` (integer) (optional)
  - Frames used to learn the reference before the gate switches on. Defaults to
    200.
- `ref_rate` (float, 0 exclusive to 1) (optional)
  - EMA rate at which accepted frames update the reference afterwards. Defaults
    to 0.01.
- `ref_floor` (float, 0 to 1 exclusive) (optional)
  - Rows whose reference is below this fraction of the brightest reference row
    are excluded from the test. Defaults to 0.25.
- `init_y`, `init_x` (integer, or both `"auto"`) (optional)
  - Initial window centre, in image pixels. Must be given together. If omitted,
    or set to `"auto"`, the window is acquired from the brightest pixel of the
    first usable frame.
- `reacquire_after` (integer, or `"auto"`) (optional)
  - Consecutive frames of zero signal inside the window before re-acquiring from
    the brightest pixel of the whole image. Defaults to 30. At that point the
    device also asserts the `AYLP_BEAM_LOST` pipeline-status flag. A valid
    centroid clears the flag so downstream control can resume. `"auto"` starts
    at 30 frames and doubles after each unexpectedly long dark interval, so a
    pulsed source teaches the hold time without embedding its duty cycle in the
    configuration.
- `acquire_seconds` (float) (optional)
  - Run with a wide acquisition window for this long before narrowing to
    `region_height`/`region_width`. Defaults to 0 (no acquisition phase). Also
    re-entered whenever the window re-acquires after losing the beam.
- `acquire_height`, `acquire_width` (integer) (optional)
  - Size of the acquisition window. Default 0, meaning the whole image.

Tracking window
---------------

Center of mass is a flux-weighted average, so anything bright in the region
contributes. A stray reflection alongside the beam does two things: it pulls the
centroid off the beam by a fixed offset, and — the one that matters — it
attenuates the response to real beam motion by `S_beam / (S_beam + S_reflection)`,
silently scaling your loop gain.

Setting `track` confines the sum to a window that follows the beam, so a
reflection outside that window never enters the sum. The image is passed down the
pipeline untouched, so a `udp_sink` placed *before* this device still shows the
whole frame, reflection and all.

```json
{
  "uri": "anyloop:center_of_mass",
    "params": {
      "region_height": 25,
      "region_width": 25,
      "threshold": 20,
      "track": true
    }
}
```

The output is normalized across the **whole image**, not the window. This is what
makes the mode usable in a control loop: the window chases the beam, so a
window-relative coordinate would sit near zero no matter where the beam actually
was, and the loop would have no error signal to act on. Because the normalization
is to the image, the setpoint stays the image centre and the error per pixel of
beam motion — hence the loop gain — does not change with window size. You can
resize the window freely without retuning the PID.

Behaviour worth knowing:

- **Acquisition.** By default the window lands on the brightest pixel of the
  first frame. That scan keeps the *first* maximum it meets in raster order, so
  when the beam and a reflection both saturate at 255 the one nearer the top of
  the frame wins — which may well be the reflection.

  `acquire_seconds` fixes this without hardcoding a position: for that long the
  sum runs over a wide window (the whole image by default), where the centroid is
  flux-weighted across everything present and therefore drawn toward whichever
  spot carries the most *light*. When the phase ends, the window narrows onto
  wherever that centroid settled, and because it re-centres every frame it then
  walks onto the dominant spot.

  This works **iff the beam carries more total flux than the reflection**. Peak
  brightness is irrelevant once both saturate; what matters is the integral above
  `threshold`. At an even flux split the window stalls between the two spots, and
  past that it converges onto the reflection. When the beam does not dominate,
  use `init_y`/`init_x` — that is deterministic and needs no acquisition phase.
- **How close the beam may get.** The window reaches `region_height/2` from its
  centre. Contamination returns once the reflection's above-threshold pixels
  reach inside that, i.e. when the beam-to-reflection separation drops below
  roughly `region/2 + reflection_radius`. Size the window as small as the beam's
  per-frame motion and spot size allow.
- **Loss of signal.** If every pixel in the window is at or below `threshold`,
  the brightest one falls short of `min_peak`, or some rows are lit while others
  are not by more than `ref_cut` allows, the last valid output is held rather
  than reporting `(0, 0)` — which downstream reads as "perfectly centred",
  not "no signal", and would let an integrator park and then lurch when the beam
  reappears. Every held frame carries `AYLP_NO_SIGNAL`, which is what `pid`'s
  `park_after` counts. After `reacquire_after` such frames the window
  additionally asserts `AYLP_BEAM_LOST` and re-acquires from the brightest pixel
  of the whole image, so a stranded window can recover. FSP uses that explicit
  status, rather than a large error or command, to hold output at zero only while
  the beam remains missing. Note that reacquisition can lock onto a reflection if
  the beam is genuinely gone.
- **Why `threshold` alone is not a beam test.** Without `min_peak` the test for
  "is the beam here" is just "is the thresholded sum nonzero", which a dark
  sensor passes trivially whenever its read noise runs a few counts above
  `threshold` — and the resulting centroid is noise, wandering the full width of
  the window. Measured on an ASI662MM with the beam blocked, `threshold` 2 and
  `min_peak` 0: 300 consecutive frames gave 300 distinct centroids spanning
  0.62 px with an RMS of 0.108 px, i.e. a fictitious error signal several times
  larger than a good closed-loop residual. With `min_peak` 40 the same 300 frames
  gave one held value. Whether the old behaviour bites depends on where
  `threshold` sits relative to the noise floor, so it can lie dormant for months
  and then appear after a camera, gain or exposure change.
- **Choosing `min_peak`.** Block the beam and read the window peak the warning
  reports, then leave several times that margin. Too high is not the safe
  direction: a beam dimmer than `min_peak` freezes the centroid, and the loop
  then drives on a stale error, which is worse than a slightly noisy one.
- **Edges.** The window is clamped to lie inside the image, so a beam pinned
  against an edge reports a coordinate near ±1 rather than reading out of bounds.

Partial beams
-------------

A rolling shutter reading out across a chopped or pulsed source gives frames in
which only some rows caught the source's on-window, so the beam is cut by a hard
horizontal edge. Neither brightness test sees this: the rows that were lit are at
full brightness, so `min_peak` passes, and a large fragment still carries plenty
of flux. The centroid of such a frame is pulled toward the surviving rows, which
is a systematic error along the shutter axis — it does not average out, and since
the chop beats against the frame rate it arrives as a slow wander or an
oscillation rather than as visible noise.

`ref_cut` asks the question directly — are some rows lit while others are not? —
by keeping a reference profile of what a whole beam looks like and dividing each
frame's row profile by it. On a whole beam every row comes out at the same ratio,
whatever that ratio happens to be, so the frame's overall brightness *cancels*
rather than merely being normalized away: a dimmer frame is uniformly dimmer. On
a cut beam the ratios split in two — near the frame's own level where the shutter
let light through, near zero where it did not — and `ref_cut` is the gap between
the dimmest row and the frame's typical (median) row that counts as split.

Because the reference is the beam's own measured profile, no assumption about
beam shape or size enters. A tight spot is handled exactly like a broad one,
which a test keyed on the profile's steepness cannot do: a Gaussian one row wide
changes by 0.61 of its peak row between adjacent rows all by itself, which no
steepness threshold can tell from a cut.

The reference is learned in two stages. For the first `ref_warmup` frames it is
the row-wise **maximum**, not the mean, because the source may already be chopping
while the warmup runs — a cut only ever removes light, so each row's largest
value over enough frames is its uncut value, while a mean would learn a blend of
whole and cut beams and sit too low. After that it is an EMA at `ref_rate` over
*accepted* frames only, so it follows slow drift in power, focus and alignment
without letting cut frames pull it down toward themselves. The gate is off during
the warmup, and the reference is thrown away whenever the window re-acquires,
since it describes a beam at the old position.

- **Choosing it.** The separation is large enough that the value is not delicate.
  Measured on a simulated σ = 2 beam: a whole beam reads 1.00, a beam dimmed to a
  quarter of the brightness the reference was learned at still reads 0.61, and
  any cut through the significant rows reads 0.00. 0.5 is a good default; lower
  it toward 0.3 if the source swings in power by more than about 4×.
- **`ref_floor` sets the sensitivity.** Only rows the reference says should carry
  real signal are judged, since elsewhere the ratio is a small number over a
  small number. Lowering the floor extends the test further into the beam's
  skirts and catches shallower cuts, at the cost of noisier ratios and more
  sensitivity to beam motion. Unlike a shape-based test this is a knob rather
  than a hard limit.
- **The one part that is not exactly scale-free.** `threshold` is subtracted
  before the profile is summed, and that subtraction is not linear in pixel
  value, so a much dimmer frame loses proportionally more of its skirts and its
  ratios tilt slightly. That is the whole reason a 4× dimming reads 0.61 instead
  of 1.00. Keeping the beam well above `threshold` keeps the test exact.
- **Fast beam motion.** The window recentres on the previous frame's centroid, so
  the beam sits at a stable place in the window and the profile is stationary
  there. A beam that moves a large fraction of its own width within one frame
  will tilt the ratios and can be held. This is safe — holding is what you would
  want on a frame you cannot trust — but it means the gate is not free during a
  violent transient.
- **Watch the reject rate.** It is reported at exit: `held N of M frames`. Every
  held frame hands the loop the previous frame's error a second time, so the
  integrator keeps winding on stale data — downstream this reads as added loop
  delay and eats phase margin. A few percent is free; a third is not, and at that
  point the fix is to stop generating partial frames (expose over a whole chop
  period, or trigger the camera off the chop) rather than to filter them.
- **Log volume.** The beam-lost and beam-back lines are rate-limited to one per
  second, because a gate that rejects alternate frames would otherwise emit
  hundreds of log lines a second from inside the control loop. Suppressed
  episodes are still counted in the exit summary.

Without `track`, the device behaves exactly as before: the region grid starts at
the top-left of the image and tiles across it, and the output is normalized per
region.
