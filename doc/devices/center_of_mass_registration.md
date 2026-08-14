# Center-of-mass registration tracking

The optional `registration` mode on `anyloop:center_of_mass` estimates rigid
image translation against a retained keyframe. It fits translation together
with an intensity gain and offset, uses keyframe gradients, and applies
per-frame Tukey robust weighting. The ordinary flux centroid remains the slow
absolute anchor only when a keyframe handoff is needed; a clipped border
centroid is not used as an anchor.

## Validation

The implementation is covered by `devices/center_of_mass_test.c`. The replay
mode was run on the 13 camera captures currently present in the repository
(256, 320, and 384 pixel ROIs). With `threshold=1` and `min_peak=100`, the
same 16,198 illuminated frames were accepted with registration enabled and
disabled:

- pooled combined position RMS: **6.23 px with registration**, **15.58 px
  without** (60.0% lower);
- pooled frame-to-frame/HF RMS: **0.090 px with registration**, **1.981 px
  without** (95.5% lower);
- per-capture combined RMS improvement: 13–90% (12 of 13 captures improved by
  at least 30%).

These are replay measurements of recorded beam images, not a claim that every
removed centroid excursion was physically false. The injected-translation
tests provide the absolute check: pooled registration error is about 0.016 px
radial RMS, recovered gain is effectively unity, and terminal bias is
negligible. Forward/backward closure is 0.110/0.092 px RMS, with a worst case
of 0.197 px; three-frame composition is 0.360 px RMS, worst 0.652 px.

## Runtime

The steady-state registration path is below 100 microseconds on the replay
captures: approximately 42–43 us median for 256x256, 53 us for 320x320, and
64 us for 384x384. The 243x243 production window is expected to be about
40–45 us. Installing or reacquiring a keyframe is a deliberately slower
one-time path because it scans and ranks the ROI; it is not representative of
the steady-state loop delay and should be measured separately when fitting
reacquisition behavior.

## Live-use status

The registration implementation passes all eight Meson suites and the full
capture replay. The canonical steering profiles now set `registration=true`
after replay validation. A live A/B run must still measure the complete camera
→ COM → controller delay before closed-loop rejection numbers are compared. Do
not copy the old profile notes that described registration as permanently
disabled due to the earlier local-keyframe failure; that failure was fixed by
retaining the original keyframe, using an unclipped anchor for handoff, and
validating absolute translation and overlap.
