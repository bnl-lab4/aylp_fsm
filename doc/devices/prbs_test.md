anyloop:prbs_test
=================

Types and units: `[T_VECTOR, U_ANY] -> [T_VECTOR, U_MINMAX]`.

Measures the loop's end-to-end delay by driving one output element with bursts
of a maximum-length pseudorandom binary sequence and cross-correlating the
command that went out against the response that came back. The device sits
where a controller normally sits — sensor upstream, DAC stage downstream — so
the measured delay covers the whole physical path: DAC write, actuator dead
time and motion, exposure, camera readout, centroid. Exactly the delay a
controller in the same spot sees, and the number `fsp` wants as
`delay`/`delay_frac`.

It answers the same question as `latency_test` and as `bode_plot`'s
`fit_tau_ms`, and is the cheapest of the three: a full measurement is a few
seconds of drive, against a Bode sweep's several minutes. Where a square-wave
step spends its energy on one edge at a time, a PRBS spreads it flat across the
band and correlation gathers all of it, so the delay comes out of a few
thousand frames instead of a few hundred edges.

Why bursts
----------

Every burst is driven from the same LFSR seed, so all of them are the *same*
stimulus and their responses can be ensemble-averaged before anything else
happens. Room tone, the rig's 5–80 Hz disturbance lines and centroid noise are
uncorrelated with the sequence, so they fall as `1/sqrt(n_bursts)` while the
response does not. The quiet tail between bursts lets the actuator return to
rest, so each burst starts from the same state, and it gives the correlogram
somewhere to show its own noise floor.

The three delays
----------------

They are three different things, and the gap between them is the plant's rise
rather than disagreement:

- **onset** — the lag where the correlation crosses
  `max(onset_frac × peak, 4 σ)`, interpolated between the two straddling lags.
  The earliest of the three, and it reads **early**, not late: an m-sequence's
  autocorrelation is a triangle a chip wide either side, so the correlogram
  starts rising one chip before the true lag. Measured against a pure
  integer-frame delay it sits a constant **0.79 chips early** at the default
  `onset_frac` of 0.2. Treat it as an early-side bound on first motion, not as
  a number to quote.
- **peak** — the lag of the largest correlation, refined to sub-frame by
  parabolic interpolation. Against a pure delay this is **exact**: measured
  within 0.003 frames for delays of 2–6 frames. Against a real plant it is the
  delay plus the actuator's rise.
- **phase slope** — the group delay from a magnitude-weighted straight-line fit
  to the phase of the cross-spectrum over `[phase_f_lo, phase_f_hi]`. Also
  exact against a pure delay (within 0.0005 frames). This is the same
  estimator, over the same band, that `bode_plot` fits, so it is the number to
  compare against a sweep and the one to hand to `fsp`.

Peak and phase slope agree on a pure delay and separate on a real plant, and
that separation is the actuator's rise rather than an error in either.

Comparing against a Bode sweep: mind the hidden frame
-----------------------------------------------------

**Compare in frames against `fit_frames`, never in ms against `fit_tau_ms`.**
The two devices reference their correlators differently, and the difference is
exactly one frame.

`bode_plot` correlates the response `y[n]` against `prev_phi`, the excitation
phase of iteration *n−1* (`bode_plot.c`, step 5), because `y[n]` cannot be a
response to the `u[n]` that goes out downstream later in the same iteration.
Its fitted phase is therefore referenced to a stimulus that is already one
frame old, and `fit_tau_ms` comes out one frame *short* of the real
command-to-response interval; `fit_frames` adds the `+1` back at the measured
rate. It is bookkeeping, not physics.

This device has nothing to add back. It correlates `resp[k+L]` against the
command actually emitted at frame `k`, so a reported lag is the whole interval
from the DAC write to the sensor value arriving. Verified against a plant that
is a pure N-frame delay: the reported lag is N, to within 0.003 frames, for
N of 2 through 6.

So on this rig, x reads `fit_tau_ms` 0.4654 ms and `fit_frames` 2.763 — and it
is 2.763 that this device should reproduce (about 0.729 ms at 3788 Hz), not
0.4654 ms. `fsp`'s `delay` + `delay_frac` is in that same hidden-frame-included
convention (2 + 0.76), so the phase-slope figure here goes straight into it.

A **centroid** of the correlation lobe is reported too, and every burst is
correlated on its own as well: the median absolute deviation over bursts is the
honest uncertainty on the ensemble numbers.

Reading a lag
-------------

A lag of L frames means the response appeared L iterations after the command
was emitted. **L = 0 is physically impossible** — the sensor value read in an
iteration was captured before that iteration's command went out — so a reported
lag below one frame means the correlation found noise, and the run should be
treated as a failure rather than a fast loop. The device warns when that
happens.

There is no hidden frame to add back, unlike `bode_plot` — see the section
above for what that is and how to compare the two.

Sizing the drive
----------------

- `amplitude` must keep the beam on the sensor. The plant reaches roughly full
  deflection during the sequence's longest run of like chips (`order` chips),
  so budget `amplitude × K × pixel_scale` px of excursion — at this rig's
  2026-08-07 Kx 0.579 units/V and a 64×64 ROI, 0.05 V is about ±0.9 px.
  Averaging is what buys SNR here, not amplitude; the shipped configs stay well
  inside the linear region on purpose, so the delay is not measured at an
  operating point the loop never visits. The ceiling is the beam rather than
  the frame — keep 3σ inside the ROI, since truncation at the extremes is
  *synchronous with the drive* and so cannot be averaged out.
- `chip_frames` should stay at 1 unless the response is too weak to see. A
  chip is a boxcar, and a chip longer than the delay smears the onset by its
  own length: at 2.3 kHz the delay is about three frames, so a 10-frame chip
  would measure almost nothing else.
- **A chopped or dropping source is handled, not assumed away.** See below.
- `quiet_frames` should be at least `max_lag`, or one burst's tail overlaps the
  next burst's start. The device warns if it isn't.
- The window (`burst + quiet`) has to be longer than the lags asked about,
  since every lag is evaluated over the same stretch of command; init refuses
  a combination that isn't.

The correlation's noise floor is not just sensor noise: the *linear* (rather
than circular) autocorrelation of a finite m-sequence has sidelobes of order
`1/sqrt(n_chips)`, a few percent of the peak. That is what the acausal lags
measure and what the onset threshold is set against.

A chopped source, and frames the sensor did not measure
-------------------------------------------------------

A centroid stage that cannot fit a frame does not report a gap — `fit_com`,
`wfs_com` and a tracked `center_of_mass` all **re-publish their previous
coordinate** and raise `AYLP_FRAME_REJECTED`. Those samples are not missing
data, they are wrong data: each carries some earlier frame's value.

This device records that flag for every frame and drops the held samples from
the ensemble mean and from every correlation. Nothing is assumed about the
duty cycle, its period, or its regularity; the live fraction is *measured* and
reported. What makes the repair cheap is that a chop is almost never locked to
the burst period, so a window position that was dark in one burst is live in
others and the ensemble mean fills itself in.

Measured on the simulated plant in `devices/prbs_test_test.c` (a 4-frame lag,
25% of frames held in 13-frame streaks), masking recovers the CW answer to
within 0.001 frames on the peak and 0.001 on rho. Correlating the held samples
instead — `use_rejected: true` — costs:

| held | streaks | peak | lobe centroid | phase slope | rho |
|------|---------|------|---------------|-------------|------|
| 0%   | —       | 4.165 | 4.703 | 4.521 | 0.866 |
| 25%  | 13 fr   | +0.013 | +0.028 | −0.24 | 0.836 |
| 51%  | 26 fr   | +0.031 | +0.293 | −0.66 | 0.732 |
| 10%  | 5 fr    | +0.015 | +0.089 | ±0.03 | 0.838 |

So the damage lands on the **phase slope** — the number you copy into an fsp
`delay`/`delay_frac` pair — and it reads *short*, not long. Held samples have a
spread of ages, so they decorrelate the response rather than shifting it
coherently, and a decorrelated high-frequency end flattens a phase-slope fit.
Below about 10% held it is all decorrelation and the sign is not even stable.
The exact size depends on where `[phase_f_lo, phase_f_hi]` sits in the
spectrum.

Two things to read in the log before quoting a delay:

- **the live fraction.** Under ~50% the ensemble mean is thin and the noise
  floor rises correspondingly. That is a source problem, not a fit problem.
- **holes.** A window position dark in *every* burst has no measurement at all.
  The correlogram skips it, but the phase slope is an FFT and cannot, so it is
  filled by interpolating its neighbours and the run warns. This is what a chop
  locked to the burst period does; change `quiet_frames` to walk the two apart.
  `resp_n` in the `_traces.dat` file is the live burst count per position, and
  0 marks an interpolated one.

Quality gates
-------------

Finding a correlation lobe is not, by itself, enough to configure a
controller. Calibration configurations can require a minimum peak/noise ratio,
a sufficiently linear phase fit, enough independent bursts, low burst-to-burst
peak scatter, adequate live-frame coverage, no interpolated holes, and
reasonable agreement between the correlation peak and phase delay. A failed
gate leaves the candidate diagnostics in the `.dat` file but suppresses the
PDF and prints `QUALITY GATES FAILED`; the delay must not be copied into FSP.

The shipped parport configurations also follow the gain ladder's central
cross-check: the calibration suite repeats each axis at 35, 50 and 75 mV and
requires the phase-delay estimates to agree within 0.2 frame. A tracker or
plant whose dynamics change with excursion therefore fails explicitly instead
of contributing one amplitude-dependent delay.

`use_rejected` exists so the two can be run against the same data and compared,
which is the only way to see what a given source is costing. It is not a way to
get a number out of a bad run.

Loop-rate jitter
----------------

The stimulus is emitted per *iteration* and the correlation is over frame
index, so an irregular loop period does not smear a delay that is itself
frame-quantized (exposure, readout, transport). The unfavourable case is a
delay fixed in **time** sampled at irregular instants, where a constant delay
maps to a varying number of frames. Simulated at that unfavourable case — a
fixed-time delay of 3.70 frames at 3788 Hz through a 0.8-frame rise, 64 bursts,
uniform jitter of the stated rms as a fraction of the period:

| jitter rms | onset | peak | phase slope | peak ρ | per-burst MAD |
|---|---|---|---|---|---|
| 5 % | 3.294 | 4.082 | 4.374 | 0.94 | 0.002 fr |
| 17 % | 3.331 | 4.129 | 4.483 | 0.89 | 0.022 fr |
| 30 % | 3.297 | 4.222 | 4.556 | 0.78 | 0.081 fr |
| 60 % | 2.802 | 4.150 | 5.255 | 0.63 | 0.143 fr |

So it survives comfortably to ~30 % rms — the phase slope moves 4 %, onset
moves 0.003 frames, and every burst still yields its own estimate. (The
absolute offsets in that table are the plant's rise and the onset bias
described above, not jitter; what jitter costs is the movement *down* a
column.) Past that it
degrades honestly rather than silently: at 60 % the correlation peak has lost a
third of its height, the leading edge smears *earlier* (onset falls below its
true value), the phase slope runs 20 % long, and only 37 of 64 bursts produce
an estimate at all. The measured jitter is reported in the log and the `.dat`
header, and the device warns above 20 %; the per-burst MAD and the peak/noise
ratio are the two numbers to read when deciding whether a run is trustworthy.

Output
------

Everything is logged. `output_file` also writes a PDF (the stimulus with the
ensemble-averaged response above, the correlogram with onset and peak marked
below) plus two `.dat` files alongside it: `<base>.dat` with the fit results in
its header and one row per lag, and `<base>_traces.dat` with the
ensemble-averaged command and response frame by frame. Set `output_file` to
`""` to log the numbers and write nothing.

Parameters
----------

- `index_cmd` (int) (optional)
  - Output element to drive. Default 1 (x, for a `[y, x]` command vector).
- `index_err` (int) (optional)
  - Input element to watch. Default 1 (x, for a `center_of_mass` `[y, x]`).
- `out_size` (int) (optional)
  - Output vector length; every other element is held at 0. Default 2.
- `amplitude` (float) (optional)
  - PRBS half-swing in command units, about `bias`. Default 0.05.
- `bias` (float) (optional)
  - Command between bursts, and before and after the test. Default 0.
- `order` (int, 5–16) (optional)
  - LFSR order; a burst is `2^order - 1` chips. Default 7 (127). The tap masks
    are checked at init by walking the whole cycle, so a short sequence — which
    would put a periodic ridge in the correlogram — cannot get through.
- `chip_frames` (int) (optional)
  - Loop iterations per chip. Default 1.
- `n_bursts` (int) (optional)
  - Bursts to average. Default 32. Noise falls as its square root.
- `quiet_frames` (int) (optional)
  - Frames parked at `bias` after each burst. Default 128.
- `warmup` (float) (optional)
  - Seconds parked at `bias` before the first burst, for sensor acquisition.
    Default 5.
- `max_lag`, `neg_lags` (int) (optional)
  - Positive lags evaluated, and acausal ones. Defaults 48 and 24. The acausal
    lags are the noise floor; keep enough of them to measure it.
- `onset_frac` (float) (optional)
  - Fraction of the peak that counts as onset. Default 0.2.
- `phase_f_lo`, `phase_f_hi` (float) (optional)
  - Band for the phase-slope fit, Hz. Defaults 20 and 200 — the top matches
    the `freq_end` of the Bode sweeps on this rig, so the two tau figures are
    fit over the same band. Pushing the top past the plant's own roll-off just
    fits noise.
- `volts_per_unit` (float) (optional)
  - The DAC `scale` of the driven channel, used only for reporting. Default 1.
- `pixel_scale` (float) (optional)
  - Pixels per response unit, `(dim - 1)/2` for a `center_of_mass` output.
    Default 1.
- `use_rejected` (bool) (optional)
  - Correlate frames the sensor flagged `AYLP_FRAME_REJECTED` as if they were
    measurements. Default false, which drops them. For comparing a run against
    itself to see what a chopped source costs — see the section above — not for
    producing a number. The run warns when it is set.
- `min_pair_frac` (float) (optional)
  - Fraction of the correlation window a lag needs in live sample pairs before
    it is fit at all; lags below it are dropped from the correlogram rather
    than allowed to win the peak on a handful of samples. Default 0.25, floored
    at 8 pairs. Only bites when the live fraction is very low or a chop is
    locked to the burst period.
- `min_peak_snr` (float) (optional)
  - Required absolute correlation peak divided by the acausal-lag RMS. Zero
    disables the final-result gate; peak detection itself always requires 4.
- `max_phase_resid_deg` (float) (optional)
  - Maximum RMS residual of the phase-slope line, in degrees. Zero disables.
- `min_phase_bins` (int) (optional)
  - Minimum number of FFT bins in the phase fit. Zero disables.
- `min_burst_frac` (float) (optional)
  - Minimum fraction of individual bursts that must independently find a
    correlation peak. Zero disables.
- `max_peak_mad` (float) (optional)
  - Maximum median absolute deviation of the per-burst peak lags, in frames.
    Zero disables.
- `min_live_frac` (float) (optional)
  - Minimum fraction of sensor frames not marked rejected. Zero disables.
- `max_holes` (int) (optional)
  - Maximum window positions with no live sample in any burst. Unlimited by
    default; strict calibration configurations set zero.
- `max_phase_peak_delta` (float) (optional)
  - Maximum absolute difference between phase-slope delay and correlation
    peak, in frames. Zero disables. Leave enough room for real plant rise.
- `output_file` (string) (optional)
  - PDF path; the `.dat` files are written alongside it. Default
    `"prbs_test.pdf"`. `""` writes nothing.
- `config` (string) (optional)
  - Free text copied into the `.dat` header.

Example
-------

See `contrib/calibration-scripts/configurations/conf_prbs_par_x.json` and `contrib/calibration-scripts/configurations/conf_prbs_par_y.json` for
complete pipelines (`asi_source → center_of_mass → prbs_test → parport_dac`).

```json
{
 "uri": "anyloop:prbs_test",
 "params": {
  "index_cmd": 1, "index_err": 1, "out_size": 2,
  "amplitude": 0.05, "order": 7, "chip_frames": 1,
  "n_bursts": 64, "quiet_frames": 128,
  "phase_f_lo": 20, "phase_f_hi": 200,
  "volts_per_unit": 1.0, "pixel_scale": 31.5,
  "output_file": "fsm_prbs_x.pdf"
 }
}
```

Tests
-----

`devices/prbs_test_test.c` drives the device against a simulated plant of known
delay and checks that the reported onset, peak and phase-slope numbers are that
delay — including on an inverting plant, under noise several times the
response, and on a plant that never moves, which must be reported as a failure
rather than as a fast loop. `ninja -C build test` runs it with no hardware.
