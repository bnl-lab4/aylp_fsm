anyloop:gain_test
=================

Types and units: `[T_VECTOR, U_ANY] -> [T_VECTOR, U_MINMAX]`.

Measures the plant's DC gain with a staircase. The device sits where a
controller (`anyloop:pid`, `anyloop:fsp`) normally sits — sensor upstream, DAC
stage downstream — and walks one output element from `low` to `high` in `step`
increments, holding each level for `dwell` seconds. The trailing `settle_frac`
of every level is averaged into one settled response, and a least-squares line
through the (command, settled response) pairs is the gain.

That gain is the same quantity `bode_plot` reports as `fit_K` and
`latency_test` as `gain`, measured from the same point in the loop — but at DC,
in about half a minute, and with the whole command range visible at once
instead of one small-signal amplitude. Use it to re-check a `K` between sweeps,
or to find where the actuator stops being linear.

What the staircase buys over a two-level step:

- **compression is visible**, as curvature in the per-level residuals. The FSM
  on this rig loses a few percent of gain above ~50 mV, and the `_K_warning` in
  `contrib/steering/configurations/support/attenuation_par_fsp.json` is a long note about how that bit an
  earlier run. The largest residual is reported as a fraction of the response
  span, and the residual panel of the PDF shows its shape.
- **`fit_range`** re-fits the small-signal region the closed loop actually
  lives in. If that slope and the full-range slope differ, the full-range one
  is not the number to put in a controller.
- **`updown`** walks back down again, which measures hysteresis — and is the
  only way to separate gain from beam drift (see below).

Stages: park at `bias` for `warmup` s so the sensor can acquire; glide to `low`
over `ramp` s; run the staircase; glide back to `bias`; then fit, log, write
the files and raise `AYLP_DONE`. Both glides are deliberate — a tracked
`center_of_mass` window only follows motion that stays inside it, so jumping
straight to `low` can lose the beam before the test starts, and jumping home at
the end throws the actuator for no reason.

Sizing the sweep to the sensor
------------------------------

**The useful range is set by the camera, not by the DAC**, and on a
whole-frame estimator it is set by the *beam*, not by the frame:

    span:  |V| < (ROI/2 - 3*sigma) / (K * pixel_scale)
    step:  step * K * pixel_scale < (window/2 - beam radius)

A tail that leaves the frame at the sweep extremes pulls the centroid inward
and reads as compression — which is exactly the curvature the sweep is there to
measure — so the binding limit is the beam's 3σ, not the ROI edge.

Worked from this rig as of 2026-08-07 (ASI290MM, 64×64, 31.5 px per unit,
σ 9.33 px, Kx 0.579 and Ky 0.8156 units/V, i.e. 18.2 and 25.7 px/V): the
ceiling is (31.5 − 28.0)/18.2 = 0.19 V on x and 0.14 V on y. The shipped
configs use ±0.12 V in 12 mV steps — ±2.2 px on x, ±3.1 px on y, a step of
0.22–0.31 px. **That step is at or below the per-sample centroid noise
(~0.27 px); what resolves it is the ~1150 settled samples per level, which put
the SEM near 0.008 px.** Read the per-level SEM in the report. Re-derive all of
it if the ROI, the optics, the beam width or `K` change — and note that `K` is
normalized, so it rescales with `(dim - 1)/2` on top of any optical change.

Levels that break those bounds anyway are recorded but excluded from the fit,
rather than quietly flattening the slope:

- a settled `|response|` past `resp_max` — the beam is running out of the
  sensor. Note this is a *runaway* guard, not the sizing limit above: the
  default 0.9 is 28 px at a 64×64 ROI, far past where a wide beam's tails
  start leaving the frame;
- a settled window whose samples are all bit-identical — a tracked
  `center_of_mass` past `min_peak`, or `fit_com` on a failed fit, holds its
  last coordinate, which otherwise reads as a beautifully quiet level. Only a
  level held *throughout* is caught; see the note on partly-held levels in
  `devices/gain_test.c`.

Two unusable levels in a row end the sweep early *once it has been in range* —
a sweep that starts wider than the sensor walks into range and its good middle
is kept, while one that leaves range at the top stops instead of driving the
beam further off.

Drift, and why it needs `updown`
--------------------------------

The beam on this rig walks several px per hour and the small-signal gain itself
wanders a few percent over one, so a slope measured over 30 s carries whatever
drift happened during those 30 s. The device fits a linear-in-time term
alongside the gain to separate them — but **that fit is only possible on an
up-and-down sweep**. In a plain rising staircase the command and the clock are
the same axis: gain and drift are perfectly collinear and no fit can tell one
from the other. The device checks for that and says so rather than inventing a
split. With `updown`, the down branch repeats every command at a later time,
which breaks the tie; the gap between the plain and the drift-carrying slope is
then the amount of drift that had been landing in the gain.

Output
------

Everything is logged. `output_file` also writes a PDF (response vs command with
the fit, plus a residual panel in px) and, alongside it, a `.dat` carrying the
fit results in its header and one row per level: command in units and in volts,
settled mean, standard deviation, standard error, sample count, timestamp,
residual, validity and branch. Set `output_file` to `""` to log the numbers and
write nothing.

Parameters
----------

- `index_cmd` (int) (optional)
  - Output element to drive. Default 1 (x, for a `[y, x]` command vector).
- `index_err` (int) (optional)
  - Input element to watch. Default 1 (x, for a `center_of_mass` `[y, x]`).
- `out_size` (int) (optional)
  - Output vector length; every other element is held at 0. Default 2.
- `low`, `high`, `step` (float) (optional)
  - Staircase endpoints and increment, in command units. Defaults -0.1, 0.1,
    0.01. `(high - low)` need not be a whole number of steps; the last level is
    clamped to `high`. Note these are *command* units — with a negative DAC
    `scale` the sweep descends in volts.
- `bias` (float) (optional)
  - Command held before and after the sweep, and the centre `fit_range` is
    measured about. Default 0.
- `dwell` (float) (optional)
  - Seconds at each level. Default 1.0.
- `settle_frac` (float) (optional)
  - Trailing fraction of each dwell that is averaged. Default 0.5. The leading
    part is thrown away, so it has to cover the plant's settling.
- `warmup` (float) (optional)
  - Seconds parked at `bias` before the sweep. Default 5. Must exceed the
    sensor's own acquisition time.
- `ramp` (float) (optional)
  - Seconds to glide into and out of the staircase. Default 2.
- `updown` (bool) (optional)
  - Sweep back down after reaching `high`, measuring each level twice and
    reporting the mean up/down offset as hysteresis. Doubles the run. Default
    false. Required for the drift fit.
- `resp_max` (float) (optional)
  - A settled `|response|` past this marks the level unusable. Default 0.9,
    i.e. the very edge of a `U_MINMAX` sensor.
- `fit_range` (float) (optional)
  - Also fit the levels within this distance of `bias`, for a small-signal
    gain. Default 0 (off).
- `volts_per_unit` (float) (optional)
  - The DAC `scale` of the driven channel, used only for reporting: it converts
    the fit into units/V, px/V and mV/px, and carries the sign of the optical
    path. Default 1.
- `pixel_scale` (float | `"auto"`) (optional)
  - Pixels per response unit, `(dim - 1)/2` for a `center_of_mass` output —
    31.5 at a 64×64 ROI. Default 1. Same convention as `attenuation_test`.
    `"auto"` takes it from the frame the source publishes (see
    `libaylp/timing.h`), using the height for `index_err` 0 and the width for 1.
    Use it whenever the source sizes its own ROI — an `asi_source` with an auto
    ROI settles its frame size at startup, so no number written here can be
    trusted to match it. Reporting only: it scales the px columns of the `.dat`
    and the plot, never the fitted gain.
- `output_file` (string) (optional)
  - PDF path; the `.dat` is written alongside it. Default `"gain_test.pdf"`.
    `""` writes nothing.
- `write_config` (string) (optional)
  - Path to a run config whose `anyloop:fsp` stage should receive the measured
    gain, so a gain run leaves the controller ready instead of leaving a number
    to be transcribed. Writes `|small-signal slope|` into the axis's `K` — the
    closed loop lives near the bias, and `K` is positive by convention here with
    the sign carried by the DAC stage's `scale`. Also leaves an
    `_auto_gain_write` note recording the fit and the previous value. Written
    through a temporary file and renamed, so an interrupted write cannot leave a
    half-written controller behind.
  - Refuses to write, and says why, unless: there is a fit; the small-signal fit
    covers ≥ 5 levels; R² ≥ 0.90; the small-signal uncertainty is within 15%;
    and small-signal and full-span gains agree within 25%. The `.dat` still has
    everything either way. Same gates the calibration suite applies.
- `write_axis` (`"x"` | `"y"`) (optional)
  - Which axis of the target config `write_config` updates. Defaults to
    `index_err` (0 → `y`, 1 → `x`).
- `config` (string) (optional)
  - Free text copied into the `.dat` header — the operating point, biases,
    camera settings, why the run was made.

Example
-------

See `contrib/calibration-scripts/configurations/conf_gain_par_x.json` and `contrib/calibration-scripts/configurations/conf_gain_par_y.json` for
complete pipelines (`asi_source → center_of_mass → gain_test → parport_dac`).

```json
{
 "uri": "anyloop:gain_test",
 "params": {
  "index_cmd": 1, "index_err": 1, "out_size": 2,
  "low": -0.12, "high": 0.12, "step": 0.012,
  "dwell": 1.0, "settle_frac": 0.5,
  "fit_range": 0.05,
  "volts_per_unit": 1.0, "pixel_scale": 31.5,
  "output_file": "fsm_gain_x.pdf"
 }
}
```

Tests
-----

`devices/gain_test_test.c` drives the device against a simulated plant of known
gain and checks that the reported slope is that gain — including under
compression, sensor saturation and beam drift. `ninja -C build test` runs it
with no hardware.
