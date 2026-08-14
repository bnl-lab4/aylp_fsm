anyloop:latency_test
====================

Types and units: `[T_VECTOR, U_ANY] -> [T_VECTOR, U_MINMAX]`.

Measures the loop's end-to-end feedback-path latency. The device sits where a
controller (e.g. `anyloop:pid`) normally sits — sensor upstream, DAC stage
downstream — and toggles one output element between `low` and `high` as a
square wave while watching one input element (typically a `center_of_mass`
coordinate) for the response. The command edge is timestamped in the same
iteration that hands the new value to the DAC stage, and the response at the
iteration that first carries the moved centroid, so the measured delay covers
the whole physical path: serial write → DAC → actuator motion → exposure →
camera readout/USB → centroid computation. Exactly the delay a controller in
the same spot would see.

The test runs in stages:

1. **warmup** — hold `low` for `warmup` seconds so the sensor side (e.g.
   `center_of_mass` acquisition) can settle.
2. **cal** — run `cal_cycles` low/high cycles. The last `settle_frac` of each
   half-cycle is averaged into a settled window mean. The step size is
   estimated from consecutive high/low window-mean differences, which cancels
   slow open-loop beam drift, and the noise floor from within-window scatter
   only, so that same drift doesn't count as noise. The test refuses to
   continue if the step is < 8 sigma of the noise — raise `|high - low|` if it
   does.
3. **measure** — for `n_steps` edges, two latencies are recorded per edge,
   both relative to a baseline re-measured from the settled tail of the
   half-cycle the edge departs from (so drift between edges cannot stale it):
   - *departure*: first of 2 consecutive samples past
     `max(4 sigma, 2% of step)` away from the baseline toward the target —
     transport delay plus first motion;
   - *50% crossing*: first sample past baseline + step/2 — delay plus half
     the rise.
   Every edge is logged as it is measured, so a wedged test is visible
   immediately.
4. **report** — log median/mean/min/max of both statistics, the mean sample
   interval, and the median departure expressed in frames; then park the
   command vector at zero and set `AYLP_DONE`.

On any fatal condition (bad input size, too few settled samples, step lost in
noise) the device parks the output at zero, still publishes it, and sets
`AYLP_DONE` itself, since the main loop treats proc errors as recoverable.

See `contrib/legacy/outdated_scripts/conf_latency_test.json` for a complete pipeline
(`asi_source → center_of_mass → latency_test → piplate_bridge`).

Parameters
----------

- `index_cmd` (int) (optional)
  - Output element to drive. Default 0.
- `index_err` (int) (optional)
  - Input element to watch. Default 1 (x, for a `center_of_mass` `[y, x]`
    vector).
- `out_size` (int) (optional)
  - Output vector length; elements other than `index_cmd` are held at 0.
    Default 4 (the dual-stage DAC map).
- `low`, `high` (float) (optional)
  - Command levels in minmax units. Defaults -0.05 and 0.05. The resulting
    beam step must stay well inside the `center_of_mass` window, but must
    also clear the 8-sigma calibration gate; if calibration reports
    `step/sigma < 8`, raise the swing.
- `period` (float) (optional)
  - Seconds spent at each level (one half-cycle). Default 0.25. Must be long
    enough for the actuator and centroid to settle well within it.
- `warmup` (float) (optional)
  - Seconds to hold `low` before calibrating. Default 6. Must exceed
    `center_of_mass`'s `acquire_seconds`.
- `cal_cycles` (int) (optional)
  - Full low/high cycles used for calibration. Default 3.
- `n_steps` (int) (optional)
  - Command edges to measure before reporting. Default 40.
- `settle_frac` (float) (optional)
  - Fraction at the end of each half-cycle treated as settled, used for the
    calibration windows and the per-edge baselines. Default 0.4.
- `results_file` (string) (optional)
  - Append the full results to this file: a run header with a wall-clock
    timestamp and the test parameters, the calibration line (step, sigma),
    one line per edge (`edge direction departure_ms half_ms`, −1 for a
    response that was never detected), and the summary statistics as `#`
    comment lines (median/mean/min/max of both latencies, sample interval,
    loop rate, median departure in frames). The file is opened in append
    mode and flushed after every line, so one file accumulates a history of
    runs and a wedged or aborted test still leaves everything measured so
    far on disk. By default nothing is written.
- `label` (string) (optional)
  - Free-text label written into the run header in `results_file`, to tell
    runs apart (e.g. `"run3"`).
