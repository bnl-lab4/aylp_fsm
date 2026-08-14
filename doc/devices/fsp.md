anyloop:fsp
===========

Types and units: `[T_VECTOR, U_ANY] -> [T_UNCHANGED, U_command]`.

Adaptive **Filtered Smith Predictor** / adaptive LQG for tip-tilt beam
stabilization. Unlike `anyloop:kalman_filter` (a predictor that pre-processes
the error for a downstream `pid`), this is a **complete controller**: error in,
command out. It **replaces both `kalman_filter` and `pid`** — running it with
either double-compensates the loop delay and pumps the vibration lines. Place
it where `pid` used to sit: after `center_of_mass` and any logging sink, before
`clamp`.

The plant seen by the loop is a gain `K`, an optional Bode-fitted
frequency-shaping biquad `H(z)`, and a possibly fractional transport delay
(camera latency + compute + DAC ZOH): `e = phi + K*H(z)*D(z)*u`, where `D`
contains `delay`/`delay_frac` and `phi` is the open-loop disturbance. Three
stages, after
Kulcsár 2006, Petit 2008, Meimon 2010, Correia 2010:

1. **Smith-predictor core.** Reconstruct the disturbance by removing our own
   delayed plant contribution: `phi_meas = e - K*H(z)*D(z)*u`. This takes the
   identified plant out of the disturbance estimate. When `plant_b`/`plant_a`
   are omitted, `H(z)=1` and this is the original scalar-gain model.
2. **Disturbance observer.** The modal path models `phi` as a sum of `n_modes`
   narrowband damped oscillators (`freqs`/`zeta`), each driven by white noise
   of variance `q`. A stationary Kalman gain (precomputed by iterating the
   Riccati recursion at init, and again on each adaptation tick) estimates the
   modal state from `phi_meas`. The `q`/`r` ratio is the LQG tuning knob.
   With `broad_order > 0`, the controller instead uses an identified
   full-band Wiener/Kalman FIR realization, covering both the broadband
   continuum and vibration peaks.
3. **delay-step prediction + minimum-variance control.** Roll the modal state
   forward `delay` steps and cancel it. With a shaped plant, the requested
   actuator-space correction is passed through the matched stable inverse:
   `u = H^-1(z)*(-phi_hat/K)`. Requiring both `H` and `H^-1` to be stable keeps
   this causal and prevents the fit from hiding an unstable plant inversion.
4. **Optional drift/transient split.** With `drift_tau > 0`, a slow EWMA state
   carries DC/random-walk pointing drift and the modal/FIR predictor sees only
   the remaining vibration. With `transient_sigma > 0`, an unexpectedly large
   prediction innovation cross-fades from predictive
   cancellation to a bounded proportional recentering servo. After the
   innovation remains quiet for `transient_hold`, the predictor is blended
   back over `transient_ramp`. With `transient_modal_q_scale > 0`, recovery
   instead uses a high-gain copy of the configured damped-mode observer to
   estimate displacement and velocity and predict the released ring-down,
   retaining a small proportional correction for model error. This gives bumps
   a direct event path without permanently placing a second controller around
   the predictive loop.

This is the one in-loop mechanism that beats the delay-limited crossover on the
**predictable** (modal) part of the disturbance. It does nothing for the
broadband floor or above-bandwidth content — a fundamental limit.

Adaptation (`adapt_period > 0`): every `adapt_period` s, refresh each mode's
`q` from its recent estimated energy, `r` from the innovation floor, nudge each
center frequency toward the locally demodulated line (capped at `adapt_df_max`
Hz/update, sized to the ~0.5 Hz/run wander), rebuild coefficients, and recompute
both the normal and event-only modal gains. The two Riccati solves are
amortized serially across frames; modal event entry is briefly inhibited until
the event gain matches the updated model. Set `adapt_period <= 0` for a fixed
FSP while first validating.

Parameters
----------
- `type`: must be `"vector"`.
- `units`: output units (e.g. `"minmax"`).
- `delay`: loop transport delay in samples. Fit it from the full-band Bode
  phase slope plus the one hidden bode frame (2026-07-16 @ 3788 Hz: x 5.62,
  y 6.37 frames), NOT from the step-departure latency, which understates the
  effective delay. Beware the units when the frame rate changes: latency in ms
  can fall while the delay in FRAMES rises. Each axis object may set its own
  `delay`/`delay_frac`, overriding this global value (the two plants differ).
- `delay_frac`: fractional remainder of the transport delay. The Smith plant
  uses a first-order Thiran all-pass, preserving command magnitude while
  matching the fractional group delay; the full-band observer blends its
  adjacent horizon predictions. Per-axis override allowed as with `delay`.
- `delay`: `"auto"` — measure the delay from the running loop instead of
  fitting it, setting both `delay` and `delay_frac`. Works globally or on one
  axis. Two of the three terms are measured: the **compute** time, published by
  the source device (only the device that blocks for the frame can tell work
  from waiting — for anything downstream the gap since its own return is
  identically the loop period minus its own duration), and the **loop period**,
  measured here. The third, `delay_auto_bias`, is the part no device in the
  loop can see: sensor integration and readout, the DAC's zero-order hold, and
  actuator lag. Measurement runs at the start of the hold and is installed as
  soon as its window closes — still inside the hold, so the broadband filter
  spends the rest of the hold training at the horizon it will run at;
  `start_delay` is raised at init if it is too short. The source must publish `libaylp/timing.h`
  timing — `asi_source` does — or init falls back, loudly, to the configured
  value. **An auto delay is an estimate, not an identification**: the compute
  term includes sinks running after the DAC, and the bias term is assumed.
  Verify with a Bode sweep before quoting an attenuation number.
- `delay_auto_bias`: frames of sensor + ZOH + actuator latency added to the
  measured compute term. Default 1.5 (one frame of sensor, half a frame of
  hold). `"auto"` recovers it from `delay_ident_ms` instead —
  `(identified delay − measured compute)`, converted at the measured rate — so
  it stays right when the frame rate moves. A frame count written by hand is
  only valid at the rate it was written for.
- `delay_ident_ms`: the last **identified** end-to-end transport delay, in
  milliseconds, from a Bode sweep or a PRBS run. Used by
  `delay_auto_bias: "auto"`. Milliseconds because that is the unit that does
  not move with the frame rate. `"auto"` means no identification is being
  asserted: the sensor path is modelled from the source's published exposure as
  `exposure/2 + 2.5 frames`. The 2.5 is measured, not assumed — removing the
  exposure half from four Bode identifications on this bench leaves 2.53–2.67
  frames across an 8× range of frame rate, i.e. ~2 frames of camera pipelining
  plus the half-frame ZOH, and the model reproduces those to within 4%. It is
  still a model: it cannot see actuator lag and it assumes this camera's
  pipelining. The calibration suite replaces it with a measured number on the
  first PRBS delay run.
- `delay_auto_max`: ceiling the estimate is clamped to, and the delay every
  ring is sized for while the real one is unknown. Default 16.
- `delay_auto_settle`, `delay_auto_window`: seconds of hold discarded before
  measuring, and seconds measured. Defaults 0.5 and 1.0.
- `fs`: loop rate (Hz); must match, so AR coefficients land on the right digital
  frequencies. `"auto"` measures it from the same startup pass as the delay and
  rebuilds everything derived from it — mode coefficients and their Riccati
  gains, the burst guard, the command filter, and the EWMA weights converted
  from time constants. Use it with a self-sizing camera ROI, where the frame
  rate is not knowable when the config is written.
- `gap_trip`: seconds between `proc()` calls past which the source is taken to
  have dropped frames. `"auto"` sets it to 1.5 measured frame periods, so a
  normal frame passes and a single missed one trips — a fixed number of seconds
  only means one thing at one frame rate.
- `clamp`: command magnitude limit — the symmetric shorthand for
  `[-clamp, +clamp]`. Default 1.0.
- `clamp_min`, `clamp_max`: the actual output bounds, for an **asymmetric**
  limit. Either one overrides the corresponding side of `clamp`, in any config
  order. Use these when the actuator's reachable range is not centred on the
  command origin: a mirror spanning −10…+10 V driven from a +2.5 V bias at
  1 V/command-unit is `clamp_min: -12.5, clamp_max: 7.5`. A symmetric limit
  there would either give up half the negative range or let the loop request a
  voltage the hardware cannot produce.
  - The window must contain 0 (init fails otherwise): the `start_delay` hold,
    the `ramp` blend, the burst guard and the non-finite fallback all drive the
    command to zero, so a window excluding it would make those unrepresentable.
    Zero exactly on a boundary is allowed.
  - **Put the limit here, not in a downstream `anyloop:clamp` stage.** The
    clamped command is also what enters the fractional-delay filter and command
    ring, so bounding it here keeps this device's plant model identical to what
    the actuator received. A clamp *after* `fsp` would cut the command without
    `fsp` knowing, and that divergence is the command-echo mismatch that drives
    regeneration. A downstream `clamp` stage matched exactly to these bounds is
    fine as a backstop — just never tighter.
- `start_delay`, `ramp`: hold the command at 0 for `start_delay` s (Kalman
  converges on the clean open-loop error), then blend to full authority over
  `ramp` s.
- `adapt_period`, `adapt_df_max`, `adapt_tau`: adaptation cadence, per-update
  frequency-correction cap, and EWMA time constant.
- `cmd_fc`: optional 2nd-order low-pass on the command with per-mode
  roll-forward phase compensation (0 = off, the default). Leave it off: the
  `fsp_sim.py` study showed the waterbed does **not** come from HF command
  noise, so the filter costs prediction horizon without buying anything. The
  real waterbed fix is the q/r scaling below.
- `broad_order`, `broad_mu`: enable the full compound-disturbance observer.
  The observer learns the delay-step conditional mean directly from the
  Smith-reconstructed disturbance using an order-`broad_order` FIR
  realization of the scalar Wiener/Kalman predictor, updated by NLMS with
  step `broad_mu`. This supplies the command instead of the modal estimate
  (the two must not be summed). Set `broad_order` to 0 for the original
  modal-only controller. Attenuation12 uses order 128 and `broad_mu=0.005`;
  attenuation11 showed that 512 taps at 0.03 preserved more high-frequency
  coefficient noise than the additional prediction accuracy justified.
- `wiener_file`: optional offline initialization for the full-band predictor.
  The UTF-8 text file must contain exactly `broad_order` non-comment rows,
  indexed from zero, with columns
  `index y_w y_w_next x_w x_w_next`. Blank lines and lines beginning with `#`
  are ignored. The offline solve must use the same `broad_order`, per-axis
  horizon, `drift_tau`, and `broad_lp` preprocessing as the runtime. After
  loading, `broad_freeze_closed` and `broad_mu` determine whether the solution
  stays fixed or continues adapting.
- `wiener_out`: optional path to dump the **learned** weights at exit, in
  exactly the format `wiener_file` reads, so a converged run can be analysed
  offline and then replayed as the next run's initialization. The dump happens
  in `fini`, which runs on both the `AYLP_DONE` and the `SIGINT` path, so an
  aborted run still saves whatever the observer had learned. The path is
  probed for writability at init, so a bad path fails immediately rather than
  after the run. Commented header lines record `broad_order`, `broad_mu`,
  `broad_lp`, `fs`, the per-axis horizon and `K`, the closed-frame count, and
  the per-axis tap norms `||w||` — the norms are the cheap convergence check,
  since a converged observer settles and a diverging one grows. Reload only
  into a run with the same `broad_order`, horizon, `broad_lp` and `drift_tau`:
  taps are meaningful only against the model they were learned under, and the
  loader's contiguous-index check cannot catch a horizon mismatch.
- `wiener_trace`, `wiener_trace_period`: optional convergence trace, one text
  line every `wiener_trace_period` seconds (default 10). Columns: elapsed
  seconds, closed-frame count, per-axis `||w||`, `||w_next||` and `||dw||`
  (the tap change since the previous sample), per-axis `drift_hat`, and the
  cumulative per-axis guard and transient counts. `wiener_out` says where the
  observer *ended up*; this says whether it had **stopped moving**, which is
  the question a `settle_time` exists to answer — `||dw||` must have flattened
  before the scored window opens. Deliberately a short line rather than a tap
  snapshot, so the periodic write cannot disturb the RT loop; it is flushed
  every sample so the trace survives an abort. With `pass_open: false` the
  observer is fed zeros during the open phase and the NLMS step is
  `mu*pe/energy` with `energy` floored at `1e-12`, so `||dw||` is exactly 0
  there — a nonzero value before the loop closes means open-loop data is
  leaking into the observer.
- `broad_lp`: observer band-limit (odd boxcar taps on the NLMS input; 0 = off).
  The full-band observer is otherwise sensitive out to Nyquist, and any K/delay
  model error leaks command echo into its input as a *predictable* HF signal
  that closed-loop NLMS learns and chases into a self-sustained ring at the
  regeneration frequency (2026-07-22: a ~380 Hz x-axis limit cycle, 4.8 px rms,
  ignited ~1 min after close with same-day bode-matched K/delay — the residual
  mismatch is amplitude-dependent, so no static tune removes it). A `broad_lp`
  boxcar has its first null at `fs/broad_lp` and attenuates the whole
  regeneration band; being linear-phase with exactly integer group delay
  `(broad_lp-1)/2`, that delay is folded into the broad prediction horizon so
  in-band cancellation timing is unchanged (11 taps at 3788 Hz: null 344 Hz,
  <1 % droop at 30 Hz). This removes at the root what the burst guard only
  reacts to; verified 180 s with zero ring vs. a latched limit cycle without it.
- `dark_predict`: keep commanding when an upstream centroid device marks a
  frame `AYLP_FRAME_REJECTED` (default true). A bank of direct multi-horizon
  NLMS predictors reads only the baseline-masked history ending at the last
  real sample. On dark frame `k`, horizon `k + delay + broad_gd` supplies the
  command directly; the adjacent horizon supplies fractional-delay blending.
  Bank output never enters the ordinary observer history, preventing a
  prediction from feeding itself or contaminating later live frames. Disable
  this to retain the previous mask-and-hold behavior.
- `dark_predict_max`: consecutive dark frames to predict before holding again.
  Zero selects `2 * (delay + broad_gd)` per axis. The bank internally includes
  the additional command horizon and adjacent fractional-delay filter.
- `dark_bank_train`: bank horizons updated per clean live frame (default 4),
  round-robin. This bounds the added per-frame work.
- `dark_train_min_real`: minimum real fraction of the ordinary FIR regressor
  required for an update (default 0.5). Baseline-masked taps remain valid
  prediction inputs but are excluded from both the update error and the
  coefficient update.
- `dark_fab_frames`: consecutive dark frames treated as a full outage, after
  which ordinary FIR training pauses until fabricated history has drained.
  Zero selects half the history length; leave it automatic for short chops.
- `drift_tau`: slow-drift EWMA time constant in seconds. `<= 0` preserves the
  original compound predictor. When enabled, the drift estimate is cancelled
  separately as constant over the command horizon and is subtracted before
  the boxcar/modal/FIR vibration observers. Start around 1--5 s: shorter values
  transfer more low-frequency motion out of the predictor, while longer values
  leave more of it in the learned vibration model.
- `drift_order`: slow internal-model order (default 1). `1` retains the EWMA
  position estimate. `2` uses an alpha-beta position/rate observer and
  propagates the rate through the command horizon. This embeds DC integral and
  ramp-drift rejection in the normal Smith predictor; it does not enable or
  depend on the event-only transient integral controller.
- `transient_mode`: `manual` (default) honors the configured values; `off`
  disables event recovery regardless of retained tuning; `auto` enables the
  conservative default detector (`transient_sigma = 6` when unset) while
  retaining the other explicit values/defaults. Use `off` in an unvalidated
  control profile rather than relying on several zero-valued fields.
- `transient_sigma`, `transient_floor`: enable large-event recovery and set its
  innovation threshold to
  `max(transient_sigma * quiet_innovation_rms, transient_floor)`, in normalized
  error units. `transient_sigma <= 0` disables the path. The quiet variance is
  learned with time constant `transient_tau` whenever no event is active.
- `transient_settle_error`: absolute raw-error band used to score output
  settling in `transient_log` (default 0.03 normalized). It does not affect the
  detector or recovery state machine. The event's output settling time is the
  last crossing of this band, while its innovation-quiet time separately says
  how long the underlying reconstructed disturbance remained unexpected.
- `transient_controller`: event recovery law. Supported values are:
  - `proportional`: `u = -(drift_hat + Kp*error)/K`.
  - `integral`: `u = -drift_hat/K + I`, a pure event-local integrator.
  - `pi`: proportional plus the event-local integrator.
  - `modal`: fast event-only modal Kalman prediction plus proportional
    correction (the previously implemented modal recovery).
  - `hybrid`: modal prediction plus proportional and integral correction.
  The default is `proportional`. For backward compatibility, an unnamed
  controller with `transient_modal_q_scale > 0` selects `modal`.
- `transient_kp`: proportional gain used by every mode except `integral`; it
  must satisfy `0 < transient_kp < 1`. Preserving `drift_hat` in every mode
  prevents an event from releasing the existing slow pointing correction.
- `transient_ki`: integral gain in 1/s (default 5). The integral contribution
  is stored directly in command units and obeys
  `dI/dt = -transient_ki*error/K - transient_i_leak*I`.
- `transient_i_leak`: integral leakage rate in 1/s (default 0). Zero gives a
  mathematically pure integrator; a positive value provides a leaky integrator.
- `transient_i_limit`: absolute bound on the integral command contribution
  (default 1 command unit). Integral modes require a positive value. The state
  resets at event entry/completion, timeout, and beam loss. It stops integrating
  during the return cross-fade, and conditional-integration anti-windup rolls
  back any update that would push an already clamped actuator farther into
  saturation.
- `transient_modal_q_scale`: enable an event-only physical ring-down observer
  by multiplying the configured modal process-noise variances by this factor
  when its Kalman gain is solved (`0` disables, preserving the proportional
  fallback). The normal modal gain and the full-band predictor are unchanged.
  During an event, the higher-gain observer estimates each damped AR(2) mode's
  displacement/velocity state, propagates it through `delay + delay_frac`, and
  uses
  `u = -(drift_hat + modal_prediction + transient_kp * error) / K`.
  Thus a repeatable push-excited sinusoid is anticipated rather than integrated,
  while the small P term covers motion outside the physical model. Start around
  10--100 and compare event settling; excessive values make measurement noise
  look like modal state. There is no separate transient plant calibration:
  Smith reconstruction, modal prediction, burst-guard placement, and event
  command scaling all use the same per-axis `K`, `delay`, and `delay_frac` as
  normal operation. After a new Bode fit, update those fields once. Online
  modal adaptation refreshes both the normal and event-only Kalman gains.
- `transient_modal_ab`: when the selected controller is `modal` or `hybrid`,
  alternate recovery per axis between proportional (odd-numbered events) and
  the selected controller (even-numbered events). The selected recovery is
  written to `transient_log`. This preserves the same-run manual-push A/B.
- `transient_tau`, `transient_hold`, `transient_ramp`: quiet innovation-variance
  time constant, required quiet time after the most recent threshold crossing,
  and bumpless cross-fade time back to predictive cancellation. The broadband
  NLMS weights freeze from event entry through the return cross-fade and resume
  on the first fully normal frame, so a push cannot rewrite the stationary
  predictor or teach the release detector its own event. Slow modal covariance
  and frequency statistics also remain frozen until recovery completes.
- `transient_shadow_mu`, `transient_shadow_tau`,
  `transient_shadow_min_duration`, `transient_shadow_hold`,
  `transient_shadow_ratio`, `transient_shadow_norm_ratio`: optional escape from
  a genuinely changed stationary environment while preserving the event-time
  freeze. At event entry FSP clones both broadband FIR horizons into a private
  shadow. The production weights remain fixed; the shadow learns with
  `transient_shadow_mu`, and its pre-update (prequential) squared error is
  averaged over `transient_shadow_tau`. It may replace the production model
  only after the event is at least `transient_shadow_min_duration` old, its RMS
  error remains below `transient_shadow_ratio` times the frozen model for
  `transient_shadow_hold`, and its coefficient change is no more than
  `transient_shadow_norm_ratio` times the production norm. Any event command
  saturation or burst-guard suppression disqualifies the candidate. After
  promotion, conservative shadow updates are copied into the active model until
  the ordinary innovation-quiet criterion releases recovery. Thus an impulse
  ring-down cannot contaminate normal NLMS, while a persistent predictable
  regime cannot deadlock forever against a frozen obsolete model. Set
  `transient_shadow_mu <= 0` to disable (the default). Each event CSV row records
  `shadow_promoted` and its held-out `shadow_error_ratio`.
- `transient_clamp` (or `transient_clamp_min`/`transient_clamp_max`): optional
  command envelope used only while an event is active. It must contain the
  normal `clamp` window. Outside an event, the normal clamp still applies. Any
  downstream backstop must match this wider envelope or it will silently break
  FSP's command echo during recovery.
- `transient_trip_command`: pre-clamp command trip used during an event. It
  defaults to the normal `trip_command`. Keep it no higher than the actuator's
  physical limit.
- `transient_max_duration`: hard event-authority timeout in seconds. If an
  event remains active longer than this, FSP returns that axis to its normal
  command envelope and inhibits another event until the innovation becomes
  quiet. It does not latch zero. Zero disables the wall-clock criterion, so an
  active event then ends only after innovation remains below the learned model
  threshold for `transient_hold` plus the `transient_ramp` handoff, or when the
  upstream sensor reports beam loss. Temporary authority reductions do not
  count as quiet innovation.
- `transient_arm_delay`: seconds after the startup hold ends (i.e. after
  `start_delay`) before event detection is enabled (default 0). During this
  interval the closed-loop predictor continues learning and the residual-
  variance EWMA is calibrated continuously; only triggering is inhibited.
  Exists because a cold (zero-weight) full-band predictor's own convergence
  transient looks exactly like the event this path is meant to catch: `pe` is
  large right at close simply because the filter hasn't learned the
  disturbance yet, and depending on what fed the predictor during the hold,
  the detector needs a post-close warm-up so the threshold is based on the
  converged closed-loop residual rather than an artificially quiet hold (e.g.
  a zeroed input during `attenuation_test`'s open phase). A run with
  `broad_order=512`,
  the `broad_mu_init`/`broad_mu_tau` schedule and this rig's disturbance
  spectrum was verified offline (replaying a recorded open-phase disturbance
  through the same filter) to reach steady prediction-error character by
  `t=300 s` after real data starts reaching the predictor; `transient_arm_delay`
  should be set with margin above that for production scoring. Because NLMS is
  intentionally frozen during an event, this warm-up must be long enough that
  a cold predictor cannot false-trigger before it has converged. Note this
  delay is measured from `start_delay` (the
  hold), not from process start -- if nothing feeds the predictor real data
  until close (e.g. `attenuation_test` with `pass_open: false`), those are the
  same instant anyway.
- `transient_log`: optional CSV path. One row is flushed after every completed
  axis-event with its recovery type, start time, output-error settling time,
  innovation-quiet time, full recovery/handoff time, peak absolute error, event
  RMS error, peak absolute clamped command, peak absolute integral contribution,
  and sample count. This is intended for repeatable manual push/ring-down trials;
  unlike an attenuation PSD, it scores the transient directly.

The full-rate error and applied-command records from
`contrib/legacy/legacy_steering_profiles/push_event_par_fsp.json` can also validate whether the modal frequencies
survive a move to a new environment:

```sh
python3 contrib/data-analysis-tools/analyze_push_events.py push_modal_events.csv \
  --frequency-test --error-aylp push_modal_error.aylp \
  --command-aylp push_modal_command.aylp \
  --config contrib/legacy/legacy_steering_profiles/push_event_par_fsp.json --strict
```

For every detected push, the analyzer uses the same `K`, integer/fractional
delay, and optional plant biquad as FSP to reconstruct
`disturbance = measured_error - K H(z) D(z) applied_command`. It estimates each
configured line in a fixed post-trigger window, then reports the median and MAD
across pushes. Frequency is the invariant here: push amplitude, direction and
duration primarily change modal amplitude and are handled by the per-event
line-prominence gate. `UNTESTED` means too few pushes visibly excited that mode;
it is not a pass. A repeatable strong residual line is a failure and a candidate
mode to add or substitute. Use `--help` to tune window, tolerance, excitation,
and residual thresholds. Error and command streams must start on the same loop
frame, as they do in the supplied push pipeline.

For primary modal identification, use the open phase of a completed attenuation
run instead. The fitter reads the pre-FSP AYLP record and obtains the open-window
boundaries and sample rate from the run configuration:

```sh
python3 contrib/data-analysis-tools/fit_fsp_modes.py atten_par_err6.aylp \
  --config contrib/steering/configurations/support/attenuation_par_fsp.json \
  --output fitted_fsp_modes.json
```

It prints per-axis `freqs`, `zeta`, and `q` with `r` normalized to 1 (only the
`q/r` ratio affects the Kalman gain). Candidate lines are selected by open-loop
power and prominence; each is fit to a Lorentzian PSD, and its fitted modal
variance is converted to the AR(2) drive variance expected by FSP using the same
stationary `Gv` formula as the controller. A `zeta resolution-limited` warning
means the open record/FFT cannot resolve that line's damping: retain a
conservative prior or obtain damping from a longer record/released ring-down
rather than treating the reported lower bound as exact. Review and copy the
fragment into attenuation, steering and push configs; the tool deliberately
does not rewrite production configurations automatically.

A conservative first hardware trial (all values remain opt-in unless included):

```
"drift_tau":       2.0,
"transient_sigma": 6.0,
"transient_floor": 0.03,
"transient_controller": "proportional",
"transient_kp":    0.25,
"transient_tau":   5.0,
"transient_hold":  0.10,
"transient_ramp":  0.25
```

Set `transient_floor` from a stable run rather than blindly keeping the example:
it should sit above ordinary single-frame centroid scatter but below the event
size that needs direct recovery. Count/logged activations should be rare.

- `broad_freeze_closed`: freeze full-band identification when the startup hold
  ends (default true, retained as a conservative legacy option). Set false for
  a continuously adaptive closed-loop observer; the supplied parport steering,
  attenuation and push configurations all do this and do not depend on fixed
  open-loop weights. Continuous closed-loop NLMS can identify leaked command as
  disturbance when `K` or the
  delay model is imperfect and thereby create positive feedback.
- `trip_error`, `trip_command`, `trip_frames`: non-latching diagnostics. During
  the startup hold FSP learns each axis's ordinary open-loop operating point.
  After closing, it warns if error magnitude exceeds the learned open-loop
  magnitude by `trip_error`, or if requested (pre-clamp) command exceeds
  `trip_command`, for `trip_frames` consecutive samples. The normal or event
  clamp still contains the command. This magnitude-envelope test permits
  successful motion from a static offset toward zero without mistaking the
  offset for runaway. Only the upstream `AYLP_BEAM_LOST` validity flag holds
  FSP at zero; tracked center-of-mass asserts it after `reacquire_after`
  consecutive frames without signal and clears it when the beam is found.
- `beam_recover_ramp`: seconds over which normal authority returns after
  `AYLP_BEAM_LOST` clears (default 0.5). On loss, dynamic observer, command-ring,
  and filter state is cleared while learned FIR weights are retained. Stale held
  centroids are neither predicted nor used for adaptation. Reacquisition then
  rebuilds live history and ramps from zero without a command step.
- `guard_ratio`, `guard_floor`, `guard_hold`, `guard_ramp`, `guard_tick`:
  non-latching **burst guard** (default ON: 4 / 0.008 / 0.25 s / 1 s / 10 s;
  `guard_ratio <= 0` disables). Any mismatch between the configured plant
  model (`K`, `delay`) and the true plant leaks command back into the Smith
  reconstruction, closing a parasitic loop whose phase crosses −180° at
  `fs / (2·(delay + delay_frac))` (~316 Hz y / ~337 Hz x at 3788 Hz). When
  the margin there goes negative — `K` drifts with the coarse bias and
  alignment — the loop emits intermittent ring bursts at that frequency
  (the 2026-07-16 "310 Hz spiral"; the 2026-07-17 run-15 rerun, where the
  bursts averaged into a broad 250–450 Hz PSD hump and tripled the closed
  RMS while the code was identical to the good morning run). Per axis, the
  guard band-passes the raw error at the regeneration frequency (Q 1.5),
  tracks a ~10 ms envelope against a ~10 s quiet baseline, and triggers when
  the envelope exceeds `guard_ratio` × max(baseline, `guard_floor`). On a
  trigger the axis's authority is cut to 0, NLMS training and the adaptation
  statistics freeze (so the ringing is never learned), the cut holds for
  `guard_hold` s (extended while the ring persists), then authority ramps
  back over `guard_ramp` s; a burst returning mid-ramp counts as a new
  event. Every `guard_tick` s of closed-loop time a ticker line reports the
  activation count and the percentage of frames spent at reduced authority,
  and a final total is logged at shutdown. The guard keeps a marginal run
  alive; **recurring activations mean `K`/`delay` no longer match the plant
  at the current operating point — re-run the Bode fits before trusting the
  attenuation numbers.**
- `gap_trip`: **stall-gap guard** (seconds; default 0.05, `<= 0` disables).
  If the wall-clock gap between consecutive frames exceeds this — e.g. the
  ASI camera stream stalled and asi_source restarted capture, a ~0.3–0.6 s
  hole during which the DAC held the last command while the error kept
  moving — every predictor history now spans the hole and re-entering at
  full authority onto stale state is exactly the transient that seeds a
  regeneration burst. An authority *cut* is wrong too: the command carries
  the loop's DC correction, and zeroing it snaps the FSM back to bias,
  re-exposing the accumulated drift as a multi-px excursion that the ramp
  then walks back slowly (observed 2026-07-17: ~6 px for ~1.3 s). Instead,
  each axis **keeps commanding the value the DAC held through the gap**
  (zero bump at resume) for `guard_hold` s while the predictor state is
  rebuilt — modal state zeroed (mode phases rotated unpredictably), NLMS
  input history reset so it never predicts or trains across the hole, Smith
  command ring rewritten with the held command — then **blends from the
  held command to the live controller command** over `guard_ramp` s.
  Learning stays frozen until the blend completes; an active upstream beam-loss
  hold overrides the blend. Works even with the burst detector disabled. Gap
  events are counted and logged separately from burst events: **bursts mean
  re-measure `K`; gaps mean the camera stalled** (check the asi_source
  recovery lines).
- `y`, `x`: per-axis objects (element 0 = y, element 1 = x), each with:
  - `K`: **signed** command→error gain. Wrong sign = positive feedback =
    runaway; verify with a push test.
  - `plant_b`, `plant_a`: optional three-coefficient numerator and denominator
    of the per-axis Bode residual `H(z)=(b0+b1*z^-1+b2*z^-2)/
    (a0+a1*z^-1+a2*z^-2)`. Both arrays must be supplied together. FSP
    normalizes `a0`, verifies that all poles and zeros lie inside the unit
    circle, uses `H` in the Smith reconstruction, and uses the matched `H^-1`
    on the cancellation command. Omit both for the original `H=1` behavior.
    Run 16 coefficients are fits to `fsm_bodex_716`/`fsm_bodey_716` after
    removing the run-15 `K`, transport delay, and the Bode helper's hidden
    frame; they cut 5–300 Hz complex plant-model RMS error from 8.80% to
    4.39% on x and 8.31% to 3.46% on y.
  - `r`: measurement-noise variance.
  - `freqs`: mode center frequencies (Hz), up to `AYLP_FSP_MAX_MODES` (8).
  - `zeta`: per-mode damping (~0.002 for a sharp line, ~0.3 for a broadband
    hump). Defaults to 0.002.
  - `q`: per-mode process-noise (**drive**) variance (array), or a scalar for
    all. Defaults to 1e-5. NOT the mode's visible energy: drive variance =
    state energy / Gv, where Gv ≈ 1/(4π·zeta·f·Ts) is huge for sharp lines.

Tuning: q/r is the waterbed dial (`fsp_sim.py` on the measured attenuation10 PSD)
---------------------------------------------------------------------------------
The rejection transfer is `R(ω) = 1 − e^{−jωd}·T(ω)`, with `T` the
measurement→command transfer. A conservation law identical to the feedback
Bode integral applies (`∫ ln|R| dω = 0` for stable causal `T`): in-band
rejection must be paid for out of band. The q/r ratio positions the design on
that tradeoff curve:

- **q/r ~ 1** (mis-scaled): the estimator tracks *everything*, `|T| ≈ 1` far
  out of band, and the delay phase mismatch alone makes
  `|R| = 2|sin(ωd/2)| ≈ 2–4` — a ~3× waterbed at 64–400 Hz. The Kalman gain
  also inflates (close-spaced mode pairs fight over the shared measurement —
  same pathology as the pid notes' "keep oscillators >2 Hz apart").
- **q/r ~ 1e-5** (physical): waterbed ~×1.5 — the same as the classical
  pid+lead loop — with the 10–30 Hz band rejected ~2× harder (sim: x 0.060 px
  vs 0.137 measured). Total RMS lands **at** the classical loop's
  (x 0.241 vs 0.244, y 0.161 vs 0.151 px): the conservation law means the FSP
  *redistributes* error (deeper lines, same waterbed) rather than shrinking
  the total on this disturbance mix.
- **q/r ~ 3e-6**: waterbed ~×1.35 at softer in-band rejection.

Corollary: don't judge the FSP by total RMS alone — its value over the
classical loop is *where* the residual sits (line depth vs broadband), plus
single-knob retuning. For a smaller total, the delay itself must shrink (QPD
front-end); no in-loop tuning escapes the conservation law.

The adaptation refreshes q from the measured mode state energy **divided by
Gv** (the stationary variance gain) — feeding state energy in directly
re-creates the q/r ~ 1 regime and its ×3 waterbed.

Validate with `contrib/legacy/outdated_scripts/attenuation_fsp.json`
(open-vs-closed A/B) and compare
the closed RMS-about-mean and the 120–400 Hz band ratio against
`attenuation10_lead_fc64` (x 0.244 px / ×1.57, y 0.151 px / ×1.32).
