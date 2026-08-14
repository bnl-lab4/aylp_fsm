anyloop:kalman_filter
=====================

Types and units: `[T_VECTOR, U_ANY] -> [T_UNCHANGED, U_UNCHANGED]`.

Adaptive linear predictor: replaces each element of the state vector with a
prediction of itself `horizon` samples in the future, so a downstream
controller acts on where the error is *going to be* instead of where it was
one loop-delay ago. Set `horizon` to the total command→error latency in
samples and the predictable part of the disturbance — vibration lines,
resonant humps, anything with autocorrelation longer than the delay — is
corrected as if the delay weren't there; the unpredictable white part passes
through untouched (the least-squares solution learns to ignore it).

Each element gets an order-`order` AR model of the signal, predicting
`s(t+horizon)` from `[s(t) … s(t-order+1)]`, with weights adapted online by
NLMS. True RLS (the exact Kalman filter on the weight vector) was tried first
and diverges on this kind of signal: with exponential forgetting on a
narrowband input, the covariance grows without bound in the unexcited
directions (covariance wind-up). NLMS is the robust stochastic approximation
of the same estimator, costs `O(order)` per element per frame, and reached
the ideal block-least-squares prediction bound on recorded loop data.

Because the weights start at zero, an instant handover would feed the
controller a near-zero error signal and leave the loop effectively open while
the model converges. The handover is therefore ramped:
`output = a·prediction + (1-a)·input`, with `a` rising linearly from 0 to
`gain` over `ramp` seconds once training starts. Training continues forever,
so the model tracks slow changes in the disturbance statistics.

Put it after the sensor (`center_of_mass`) and any logging/UDP sinks (so they
record the *real* error), before `anyloop:pid`.

Parameters
----------

- `type` (string) (required)
  - Must be "vector".
- `units` (string) (optional)
  - Output units. Defaults to null (unchanged from input).
- `order` (int) (optional)
  - AR model order: taps of history per element. Default 60 (25 ms at
    2404 Hz). More taps capture longer correlations but adapt slower.
- `horizon` (int) (optional)
  - Samples ahead to predict. Set to the loop's total command→error latency.
    Default 5 (≈2.1 ms at 2404 Hz, the FSM loop's measured delay).
- `mu` (float) (optional)
  - NLMS step size, `0 < mu < 2`. Smaller = lower steady-state misadjustment
    but slower convergence/tracking. Default 0.03 (offline optimum on
    recorded FSM data was ~0.02; 0.05 tracks faster at a small accuracy
    cost).
- `gain` (float) (optional)
  - Final blend for the x axis (elements ≥ 1), `0..1`. 1 = output the full
    prediction; 0 = passthrough (training continues). Default 1.
- `gainy` (float) (optional)
  - Blend for element 0 (y axis). Defaults to `gain`, same convention as
    `pid`.
- `ramp` (float) (optional)
  - Seconds over which the blend rises from 0 to `gain` after training
    starts. Default 10.
- `clamp` (float) (optional)
  - Output magnitude limit; predictions can overshoot on transients.
    Default 1.0.
- `start_delay` (float) (optional)
  - Seconds to pass the input through untouched *without training* at
    startup. Match the pid/DAC `start_delay` so the open-loop settling
    transient doesn't poison the weights. Default 0.

Results on the FSM tip/tilt loop (2026-07-10)
---------------------------------------------

Offline, the closed-loop error was ~50% linearly predictable at the 2.1 ms
loop delay (x 0.69→0.40 px, y 0.27→0.15 px ideal prediction residual), which
bounds what any causal controller can add at that latency. In the loop
(60 s A/B, steady state):

- x: 0.91 → 0.61 px RMS with the predictor alone; 0.59 px combined with the
  pid's internal-model line oscillators (they coexist — the oscillators
  adapt on whatever the predictor leaves). The 20.24/22.08 Hz lines dropped
  ~5×, 59.63 Hz ~2.4×, and the 0.3–20 Hz band halved. The >200 Hz floor
  rises slightly (0.07→0.09 px): prediction can't help white noise and the
  loop actuates a bit of it.
- y: prediction **hurt** (0.23→0.27 px, the 30.8 Hz line pumped 3×) — y's
  high-gain loop (~62 Hz crossover) reacts badly to the predictor's added
  phase dynamics. Run y as passthrough (`gainy: 0`); its dominant line is
  already handled by the pid oscillator.

The in-loop gain is smaller than the offline bound because feeding
predictions back changes the loop transfer (the predictor and controller
co-adapt); treat the offline bound as a ceiling, not a target.

The benefit scales with how line-dominated the disturbance is, and the lab's
vibration environment is nonstationary (line powers change several-fold
within hours): in a later 120 s verification against a quiet-evening
baseline, closed-loop x jitter matched open-loop (0.63 px both) — the
predictor's line kills were offset by the feedback loop's own waterbed
amplification of the broadband floor. Judge any retuning on ≥60 s runs with
a back-to-back baseline, scored from t > 25 s (settle transient + ramp).
