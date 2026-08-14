anyloop:pid
===========

Types and units: `[T_VECTOR|T_MATRIX, U_ANY] -> [T_UNCHANGED, U_UNCHANGED]`.

This device applies [PID](https://en.wikipedia.org/wiki/PID_controller) error
correction, taking the pipeline state as error input and replacing it with the
calculated correction.

Let `x_0` be the current input, `x_1` be the input from last iteration of the
loop, `dt` be the time in seconds that has elapsed since the previous iteration,
`x_acc` be the accumulated total errors since the loop was started, and `y` be
the output. Then this PID device more or less implements:

```
y = - p*x_0 - i*x_acc - d*(x_0-x_1)/dt
```

More precisely, it implements the above equation but clamping the total
correction to the `clamp` parameter in magnitude, and bounding the accumulator
so that the integral term `i*x_acc` can reach that same limit without ever
winding past it. See the anti-windup note under `clamp` below.

Parameters
----------

- `type` (string) (required)
  - "vector" if we are expecting `T_VECTOR` input, and "matrix" if we are
    expecting `T_MATRIX` input.
- `units` (string) (optional)
  - Output units, e.g. "V", "rad", etc. Defaults to null (unchanged from input).
- `p` (float) (optional)
  - Coefficent for proportional correction (see equation above). Defaults
    to 1.0.
- `i` (float) (optional)
  - Coefficent for integral correction (see equation above). Default 0.0.
- `d` (float) (optional)
  - Coefficent for derivative correction (see equation above). Default 0.0.
- `clamp` (float) (optional)
  - What to clamp the correction to in magnitude. Defaults to 1.0.
  - Anti-windup: the accumulator reaches the output only as `i*x_acc`, so it is
    bounded to `clamp/|i|` rather than to `clamp`. Bounding `x_acc` itself at
    `clamp` (as this device used to) lets the integral term range over
    `±i*clamp` — `i` times further than the output can ever express — so after a
    saturating transient the accumulator has to unwind all that excess before the
    correction even leaves the rail. At `i=300` that measured ~870 ms of dead
    time; with the bound it is one sample. On top of that, when the correction is
    already saturated, a sample that would push it further into the rail is not
    integrated at all (the `g` leak still applies).
- `start_delay` (float) (optional)
  - Seconds to hold the correction at zero, with the accumulator kept empty,
    before the loop is allowed to act. Measured from the first iteration.
    Defaults to 0.
  - Use this when something slow has to settle first — a coarse stage walking the
    beam toward centre, or `center_of_mass`'s wide acquisition phase. Holding the
    actuator alone is not enough: the integrator would keep winding on the large
    startup error behind the parked output and arrive at its rail the instant the
    hold expired. If you also park DAC channels with `piplate_bridge`'s
    `start_delay`, make this at least as long.

Coarse channels (dual-stage mode)
---------------------------------

With `"coarse": true` (vector type only), the device drives a two-stage
actuator from a single 2-element `[y, x]` error vector (e.g. from
`center_of_mass`). The output becomes 4 elements:

```
[0] x_fine   [1] y_fine   [2] x_coarse   [3] y_coarse
```

The fine channels are the ordinary PID described above (base params are the x
axis, `py`/`iy`/… the y axis, as usual). Each coarse channel is an independent
PID with its own gains, gated by a hysteresis comparator on the **voltage its
fine channel is commanding** — i.e. it engages when the fine actuator is
running out of travel. The device reconstructs that voltage as
`fine_offset + cmd * fine_scale`, so those params must mirror the DAC stage's
(`piplate_bridge`'s) `scale`/`offset` for the fine channels:

- **Activate** when the fine channel's voltage leaves
  `[coarse_on_low_v, coarse_on_high_v]` (default `[0.5, 3.5]` V).
- **Deactivate** when the fine channel returns within
  `[coarse_off_low_v, coarse_off_high_v]` (default `[1, 3]` V). On
  deactivation the channel's output **freezes** at its current value rather
  than returning to zero: the coarse stage has taken over the offset that was
  railing the fine channel and must hold that position — snapping back to bias
  would immediately rail the fine channel again. The gap between the on and
  off windows is what keeps the gate from chattering; the off window must sit
  inside the on window.

Downstream, pick the outputs by element with `piplate_bridge`'s `index` param,
e.g. `"index": [1, 0, 3, 2]` for a y/x fine, Y/X coarse DAC channel order.

Coarse-mode parameters:

- `coarse` (boolean) (optional)
  - Enable the 4-channel dual-stage output. Defaults to false.
- `coarse_p`, `coarse_i`, `coarse_d` (float) (optional)
  - PID gains for the x coarse channel, same convention as `p`/`i`/`d`. All
    default to 0.0, so an enabled but untuned coarse channel commands nothing.
- `coarse_g` (float) (optional)
  - Integrator leak for the x coarse channel, like `g`. Defaults to 1.0.
- `coarse_py`, `coarse_iy`, `coarse_dy`, `coarse_gy` (float) (optional)
  - The y coarse channel's gains; default to the corresponding x value.
- `fine_scale`, `fine_offset` (float) (required in coarse mode)
  - Volts per pipeline unit and volts at zero command for the **x fine**
    channel, exactly as configured on the DAC stage driving it. Used only to
    reconstruct the fine channel's commanded voltage for the gate. A zero
    `fine_scale` is a config error.
- `fine_scaley`, `fine_offsety` (float) (required in coarse mode)
  - Same, for the y fine channel.
- `coarse_on_low_v`, `coarse_on_high_v` (float) (optional)
  - Activation window edges, in volts: the coarse channel turns on when its
    fine channel's voltage drops to `coarse_on_low_v` or below, or rises to
    `coarse_on_high_v` or above. Defaults 0.5 and 3.5.
- `coarse_off_low_v`, `coarse_off_high_v` (float) (optional)
  - Deactivation window, in volts: the coarse channel turns off (freezing its
    output) once the fine channel's voltage is back within
    `[coarse_off_low_v, coarse_off_high_v]`. Defaults 1 and 3. Values outside
    the on window are clamped to it (with a warning), which degenerates to no
    hysteresis on that edge.

The coarse channels share the fine channels' `clamp` and `start_delay` (during
the startup hold they are parked at zero and cannot trigger), and use the same
anti-windup scheme. They do not apply the fine channels' deadband, derivative
filter, or lead compensator.

Narrowband line rejection (internal model)
------------------------------------------

With `lines`/`linesy` set (vector type only, 2-element `[y, x]` error), the
device adds one adaptive quadrature oscillator per listed disturbance line to
the corresponding axis's output. Each oscillator's two weights integrate the
error demodulated at the line frequency (filtered-x LMS), which realizes a
resonator with poles on the unit circle: infinite loop gain exactly at `f`, so
a *stable sinusoidal* disturbance is rejected even above the loop's crossover,
at the cost of a notch in the sensitivity function only `~g*|P·S|/pi` Hz wide
instead of a broadband waterbed hit. Weights hold (no adaptation) while the
total output is railed at `clamp`, and are zeroed during `start_delay`.

- `lines` / `linesy` (float or array) (optional)
  - Line frequencies to reject, in Hz, for the x axis (element 1) / y axis
    (element 0). Up to 8 per axis. Take these from a high-resolution PSD of
    the closed-loop error; being within ~1/4 of the notch width is enough.
- `line_gain` / `line_gainy` (float or array) (optional)
  - Adaptation gain per line, in 1/s; a list shorter than `lines` extends
    with its last element, and `line_gainy` falls back to `line_gain`'s last
    element. Default 20. Sets the rejection notch width `~g*|P·S|/pi` Hz —
    wider converges faster and tolerates line-frequency wander, but the
    notch's own local waterbed amplifies neighbours. Keep gains low (≲20)
    for lines near or above crossover, and keep oscillators far enough apart
    that their notches don't overlap, or they fight each other.
- `line_phase` / `line_phasey` (float or array) (optional)
  - Demodulation phase per line, in degrees: the phase of the command→error
    path at that frequency, i.e. `∠P(f) + ∠S(f)` (plant phase from a bode
    sweep, plus the sensitivity function's phase for the gains in use). With
    a phase error the adaptation converges slower; past ~90° it goes
    unstable and *pumps* the line, so measure — don't guess. Changing `i`
    (or `iy`) moves `∠S` and thus these phases.
- `line_delay` (float) (optional)
  - Total command→error loop delay in seconds, used to default any phase not
    given in `line_phase` to the pure-delay lag `-360*f*line_delay` degrees.
    Good enough well above crossover; below ~1.5× crossover the sensitivity
    function's phase advance is significant, so give an explicit phase.
    Default 0.

Validated 2026-07-10 on the FSM loop (see `contrib/legacy/legacy_steering_profiles/steering_tuned.json`):
kills targeted stable lines 2–9× in PSD power. It cannot help broadband
disturbance, lines whose frequency wanders faster than the notch tracks, or
motion that isn't correctable from that axis's actuator (cross-axis leakage).
