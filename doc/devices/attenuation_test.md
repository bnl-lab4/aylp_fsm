anyloop:attenuation_test
========================

Types and units: `[T_VECTOR, U_ANY] -> [T_UNCHANGED, U_UNCHANGED]`.

This device measures the closed-loop attenuation of a feedback loop by
comparing error spectra with the loop open and closed, in one automated run.
Place it on the error signal, *before* the controller (and before any
predictor like `anyloop:kalman_filter`, so the raw error is what gets
recorded).

The run has four phases, timed from the first loop iteration:

1. **start delay** (`start_delay` seconds): the device zeroes the vector every
   iteration and records nothing. Use this to cover other devices'
   `start_delay` holds and sensor acquisition.
2. **open** (`open_time` seconds): the device records the incoming error, then
   zeroes the vector before passing it on. A zero error keeps an integrating
   controller (and its line oscillators) parked at zero output, so the
   actuator holds its bias and the loop is effectively open — while every
   downstream device still sees a normal frame cadence.
3. **settle** (`settle_time` seconds): the vector passes through untouched —
   the loop is now closed — but nothing is recorded, giving the controller
   time to acquire and adaptive filters time to re-converge.
4. **closed** (`closed_time` seconds): the vector passes through and is
   recorded.

When the closed phase ends, the device computes a Welch PSD (Hann window, 50 %
overlap, per-segment mean removal) of each vector element for both phases,
using each phase's sample rate as measured from its own timestamps. The two
rates may legitimately differ (closing the loop adds actuator writes and
usually slows the iteration), so the open PSD is re-gridded onto the closed
phase's frequency axis before comparison. Two details keep the dB values
exact rather than approximate: the open segments span the same *time* as a
closed segment (so both windows have the same resolution bandwidth in Hz and
an off-bin line scallops identically in both PSDs, cancelling in the ratio),
and the open FFTs are zero-padded ≥ 8× before linear interpolation (so a
narrow line's peak is not flattened between grid points). Verified against a
synthetic line with known attenuation: ≤ 0.005 dB error, vs up to 0.26 dB
with plain whole-bin interpolation. Only if the closed phase somehow ran
*faster* are bins above the open Nyquist clamped (with a warning).
It writes the results next to the PDF as a `.dat` file — columns: frequency
(closed-phase axis), then `psd_open`, `psd_closed`, `atten_dB` per element,
where `atten_dB = 10*log10(psd_closed/psd_open)` (negative = attenuated,
positive = amplified/waterbed) — then shells out to python3 (numpy +
matplotlib required) to render the PDF: per element, the open/closed PSD
overlay on top, the attenuation curve in the middle with the strongest
attenuated and amplified frequencies annotated, and at the bottom the error
variance contributed by each 5 Hz band (the PSD integrated over the band, open
vs closed), with the top closed-loop bands tagged with their share of the
total closed-loop error variance. A second page reports the run's start date
and time, the `config` summary, the mean (pointing offset about the setpoint)
and RMS-about-mean (jitter) for each element open vs closed, and a table of the
RMS each 10 Hz band contributes from 0-300 Hz — all scaled to pixels by
`pixel_scale`. Finally it raises `AYLP_DONE` to stop the loop.

Parameters
----------

- `open_time` (float) (required)
  - Seconds of open-loop error to record.
- `closed_time` (float) (required)
  - Seconds of closed-loop error to record.
- `start_delay` (float) (optional)
  - Seconds to hold the loop open before the open recording starts.
    [default: 0]
- `settle_time` (float) (optional)
  - Seconds to let the loop settle after closing, before the closed recording
    starts. Make this generous: recording during acquisition corrupts the
    closed spectrum. [default: 5]
- `nfft` (integer) (optional)
  - Welch segment length; must be a power of 2 and >= 64. Frequency resolution
    is `sample_rate/nfft`. Automatically halved (with a warning) if a
    recording turns out shorter than one segment. [default: 4096]
- `labels` (string) (optional)
  - Comma-separated names for the vector elements, used in the plot titles
    (e.g. `"y,x"` for a center-of-mass error). [default: element indices]
- `pixel_scale` (float) (optional)
  - Multiplier converting one output unit to pixels, used only for the report
    on the PDF's second page (start time, per-band RMS table, and per-axis
    mean/jitter). For an `anyloop:center_of_mass` source, which emits
    `-1 + 2*px/(dim-1)` across the frame, this is `(dim-1)/2` (e.g. `31.5` at
    64x64). Leave at the default to report in raw units. [default: 1]
- `config` (string) (optional)
  - Free-text summary of the loop configuration, printed verbatim on the PDF's
    second page. Newlines (`\n`) start new lines. Use it to record the pid /
    kalman / line-rejection settings the run was taken under. [default: none]
- `output_file` (string) (optional)
  - PDF filename; the `.dat` file takes the same name with the extension
    replaced. [default: "attenuation.pdf"]

See `contrib/legacy/outdated_scripts/attenuation_test.json` for a full example built around the tuned
tip/tilt steering loop.
