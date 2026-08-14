anyloop:udp_sink
================

Types and units: `[T_ANY, U_ANY] -> [T_UNCHANGED, U_UNCHANGED]`.

This device writes the current pipeline state to a UDP port as an AYLP file. See
[filetype.md](../filetype.md) for documentation on the AYLP file format.

Parameters
----------

- `ip` (string) (required)
  - The IP address to send the data to.
- `port` (string) (required)
  - The port to send the data to.
- `decimation` (integer) (optional)
  - Send only every `decimation`'th iteration; skipped iterations return
    immediately without copying or writing anything. Defaults to 1 (send every
    iteration).
  - Use this when the loop runs much faster than whatever is reading the port.
    A viewer that redraws at tens of Hz cannot consume a multi-kHz frame stream;
    without decimation the kernel drops the excess anyway, but the loop still
    pays for the copy and the `writev` on every iteration. Setting
    `decimation` to roughly `loop_rate / viewer_rate` keeps the traffic (and the
    syscall) off the hot path.
- `reduce` (string) (optional)
  - What to send: `"none"` (the default) sends the data itself; `"stats"` sends
    `[peak, mean, saturated_fraction]` as three doubles instead, typed as a
    `T_VECTOR` in `U_COUNTS` with `log_dim` 3x1.
  - This exists so a viewer can follow a scalar summary of a large block at the
    *full loop rate*. A 248x248 uchar frame is 61.5 kB, which is 233 MB/s at
    3788 Hz and can be neither sent nor decoded every iteration, so a raw frame
    sink has to be decimated and its trace updates at `loop_rate / decimation`.
    The reduction is 24 bytes, i.e. ~91 kB/s at 3788 Hz, so it needs no
    decimation at all.
  - It is not free: the reduction is one pass over the data. Measured on a
    248x248 uchar frame, `udp_sink` proc goes from 5.5 us (raw, no reduction) to
    16.1 us. That is ~0.6% of a 2.5 ms period but ~6% of a 264 us one, so weigh
    it before adding this to a config whose numbers have to stay comparable.
  - `saturated_fraction` is the share of samples at the maximum the sample type
    can represent. That is only meaningful for uchar (camera) data; it is
    reported as 0 for double data, which has no such ceiling. It matters because
    a peak pinned at 255 hides everything above it, leaving `mean` as the only
    honest brightness reading.
  - The reduction is lossy and there is no way back to the samples from it. Use
    a raw sink when you want the data itself.
  - Note that an undecimated sink with no listener bound will log a send error
    (`Connection refused`) on *every* iteration, since the socket is
    `connect()`ed. That is loud enough to matter for a long unattended run --
    start the receiver first, or decimate.

