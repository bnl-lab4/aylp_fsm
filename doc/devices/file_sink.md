anyloop:file_sink
=================

Types and units: `[T_ANY, U_ANY] -> [T_UNCHANGED, U_UNCHANGED]`.

This device writes the current pipeline state to an AYLP file. See
[filetype.md](../filetype.md) for documentation on the AYLP file format.

Parameters
----------

- `filename` (string) (required)
  - The filename to write the pipeline data to. An existing recording under this
    name is preserved rather than overwritten — see "Runs are never
    overwritten" below.
- `flush` (boolean) (optional)
  - Whether or not to flush the output every iteration. Setting this may
    slightly hinder performance, but will ensure that anything reading the file
    is not waiting for a buffered write.

Runs are never overwritten
--------------------------

This device used to open its file with `"wb"`, which truncates. Every start
silently destroyed the previous run's recording, and nothing said so — you found
out by going to look for a run that was no longer there. Losing a measurement
that took minutes of beam time to a restart is not a recoverable mistake.

Now, if the configured filename already holds a non-empty regular file, that file
is renamed to `base.NNNN.ext` — the lowest free four-digit sequence number,
inserted before the extension — and a fresh file is opened under the configured
name. So a directory accumulates:

```
frames.aylp        <- the run happening now
frames.0001.aylp   <- the first run
frames.0002.aylp   <- the second
```

The old run is moved aside rather than the new run being diverted to a new name.
That ordering is deliberate: the configured filename always holds the *current*
run, so anything that reads a recording back by the name in the config —
`contrib/calibration-scripts/tools/find_roi.py`, the analysis scripts — keeps working and keeps seeing the
newest data, while every earlier run survives beside it.

An existing but *empty* file is reused in place, so runs that fail during startup
do not litter the directory with numbered nothings. The rename is logged at INFO,
and if it fails (a read-only directory, say) the run continues with a warning
rather than aborting — losing the old file is bad, refusing to take the
measurement in front of you is worse.

Note that nothing here deletes anything. Long recordings at loop rate are large
(a 248×248 frame is 61.5 kB, so 400 Hz is 88 GB/hour); pruning old sequence
numbers is left to you, on purpose.

