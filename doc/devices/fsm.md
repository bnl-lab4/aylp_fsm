anyloop:fsm
===========

Types and units: whatever its substages add up to — see [Typechecking](#typechecking).
For the usual steering stack that is
`[T_MATRIX_UCHAR, U_COUNTS] -> [T_VECTOR, U_V]`.

The whole fast steering mirror loop as one pipeline stage: sensor → controller →
clamp → actuator, each an ordinary anyloop device that this one owns and runs.
A config that spelled those out as four pipeline entries becomes a camera and
this:

```json
{
	"pipeline": [
		{
			"uri": "file:build/asi_source.so",
			"params": {"camera_name": "ASI290MM", "exposure": 1000}
		},
		{
			"uri": "anyloop:fsm",
			"params": {
				"sensor": {
					"threshold": "auto",
					"track": true,
					"registration": true
				},
				"controller": {
					"type": "vector", "units": "V",
					"delay": "auto", "fs": "auto",
					"broad_order": 512, "broad_mu": 0.03,
					"x": {"K": 0.19}, "y": {"K": 0.24}
				},
				"clamp": {"min": -0.4, "max": 0.4},
				"actuator": {
					"pci": "0000:05:00.2",
					"channel": [0, 1], "index": [0, 1],
					"scale": [0.9, -0.9], "range": "+-5"
				}
			}
		}
	]
}
```

The substages are the real devices — `center_of_mass`, `fsp`, `clamp`,
`parport_dac` here — initialized through the same `init_device()` the top level
uses, with their own documentation and their own parameters. This device adds
three things and no signal processing of its own:

 1. the substages are typechecked against each other at init, and the composite
    publishes the resulting types outward, so it still slots into a pipeline
    that is itself typechecked;
 2. the mirror is **parked** when the loop exits, including on Ctrl-C and on a
    fatal error;
 3. an error out of a substage that was supposed to change the pipeline type
    stops the loop instead of being swallowed.


Parameters
----------

- `sensor`, `controller`, `clamp`, `actuator` (see below) (at least one required)
  - The substages, run in **that order** regardless of the order the keys appear
    in the json object. A role that is absent is simply not run: an open-loop
    recording is `sensor` alone, a replay is `controller` and `actuator`.
- `park_on_exit` (boolean) (optional)
  - Whether to write one park command through the actuator when the loop exits.
    Defaults to `true`. Ignored (with a warning) if there is no actuator stage.
- `park_value` (float) (optional)
  - The command every element is parked at. Defaults to `0`, which puts each
    actuator channel at its own configured `offset` — the bias the plant gain
    `K` was measured at, not necessarily zero volts.

### The value of a role

Each role takes one of four forms:

| form | meaning |
| --- | --- |
| `"fit_com"` | that device, no params. A bare name is short for the `anyloop:` scheme; anything with a scheme (`file:/opt/foo/bar.so`) is passed through as-is, so plug-in devices work here too. |
| `{"uri": ..., "params": {...}}` | that device with those params |
| `{...}` (no `uri` key) | params for this role's default device |
| `[...]` | an array of any of the above, run in order |

The default device for each role is the one this bench runs:

| role | default |
| --- | --- |
| `sensor` | `anyloop:center_of_mass` |
| `controller` | `anyloop:fsp` |
| `clamp` | `anyloop:clamp` |
| `actuator` | `anyloop:parport_dac` |

`clamp` is given `"type": "vector"` unless you say otherwise, since inside here
the thing being clamped is always the command vector.

### Arrays are where the taps go

Folding the loop into one stage would otherwise cost you the `udp_sink` and
`file_sink` entries that sit *between* the stages — the live
`watch_error_ts.jl` view, the error recordings. Put them in the array for the
role they follow. This records and publishes the measurement, then controls on
it:

```json
"sensor": [
	{"uri": "anyloop:center_of_mass", "params": {"threshold": "auto"}},
	{"uri": "anyloop:file_sink", "params": {"filename": "data/err.aylp"}},
	{"uri": "anyloop:udp_sink", "params": {"ip": "127.0.0.1", "port": "64732"}}
]
```

A sink in the `sensor` array sees the measured error; one in the `controller` or
`clamp` array sees the command. Note that a tap does **not** see the park write
— only the actuator does.


Parking
-------

When the loop ends, before anything is closed, one command of `park_value` is
written through the last `actuator` stage. That is the whole mechanism: the
actuator's own `scale`/`offset` decide what voltage a zero command is, so the
mirror goes to the same place `parport_dac`'s `start_delay` hold puts it at
startup.

This matters because `parport_dac_fini()` deliberately leaves its outputs where
they are — the mirror is presumably still holding a beam — and parking has so
far been a separate `anyloop` run of
`contrib/calibration-scripts/configurations/conf_parport_dac_park.json`, which
`run_steering.py` fires *before* the loop. That covers a planned run and not the
case that actually strands the mirror at a rail: a Ctrl-C or a fatal error part
way through one. `fsm` parks on those too.

The park write is skipped, with a warning, if no command vector ever reached the
actuator (an init failure, or a run that ended before the first iteration
completed). Set `park_on_exit: false` if you want the mirror left exactly where
the last command put it — for example while measuring a settling transient.


Typechecking
------------

At init the substages' declared types are chained together exactly as
`anyloop.c` chains the top-level pipeline, and the result becomes this device's
own `type_in`/`units_in`/`type_out`/`units_out`. While the running type is still
`AYLP_T_UNCHANGED` the data in flight is whatever the composite was handed, so
those stages' requirements become the composite's requirements; once a stage
declares a concrete output type the rest of the chain is checked here and now.

So `center_of_mass` → `fsp` → `clamp` → `parport_dac` reports
`[T_MATRIX_UCHAR, U_COUNTS] -> [T_VECTOR, U_V]` outward, and a chain that cannot
work — a vector `clamp` behind a `matrix_uchar` source, say — is refused at
startup with the same message the top level would have given, rather than at
3 kHz.


Errors
------

A substage returning nonzero is treated the way `anyloop.c` treats a top-level
device: recoverable (logged, loop continues) unless the stage was supposed to
change the pipeline type, in which case everything downstream of it — the
command about to go on the mirror — is garbage.

At the top level anyloop exits outright in that case. Here the mirror is holding
a beam, so instead the loop is asked to stop cleanly (`AYLP_DONE`), the
remaining stages are skipped for that iteration so the bad command never
reaches the actuator, and the park write happens on the way out.


What this does not do
---------------------

- **No beam-loss policy of its own.** The sensors raise `AYLP_BEAM_LOST` and
  `AYLP_NO_SIGNAL`; `fsp`'s gap handling and `pid`'s `park_after` already act on
  them. A second, outer policy here would fight them.
- **No per-stage profiling.** `anyloop -p` profiles top-level devices, so the
  composite shows up as one row. Split the stages back out into the pipeline for
  a session where you care which one costs what — the parameters are identical
  either way.


On a stock anyloop
------------------

`fsm` is the most useful entry point in this plug-in
([../../README.md](../../README.md)): loaded as
`file:/usr/local/lib/anyloop/fsm.so` it resolves its substages through the
plug-in's own registry, so the `anyloop:` names in its params keep working even
on a host whose built-in registry has never heard of `fsp` or `parport_dac`, and
the config stays independent of the install path. It also clears the
measurement-validity flags on any iteration where `AYLP_DONE` is already set,
which is what keeps a pre-flags host from refusing to exit — see the caveats in
that README.
