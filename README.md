Anyloop plugin for fast steering mirror control
===============================================

The fast-steering stack from the FSL bench, built as one plug-in that a **stock
anyloop can load**. No patch to anyloop, no rebuild of it: build this, install
it, point a config at it.

The headline device is `anyloop:fsm`, which runs a whole steering loop —
sensor → controller → clamp → actuator — as a single pipeline stage, so a
config needs only a camera and this. The devices it composes ship alongside it
and can be used on their own:

| group | devices |
| --- | --- |
| sensors | `center_of_mass`, `wfs_com`, `fit_com` |
| controllers | `fsp` (adaptive filtered Smith predictor / LQG), `pid`, `kalman_filter`, `biquad_filter` |
| actuators | `parport_dac` (TI DAC81404 over a PCIe parallel card), `piplate_bridge`, `kdc101`, `tic_t834` |
| identification | `bode_plot`, `gain_test`, `prbs_test`, `latency_test`, `attenuation_test`, `spiral_search` |
| plumbing | `clamp`, `delay`, `file_sink`, `logger`, `matmul`, `poke`, `remove_piston`, `stop_after_count`, `test_source`, `udp_sink`, `udp_source`, `vonkarman_stream` |

Per-device documentation is in [doc/devices](doc/devices); start with
[fsm.md](doc/devices/fsm.md).


libaylp dependency
------------------

**None — and that is a deliberate difference from the other `aylp_*` plug-ins.**
They ask you to symlink `libaylp` from an anyloop checkout. That will not work
here: these devices need a libaylp *newer* than the one upstream anyloop ships
(`timing.c`, and `anyloop.h`'s measurement-validity flags), so an upstream
symlink would not compile. The `libaylp` directory in this repo is that newer
copy. Do not symlink over it.

Only `timing.c` is compiled into the plug-in. Everything else in `libaylp` is
here so the unit tests can link without an anyloop checkout; at run time those
functions come from the host.


Building
--------

```sh
meson setup build
meson compile -C build
meson test -C build      # no hardware needed
```

You need `meson`, `pkgconf`, `json-c`, `gsl` and `libusb`.


Installing
----------

```sh
sudo meson install -C build
```

That puts `aylp_fsm.so` plus one symlink per device into `<libdir>/anyloop` —
`/usr/local/lib/anyloop` on Arch, `/usr/local/lib/x86_64-linux-gnu/anyloop` on
Debian. To drop it into an anyloop install elsewhere, stage and copy:

```sh
meson install -C build --destdir /tmp/stage
sudo cp -a /tmp/stage/usr/local/lib/*/anyloop/. /path/to/their/lib/anyloop/
```

Use `cp -a`, not `cp`: the per-device names are symlinks to the one library, and
dereferencing them turns 1.5 MB into 40 MB of identical copies.


Using it
--------

Every device is addressed by the `file:` scheme, named after the symlink:

```json
{"uri": "file:/usr/local/lib/anyloop/fsp.so", "params": {"...": "..."}}
```

`aylp_fsm.so` and `fsm.so` are the same device, so the whole loop in one stage
looks like this. The `anyloop:` names inside it are resolved by this plug-in's
own registry, which is why they keep working on a host that has never heard of
`fsp`, and why they need no install path:

```json
{
	"uri": "file:/usr/local/lib/anyloop/aylp_fsm.so",
	"params": {
		"sensor": {"threshold": "auto", "track": true},
		"controller": {"delay": "auto", "broad_order": 512},
		"clamp": {"min": -0.4, "max": 0.4},
		"actuator": {"pci": "0000:05:00.2", "channel": [0, 1]}
	}
}
```

To convert a config that used built-in devices, `to_bundle_config.py` rewrites
the top-level URIs and finds the install prefix on its own:

```sh
./to_bundle_config.py old_config.json -o new_config.json
```

[contrib/conf_fsm_example.json](contrib/conf_fsm_example.json) runs on a
synthetic source, so it needs no camera and no DAC.


How it works
------------

**One library, many names.** A host derives the init symbol it looks for from
the *file name* it was given: `file:.../fsp.so` makes it `dlsym("fsp_init")`.
The device init functions are already named `fsp_init`, `center_of_mass_init`
and so on, so a symlink per device is the whole mechanism — one mapped copy of
the code, one copy of its state, thirty names.

**libaylp comes from the host at run time.** `log_impl`, `xmalloc`,
`aylp_type_from_string`, `get_contiguous_bytes`, the thread pool and the
pretty-printers are left undefined in the plug-in and resolve to the host
executable, which exports them because anyloop is linked with
`export_dynamic`. So the plug-in shares the host's log level (`anyloop -l TRACE`
reaches it) and its thread pool rather than running a second copy of each.

**Symbol scope is the trick.** The host's dynamic symbols take *precedence* over
anything `dlopen` brings in. Left alone, this plug-in's `init_device()` would
resolve to the **host's** — whose registry has never heard of `fsp` — and every
substage lookup inside `fsm` would fail, or worse, land silently on the stock
implementation of a device this replaces. Two linker flags prevent that:

- `aylp_fsm.map` restricts the dynamic symbol table to `*_init` and the timing
  API. Everything else, `init_device` included, is local.
- `-Bsymbolic-functions` binds references from inside the plug-in to its own
  definitions at link time, covering the `*_init` symbols that must stay
  visible for the host to find them.

Check with `nm -D --defined-only build/aylp_fsm.so`: 31 `*_init` symbols and 5
`aylp_timing_*`/`aylp_frame_geometry_*`, nothing else.


What the host must provide
--------------------------

Any anyloop from roughly 2026 onward. Concretely, these must be in its dynamic
symbol table (`nm -D --defined-only $(which anyloop)`):

```
alloc_check  aylp_type_from_string  aylp_units_from_string  get_contiguous_bytes
log_get_level  log_impl  pretty_matrix  pretty_matrix_uchar  pretty_vector
task_enqueue  task_runner  xcalloc  xfree_impl  xmalloc  xrealloc  xstrdup
```

`struct aylp_device`, `aylp_state` and `aylp_header` must also match. Nothing
here changed them — the `anyloop.h` additions are new enum values and one inline
helper — so the ABI is anyloop's own.


Caveats
-------

**A host built before the measurement flags can refuse to exit.** anyloop's loop
runs `while (!sigint_received && state.header.status ^ AYLP_DONE)`, which is
correct only while `AYLP_DONE` is the only flag in the word. These sensors raise
`AYLP_NO_SIGNAL`/`AYLP_BEAM_LOST`/`AYLP_FRAME_REJECTED`, and with one of those
still set the xor is nonzero, so the host ignores the device that just asked it
to stop — an infinite loop that looks exactly like a hang. Reproduced: a dark
frame plus `min_peak` plus a stock `anyloop:stop_after_count` never terminates.

Two things here fix it, and you want at least one:

- Use **this repo's** `stop_after_count.so` (also `bode_plot`,
  `spiral_search`). They call `aylp_set_done()`, which clears the measurement
  flags as it sets `AYLP_DONE`. `to_bundle_config.py` rewrites
  `anyloop:stop_after_count` for exactly this reason — don't "simplify" it back
  to the built-in.
- Keep an **`fsm`** stage in the pipeline. It clears the flags at the end of any
  iteration where `AYLP_DONE` is already set, so a host's own stop device gets
  through one iteration later instead of never.

**`fsp`'s `"delay": "auto"` needs a source that publishes loop timing.** The
channel is `libaylp/timing.h`, exported from this plug-in so a camera plug-in
can link against it and publish. A camera that does not — including
[aylp_asi](https://github.com/bnl-lab4/aylp_asi) — leaves nothing to read, and
`delay` has to be given a number.

**Hardware devices still need what they always needed.** `parport_dac` maps a
PCIe BAR and unbinds `parport_pc`, so it wants root; `kdc101`/`tic_t834` want
their USB permissions.

**The numbers in the docs are one bench's.** `doc/devices/*.md` quote measured
gains, delays and a PCI slot from the rig this was developed on. The physics
carries over; the constants do not. Those pages also cite scripts and configs
that live in the fork rather than here.


What has and has not been tested
--------------------------------

Read this before trusting it with beam time. The load path is well exercised;
the hardware path is not exercised at all.

**Verified**, against a clean build of anyloop `6d801f8` driven from a converted
config:

- the plug-in loads, and the symbol scoping holds — nothing leaks to the host's
  `init_device` or to its own devices;
- `fsm` resolves `anyloop:fsp` on a host with no such device, getting as far as
  fsp's own parameter validation;
- the shadowing devices really are these: `center_of_mass` accepted `track`,
  `pid` accepted `park_after`, `test_source` accepted `kind: noise` — all
  parameters a stock device rejects as unknown;
- park-on-exit fires through the plug-in;
- the `AYLP_DONE` hazard above reproduces, and both mitigations clear it;
- install layout, symlinks, RUNPATH and `DT_NEEDED`;
- the cross-plug-in timing channel: a camera linked against the equivalent
  library in the fork bound `aylp_timing_publish` to the copy `fsp` reads,
  confirmed with `LD_DEBUG=bindings`. Same mechanism, but this exact pair of
  files has not been retested;
- `meson test` — eight suites, no hardware.

**Not tested. Do not assume these work:**

1. **No hardware device has ever been run through the plug-in.**
   `parport_dac` wants root, the PCIe card and a mirror that will physically
   move; the camera path is proven to load and bind but has never grabbed a
   frame through it; `kdc101`, `tic_t834` and `piplate_bridge` have not been
   exercised at all.
2. **`fsp` has never completed init through the plug-in** — it was fed a
   deliberately incomplete parameter set. Its arithmetic has no coupling to the
   host, so the risk is low, but low is not zero.
3. **The identification devices** — `gain_test`, `prbs_test`, `bode_plot`,
   `attenuation_test`, `latency_test`, `spiral_search` — have never been loaded
   through the plug-in, though their unit tests pass.
4. **One host commit only.** anyloop `6d801f8` exactly. Older installs need
   `get_contiguous_bytes`, the thread pool and the type/units helpers to already
   be present; nothing earlier has been tried.
5. **One machine only.** x86-64, `dlopen@GLIBC_2.34`, built against the system
   gsl/json-c/libusb. Elsewhere, build it there.

Closing 1 through 3 takes one bench pass: run a converted profile against a
stock anyloop with the camera and the DAC live.


Provenance
----------

Extracted from the FSL bench fork of anyloop
([leites-lev/anyloop](https://github.com/leites-lev/anyloop)) at `836a58d`,
which is where the devices are developed against real hardware. This repo is a
snapshot for deployment, not a second line of development: fixes belong
upstream in the fork first.

Licensed under the Mozilla Public License 2.0, like anyloop itself.
