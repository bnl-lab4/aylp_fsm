#!/usr/bin/env python3
"""Point a config at the drop-in bundle instead of at built-in devices.

A stock upstream anyloop has never heard of `anyloop:fsp`, `anyloop:fit_com` or
`anyloop:parport_dac`, and its `anyloop:center_of_mass` is the version from
before this fork. This rewrites the pipeline's device URIs to load them out of
the bundle instead:

    "uri": "anyloop:fsp"   ->   "uri": "file:/usr/local/lib/anyloop/fsp.so"

Only TOP-LEVEL pipeline entries are rewritten. URIs nested inside an
`anyloop:fsm` stage are deliberately left alone: fsm resolves its substages
through the bundle's own registry, which is the same set of devices, so
rewriting them would only make the config longer and pin it to an install path.

usage:
    to_bundle_config.py run_1000hz_20pct.json > run_dropin.json
    to_bundle_config.py --prefix ~/.local/lib/anyloop conf.json -o out.json
"""

import argparse
import json
import sys
from pathlib import Path

# Devices the bundle provides, i.e. the names that exist as <name>.so symlinks
# next to aylp_fsm.so. Keep in sync with `bundle_devices` in meson.build.
BUNDLE_DEVICES = {
    "center_of_mass", "wfs_com", "fit_com", "kdc101", "tic_t834", "clamp",
    "biquad_filter", "kalman_filter", "delay", "file_sink", "logger",
    "matmul", "pid", "fsp", "fsm", "poke", "remove_piston",
    "stop_after_count", "test_source", "udp_sink", "udp_source",
    "vonkarman_stream", "piplate_bridge", "parport_dac", "bode_plot",
    "latency_test", "gain_test", "prbs_test", "attenuation_test",
    "spiral_search",
}

# Plug-ins that are built beside the bundle rather than inside it (they pull in
# an external SDK). A config naming one by some other path gets re-pointed at
# the install prefix, since that is the copy linked against the bundle's timing
# channel.
SIDECAR_PLUGINS = {"asi_source", "aylp_alsa"}


# Where the bundle usually lands. meson's libdir is not the same everywhere --
# lib on Arch, lib/x86_64-linux-gnu on Debian -- so look rather than guess.
PREFIX_GLOBS = [
    "/usr/local/lib/anyloop", "/usr/local/lib/*/anyloop",
    "/usr/lib/anyloop", "/usr/lib/*/anyloop",
]


def find_prefix():
    """Find an installed bundle, or exit telling the caller to name one."""
    for pattern in PREFIX_GLOBS:
        for path in sorted(Path("/").glob(pattern.lstrip("/"))):
            if (path / "aylp_fsm.so").exists():
                return str(path)
    sys.exit("No installed aylp_fsm.so found in " + ", ".join(PREFIX_GLOBS)
        + "\nPass --prefix with the directory you installed the bundle into.")


def rewrite_uri(uri, prefix):
    """Return (new_uri, note) for one top-level pipeline URI."""
    if uri.startswith("anyloop:"):
        name = uri[len("anyloop:"):]
        if name in BUNDLE_DEVICES:
            return f"file:{prefix}/{name}.so", f"{uri} -> bundle"
        return uri, None  # a device the bundle doesn't carry; leave it
    if uri.startswith("file:"):
        path = Path(uri[len("file:"):])
        name = path.name.split(".")[0]
        if name in SIDECAR_PLUGINS and path.parent != Path(prefix):
            return f"file:{prefix}/{path.name}", f"{uri} -> {prefix}"
        return uri, None
    return uri, None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config", help="config file to convert")
    ap.add_argument("-p", "--prefix",
        help="where the bundle is installed (default: look for it)")
    ap.add_argument("-o", "--output", help="write here instead of stdout")
    args = ap.parse_args()

    prefix = (args.prefix or find_prefix()).rstrip("/")
    conf = json.loads(Path(args.config).read_text())
    pipeline = conf.get("pipeline")
    if not isinstance(pipeline, list):
        sys.exit(f"{args.config}: no pipeline array")

    notes = []
    for entry in pipeline:
        uri = entry.get("uri")
        if not isinstance(uri, str):
            continue
        new_uri, note = rewrite_uri(uri, prefix)
        entry["uri"] = new_uri
        if note:
            notes.append(note)

    text = json.dumps(conf, indent="\t") + "\n"
    if args.output:
        Path(args.output).write_text(text)
    else:
        sys.stdout.write(text)

    for note in notes:
        print(note, file=sys.stderr)
    print(f"{len(notes)} of {len(pipeline)} stage(s) rewritten",
        file=sys.stderr)


if __name__ == "__main__":
    main()
