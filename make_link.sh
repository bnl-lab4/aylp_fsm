#!/bin/sh
# Point a per-device plug-in name at the bundle.
#
# A stock anyloop derives the init symbol it looks for from the FILE NAME it was
# given: `file:/opt/anyloop/fsp.so` makes it dlsym("fsp_init"). So every device
# in the bundle needs its own name on disk, and the cheapest way to give it one
# is a symlink -- one mapped copy of the library, one copy of its state, N names.
#
# usage: make_link.sh <bundle.so> <link.so>
set -e
ln -sf "$(basename "$1")" "$2"
