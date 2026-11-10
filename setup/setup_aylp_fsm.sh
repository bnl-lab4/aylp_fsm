#!/usr/bin/env bash
# Ubuntu/Debian, Intel/AMD 64-bit. Run as your normal user:
#   bash setup_aylp_fsm.sh
# No cleanup deletion, no DAC commands, and no changes to your existing fork.
set -Eeuo pipefail
trap 'echo "Setup stopped at line $LINENO. Fix the error above and rerun this script." >&2' ERR
[[ $(uname -s) == Linux && $(uname -m) == x86_64 && $EUID != 0 ]] || {
  echo 'Use 64-bit Intel/AMD Linux and run as your normal user.'; exit 1;
}

# 1. Install the compiler, plugin dependencies, and calibration/report libraries.
sudo apt-get update
sudo apt-get install -y git build-essential meson ninja-build pkg-config \
  libjson-c-dev libgsl-dev libusb-1.0-0-dev python3 python3-numpy \
  python3-scipy python3-matplotlib poppler-utils pciutils

# 2. Keep sources and installed programs together in a permanent, private folder.
# ZWO's library is loaded from its source folder: keep this folder after setup.
base="$HOME/aylp-fsm-system"
prefix="$base/install"
plugins="$prefix/lib/anyloop"
mkdir -p "$base"

# Fetch fixed versions for repeatable setup. Existing checkouts are never reset.
checkout() {
  local name=$1 revision=$2 dir="$base/$1" fresh=0
  if [[ ! -e "$dir" ]]; then
    git clone "https://github.com/bnl-lab4/$name.git" "$dir"
    fresh=1
  fi
  [[ -d "$dir/.git" ]] || { echo "$dir is not a Git checkout."; exit 1; }
  [[ -z $(git -C "$dir" status --porcelain --untracked-files=no) ]] || {
    echo "Edited source files in $dir; existing files were preserved."; exit 1;
  }
  if [[ $(git -C "$dir" rev-parse HEAD) != "$revision" ]]; then
    # Only a newly cloned, untouched folder may be switched to the pinned revision.
    if [[ $fresh == 0 ]]; then
      echo "Unexpected source version in $dir; existing files were preserved."; exit 1
    fi
    git -C "$dir" checkout --detach "$revision"
  fi
  printf '%s\n' "$revision" > "$dir/.setup-version"
}
checkout anyloop 6d801f807114f357c2115454cc97590732b0f1d0
checkout aylp_fsm ceed616663d6eb332fd2c7c2be0346095fb1dc3a
checkout aylp_asi 376daded09e3d7a7e9606ea5044ff689da9e7834

# 3. Only the CAMERA plugin needs headers linked from the Anyloop host.
# Keep aylp_fsm's own libaylp directory intact, as required by its README.
if [[ ! -e "$base/aylp_asi/libaylp" && ! -L "$base/aylp_asi/libaylp" ]]; then
  ln -s "$base/anyloop/libaylp" "$base/aylp_asi/libaylp"
fi
[[ $(readlink -f "$base/aylp_asi/libaylp") == "$base/anyloop/libaylp" ]] || {
  echo 'The camera libaylp link points elsewhere; existing files were preserved.'; exit 1;
}

# 4. Build and install into the private folder. Reruns reuse each build directory.
build_project() {
  local dir="$base/$1"
  if [[ -f "$dir/build/meson-private/coredata.dat" ]]; then
    meson setup --reconfigure "$dir/build" "$dir" --prefix "$prefix" --libdir lib
  else
    meson setup "$dir/build" "$dir" --prefix "$prefix" --libdir lib
  fi
  meson compile -C "$dir/build"
  meson test -C "$dir/build" --print-errorlogs
  meson install -C "$dir/build"
}
build_project anyloop
build_project aylp_fsm
build_project aylp_asi
test -f "$plugins/aylp_fsm.so"
test -f "$plugins/aylp_asi.so"

# 5. Install the camera USB rules. Its ZWO library is already bundled in aylp_asi.
sudo install -m 644 "$base/aylp_asi/libasi/lib/asi.rules" /etc/udev/rules.d/99-asi.rules
sudo udevadm control --reload-rules
python3 -c 'import ctypes, sys; ctypes.CDLL(sys.argv[1]); print("ZWO library loads successfully.")' \
  "$base/aylp_asi/libasi/lib/x64/libASICamera2.so"

# 6. Convert and run the synthetic example: tests plugin loading, with NO hardware.
python3 "$base/aylp_fsm/to_bundle_config.py" --prefix "$plugins" \
  "$base/aylp_fsm/contrib/conf_fsm_example.json" -o "$base/smoke.json"
timeout 20 "$prefix/bin/anyloop" "$base/smoke.json"

# 7. Save a launcher so you always use this installation, not another Anyloop.
printf '#!/usr/bin/env bash\nexec %q "$@"\n' "$prefix/bin/anyloop" > "$base/run-anyloop.sh"
chmod +x "$base/run-anyloop.sh"
echo "Setup and synthetic test passed. Launcher: $base/run-anyloop.sh"
echo "Plugins: $plugins"
echo 'Reconnect the camera to apply USB rules. No DAC output was commanded.'
echo 'For live hardware: configure the camera and PCI address, then calibrate.'
echo 'aylp_asi does not publish timing: use calibrated numerical fs/delay settings.'
if [[ -r /sys/kernel/security/lockdown ]]; then
  cat /sys/kernel/security/lockdown
  echo 'The direct PCI DAC backend requires [none] and administrator access.'
fi
