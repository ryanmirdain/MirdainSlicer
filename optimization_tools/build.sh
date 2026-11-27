#!/bin/bash
# Incremental RelWithDebInfo (-O2 -g) build of orca-slicer in the Linux builder container,
# then copy the binary to build/optimization-builds/<label>/orca-slicer.
# Usage: optimization_tools/build.sh <label> [jobs]
set -euo pipefail
label=$1; jobs=${2:-6}
root=$(cd "$(dirname "$0")/.." && pwd)
build() {
  docker run --rm --name "mirdain-orca-opt-build-$$" --user "$(id -u):$(id -g)" \
    --mount type=bind,src="$root",dst=/__w/OrcaSlicer/OrcaSlicer -w /__w/OrcaSlicer/OrcaSlicer \
    mirdain-orca-builder:opt cmake --build build --config RelWithDebInfo --target orca-slicer --parallel "$jobs"
}
# Changed C++ sources whose RelWithDebInfo object is older than the source. Ninja in the
# container once saw a stale mtime through the virtiofs share right after an edit and skipped
# the file, producing a binary without the change.
stale() {
  cd "$root"
  git diff --name-only -- 'src/*.cpp' | while read -r f; do
    o=$(find build/src -path '*RelWithDebInfo*' -name "$(basename "$f").o" -newer "$f" | head -1)
    [ -z "$o" ] && echo "$f"
  done
}
build
if [ -n "$(stale)" ]; then
  echo "Stale objects after build, rebuilding: $(stale)"
  sleep 2; build
  if [ -n "$(stale)" ]; then echo "Still stale: $(stale)" >&2; exit 1; fi
fi
mkdir -p "$root/build/optimization-builds/$label"
cp "$root/build/src/RelWithDebInfo/orca-slicer" "$root/build/optimization-builds/$label/orca-slicer"
echo "$root/build/optimization-builds/$label/orca-slicer"
