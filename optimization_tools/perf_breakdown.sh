#!/bin/bash
# Break down a saved perf run by function. Usage:
#   optimization_tools/perf_breakdown.sh <run label, e.g. attempt2-perf> <job prefix> <binary label> <match substring>
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
R=$root/build/ultimate-runtime-perf/squashfs-root
J=$(ls -d "$root"/build/optimization-runs/$1/$2*)
docker run --rm --privileged -v "$J":/job -v "$root/optimization_tools":/tools:ro -v "$root/build/optimization-builds/$3":/binary:ro \
  --mount type=bind,source="$R",target=/opt/orca-ro,readonly --mount type=bind,source="$R",target=/run/host_virtiofs"$R",readonly \
  --entrypoint sh mirdain-orca-profiler:opt -c "mkdir -p /opt/orca-local/bin && cp /binary/orca-slicer /opt/orca-local/bin/ &&
    perf-any script -i /job/perf.data --no-inline -F comm,tid,ip,sym 2>/dev/null | python3 /tools/perf_callers.py '$4'"
