#!/bin/bash
# Old vs new web slicer timings across meshes; see web_benchmark.py. Writes build/web-benchmark/results.jsonl.
#   A head_old: committed server code + old Orca runtime (ef8d293)
#   B work_old: working-tree server code + old Orca runtime
#   C work_new: working-tree server code + Orca runtime in the backend image
# Needs build/web-slicer-runtime/{orca-ef8d293,server-head} (see OPTIMIZATION_PROGRESS.md).
set -euo pipefail
fork=/Users/ryanpineda/Fork
rt=$fork/MirdainSlicer/build/web-slicer-runtime
out=$fork/MirdainSlicer/build/web-benchmark
report=0bb7f9560f8c46adaa8a6b0e4eb34529   # 3DBenchy, unified mode, 10 degrees
meshes="template"
for m in 0709 0e5d 3e7a 4557 50ae 5e0b 7a3f bd62 d199 f4c5; do
  meshes+=" /meshes/$(ls $fork/experiments/uploaded-3mfs | grep "^$m" | head -1)"
done
warmup=/meshes/$(ls $fork/experiments/uploaded-3mfs | grep "^a63d" | head -1)
run() {  # label server_dir runtime_dir(or empty)
  local args=(-v "$fork/MirdainSlicer/optimization_tools:/bench:ro" -v "$out:/out"
              -v "$fork/experiments/uploaded-3mfs:/meshes:ro" -v "$2:/app/server:ro")
  [ -n "$3" ] && args+=(-v "$3:/opt/orca:ro")
  docker compose -f "$fork/mirdain_web_slicer/compose.yaml" run --rm -T "${args[@]}" backend \
    sh -c "python /bench/web_benchmark.py warmup /out/warmup.jsonl $report $warmup >/dev/null && \
           python /bench/web_benchmark.py $1 /out/results.jsonl $report $meshes"
}
mkdir -p "$out"
for pass in 1 2; do
  run head_old "$rt/server-head/server" "$rt/orca-ef8d293"
  run work_old "$fork/mirdain_web_slicer/server" "$rt/orca-ef8d293"
  run work_new "$fork/mirdain_web_slicer/server" ""
done
