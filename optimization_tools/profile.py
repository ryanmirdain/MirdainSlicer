"""Profile a cached nonplanar Orca job against a given orca-slicer binary.

    python3 optimization_tools/profile.py JOB_NAME BINARY_LABEL [--repeat N] [--perf]

JOB_NAME is a folder in optimization_jobs/ (a prefix is enough); BINARY_LABEL is a folder in
build/optimization-builds/ (made by build.sh). The Orca arguments stored in the job's
section-01/command.json are replayed in the privileged mirdain-orca-profiler:opt image with the
ultimate-runtime-perf AppImage runtime mounted at /opt/orca and the chosen binary bind-mounted over
/opt/orca/bin/orca-slicer.

Timing runs (--repeat N, default 1) register uprobes on the phase functions (container_run.py
PHASES) and report, per phase, the summed inclusive wall time of its calls (min over repeats).
--perf adds one extra sampling run (perf, dwarf call graphs) and writes perf-report.txt.
Outputs: build/optimization-runs/<label>/<job>/ (G-code, logs, probes.json, perf-report.txt).
A summary line is appended to build/optimization-runs/results.jsonl.
"""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
JOBS = ROOT / 'optimization_jobs'
BUILDS = ROOT / 'build' / 'optimization-builds'
RUNS = ROOT / 'build' / 'optimization-runs'
RUNTIME = ROOT / 'build' / 'ultimate-runtime-perf' / 'squashfs-root'
IMAGE = 'mirdain-orca-profiler:opt'


def gcode_hash(path):
    """Hash of the G-code without comment lines (header, timestamps, metadata)."""
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for line in stream:
            if not line.startswith(b';'):
                digest.update(line)
    return digest.hexdigest()


def run(job, label, mode, out=None, local_io=False):
    target = RUNS / (out or label) / job
    shutil.rmtree(target, ignore_errors=True)
    shutil.copytree(JOBS / job, target, ignore=lambda folder, names: [n for n in names if n == 'source.txt'])
    section = target / 'section-01'
    (section / 'config').mkdir(exist_ok=True)
    saved = json.loads((section / 'command.json').read_text())
    arguments = saved[next(i for i, a in enumerate(saved) if a.startswith('mirdain_web_slicer-backend')) + 1:]
    command = ['docker', 'run', '--rm', '--init', '--name', f'orca-opt-{label}-{job[:8]}-{time.time_ns()}',
               '--network', 'none', '--user', '0:0', '--privileged', *(['-e', 'LOCAL_IO=1'] if local_io else []), '--workdir', '/job/section-01',
               '--mount', f'type=bind,source={target},target=/job',
               '--mount', f'type=bind,source={RUNTIME},target=/opt/orca-ro,readonly',
               '--mount', f'type=bind,source={BUILDS / label},target=/binary,readonly',
               '--mount', f'type=bind,source={RUNTIME},target=/run/host_virtiofs{RUNTIME},readonly',
               '--mount', f'type=bind,source={ROOT / "optimization_tools"},target=/tools,readonly',
               '--entrypoint', 'python3', IMAGE, '/tools/container_run.py', mode, '--',
               'xvfb-run', '-a', '/opt/orca-local/AppRun', *arguments]
    started = time.perf_counter()
    with (section / 'orca.log').open('w') as log:
        code = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT).returncode
    wall = time.perf_counter() - started
    gcode = section / 'plate_1.gcode'
    result = {'exit_code': code, 'container_wall': round(wall, 3),
              'gcode_sha256': gcode_hash(gcode) if gcode.exists() else None}
    if mode == 'probes' and (target / 'probes.json').exists():
        probes = json.loads((target / 'probes.json').read_text())
        phases = probes['phases']
        firsts = [p['first'] for p in phases.values() if p['first'] is not None]
        lasts = [p['last'] for p in phases.values() if p['last'] is not None]
        result['phases'] = {k: round(v['inclusive'], 3) for k, v in sorted(phases.items())}
        result['calls'] = {k: v['calls'] for k, v in sorted(phases.items())}
        # From 3MF read to the end of unwarping: the measured pipeline window.
        result['probe_window'] = round(max(lasts) - min(firsts), 3) if firsts else None
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('job')
    parser.add_argument('label')
    parser.add_argument('--perf', action='store_true')
    parser.add_argument('--repeat', type=int, default=1)
    parser.add_argument('--local-io', action='store_true',
                        help='run on a container-local copy of the job (no virtiofs I/O in the timed run)')
    args = parser.parse_args()
    job = next(p.name for p in sorted(JOBS.iterdir()) if p.name.startswith(args.job))
    results = [run(job, args.label, 'probes', local_io=args.local_io) for _ in range(args.repeat)]
    summary = {'job': job, 'label': args.label, 'repeat': args.repeat, 'local_io': args.local_io,
               'exit_codes': [r['exit_code'] for r in results],
               'gcode_sha256': results[-1]['gcode_sha256'],
               'gcode_hashes_agree': len({r['gcode_sha256'] for r in results}) == 1,
               'container_wall_min': min(r['container_wall'] for r in results),
               'probe_window_min': min((r['probe_window'] for r in results if r.get('probe_window')), default=None),
               'probe_windows': [r.get('probe_window') for r in results],
               # Runs without probe data (rare perf hiccup) are left out, not counted as 0.
               'phases_min': {k: min(r['phases'][k] for r in results if k in r.get('phases', {}))
                              for k in {k for r in results for k in r.get('phases', {})}},
               'calls': results[-1].get('calls')}
    if args.perf:
        # Keep the timing run's G-code; the sampling run writes into its own folder.
        sample = run(job, args.label, 'sample', args.label + '-perf')
        perf_dir = RUNS / (args.label + '-perf')
        summary['perf_report'] = str(perf_dir / job / 'perf-report.txt')
        summary['perf_exit_code'] = sample['exit_code']
    RUNS.mkdir(parents=True, exist_ok=True)
    with (RUNS / 'results.jsonl').open('a') as stream:
        stream.write(json.dumps({**summary, 'at': time.strftime('%Y-%m-%dT%H:%M:%S')}) + '\n')
    print(json.dumps(summary, indent=1))
