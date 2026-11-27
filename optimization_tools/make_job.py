"""Create a cached nonplanar Orca job from an uploaded mesh.

Run with the ultimate_slicer venv:
    ../experiments/ultimate_slicer/.venv/bin/python optimization_tools/make_job.py [--mesh PATH | --random] [--seed N]

Runs the experiments/ultimate_slicer unified pipeline (defaults) through step 14 and copies the
step-14 job folder (printer/process/filament.json + section-01/{layers.json, model.3mf, command.json})
to optimization_jobs/<mesh stem>/. Orca outputs from the pipeline's own runtime are not kept.
Prints the job folder on success. Meshes the pipeline rejects or fails on are reported and skipped.
"""

import argparse
import random
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
JOBS = HERE.parent / 'optimization_jobs'
PIPELINE = HERE.parents[1] / 'experiments' / 'ultimate_slicer'
MESHES = HERE.parents[1] / 'experiments' / 'uploaded-3mfs'
sys.path.insert(0, str(PIPELINE))

import batch_run  # noqa: E402
import main  # noqa: E402

KEEP = {'layers.json', 'model.3mf', 'command.json'}


def make(path):
    upload = batch_run.mesh_upload(path)[0]
    step = main.STEP_MODULES['14_slice_sections']
    original, captured = step.run, {}

    def capture(data):
        captured['result'] = yield from original(data)
        return captured['result']
    step.run = capture
    try:
        for event in main.run_pipeline(upload, final_output_only=True):
            if event['type'] == 'error':
                raise RuntimeError(event.get('message', event))
            if event['type'] == 'step_complete' and event['step'] == '14_slice_sections':
                break
    finally:
        step.run = original
    result = captured['result']
    section = result.slice_sections[0]
    if not section.get('gcode_file'):
        raise RuntimeError(f"native slice failed: {section.get('error')}")
    source = Path(result.slice_directory)
    target = JOBS / path.stem
    (target / 'section-01').mkdir(parents=True)
    for name in ('printer.json', 'process.json', 'filament.json'):
        shutil.copyfile(source / name, target / name)
    for name in KEEP:
        shutil.copyfile(source / 'section-01' / name, target / 'section-01' / name)
    (target / 'source.txt').write_text(f'{path}\npipeline slice dir: {source}\n')
    return target


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mesh')
    parser.add_argument('--random', action='store_true')
    parser.add_argument('--seed', type=int)
    args = parser.parse_args()
    if args.mesh:
        print(make(Path(args.mesh)))
        sys.exit()
    files = sorted(p for p in MESHES.iterdir() if p.suffix.lower() in batch_run.FORMATS)
    rng = random.Random(args.seed)
    rng.shuffle(files)
    for path in files:
        if (JOBS / path.stem).exists():
            continue
        try:
            print(make(path))
            break
        except Exception as error:  # report and move on to another mesh
            print(f'skip {path.name}: {error}', file=sys.stderr)
