"""Time full web slicer slices (layer preparation + Orca + packaging) for a set of meshes.

Runs inside the mirdain_web_slicer backend image, through the same SlicingSession.slice +
run_slice path the slice endpoint uses. A saved feedback report supplies the session template
(printer, presets, nonplanar mode); each mesh replaces the template's object. Appends one JSON
line per slice to --out with the wall time and the backend's own stage metrics.

    docker compose run --rm -T -v <this dir>:/bench:ro [-v <server>:/app/server:ro] \
        [-v <runtime>:/opt/orca:ro] backend python /bench/web_benchmark.py LABEL OUT REPORT MESH...

A mesh argument of `template` slices the report's own object.
"""
import argparse
import io
import json
import os
import shutil
import sys
from pathlib import Path
from time import monotonic, time
from uuid import uuid4

sys.path.insert(0, '/app')
from starlette.datastructures import UploadFile  # noqa: E402

from server.observability.feedback import REPORTS  # noqa: E402
from server.sessions.files import remove_session_files, session_folder  # noqa: E402
from server.sessions.store import TTL, SlicingSession, key, store  # noqa: E402
from server.slicing.engine import import_model_batch, run_slice  # noqa: E402

METRICS = Path(os.environ['SLICER_NONPLANAR_METRICS_DIRECTORY'])


def slice_mesh(report_id, mesh):
    report = json.loads((REPORTS / f'{report_id}.json').read_text())
    state = report['session_state']
    slice_id = next(iter(state['slices']))
    state['slices'][slice_id]['plate_cache'] = {}
    session_id = uuid4().hex
    shutil.copytree(REPORTS / report_id, session_folder(session_id), dirs_exist_ok=True)
    store.set(key(session_id), json.dumps(state), ex=TTL)
    try:
        session = SlicingSession.load(session_id)
        if mesh != 'template':
            for obj in list(session.get_slice(slice_id)['objects']):
                session.delete_object(slice_id, obj['id'], session.state['revision'])
            path = Path(mesh)
            upload = UploadFile(file=io.BytesIO(path.read_bytes()), filename=path.name)
            import_model_batch(session, slice_id, session.state['revision'], [(upload, 1)])
        session = SlicingSession.load(session_id)
        started_at, started = time(), monotonic()
        for _ in session.slice(slice_id, session.state['revision'], run_slice):
            pass
        wall = monotonic() - started
    finally:
        store.delete(key(session_id))
        remove_session_files(session_id)
    runs = [d for d in METRICS.iterdir() if d.stat().st_mtime >= started_at - 1]
    stages = {}
    if runs:
        for line in (max(runs, key=lambda d: d.stat().st_mtime) / 'pipeline.jsonl').read_text().splitlines():
            record = json.loads(line)
            stages[record['stage']] = round(record['seconds'], 3)
    return {'wall': round(wall, 3), 'stages': stages}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('label')
    parser.add_argument('out', type=Path)
    parser.add_argument('report')
    parser.add_argument('meshes', nargs='+')
    args = parser.parse_args()
    for mesh in args.meshes:
        name = 'template' if mesh == 'template' else Path(mesh).stem[:8]
        try:
            row = {'label': args.label, 'mesh': name, **slice_mesh(args.report, mesh)}
        except Exception as error:
            row = {'label': args.label, 'mesh': name, 'error': f'{type(error).__name__}: {error}'[:300]}
        print(json.dumps(row), flush=True)
        with args.out.open('a') as output:
            output.write(json.dumps(row) + '\n')


if __name__ == '__main__':
    main()
