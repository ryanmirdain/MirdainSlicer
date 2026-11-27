"""Runs INSIDE the privileged profiler container (see profile.py). Do not run on the host.

    python3 container_run.py MODE -- orca args...

MODE "probes": registers uprobes (entry and return) on the nonplanar phase functions in
/opt/orca-local/bin/orca-slicer, records them with perf (cheap: a few hundred events), and writes
/job/probes.json with each phase's first entry, last return and summed inclusive time.
MODE "sample": perf record -F 199 --call-graph dwarf over the run, then writes /job/perf-report.txt.
"""

import json
import re
import subprocess
import sys

BINARY = '/opt/orca-local/bin/orca-slicer'
PERF = '/usr/local/bin/perf-any'
TRACING = '/sys/kernel/tracing'
# label -> exact demangled name prefix (up to the argument list).
PHASES = {
    'read_3mf': 'Slic3r::Model::read_from_file(',
    'layerstack_parse': 'Slic3r::Nonplanar::LayerStack::parse(',
    'layerstack_read': 'Slic3r::Nonplanar::LayerStack::read(',
    'print_process': 'Slic3r::Print::process(',
    'object_slice': 'Slic3r::PrintObject::slice()',
    'np_mesh_index': 'Slic3r::Nonplanar::MeshSlicer::MeshSlicer(',
    'np_cross_sections': 'Slic3r::Nonplanar::MeshSlicer::cross_sections(',
    'np_slice_volume': 'Slic3r::Nonplanar::slice_volume(',
    'perimeters': 'Slic3r::PrintObject::make_perimeters()',
    'prepare_infill': 'Slic3r::PrintObject::prepare_infill()',
    'infill': 'Slic3r::PrintObject::infill()',
    'support': 'Slic3r::PrintObject::generate_support_material()',
    'np_clip_extrusions': 'Slic3r::Nonplanar::clip_extrusions(',
    'export_gcode': 'Slic3r::Print::export_gcode(',
    'gcode_do_export': 'Slic3r::GCode::do_export(',
    'np_unwarp': 'Slic3r::Nonplanar::unwarp_gcode(',
    'gcode_process_file': 'Slic3r::GCodeProcessor::process_file(',
}


def stage():
    """Build /opt/orca on the container's own filesystem: uprobes fail (EIO) on virtiofs files and
    perf cannot resolve them by their host path. Runtime entries are symlinked; the binary is copied."""
    import os
    import shutil
    os.makedirs('/opt/orca-local/bin')
    for name in os.listdir('/opt/orca-ro'):
        if name not in ('bin', 'AppRun'):
            os.symlink(f'/opt/orca-ro/{name}', f'/opt/orca-local/{name}')
    shutil.copy2('/opt/orca-ro/AppRun', '/opt/orca-local/AppRun')
    shutil.copy2('/binary/orca-slicer', BINARY)


def symbols():
    listing = subprocess.run(['nm', BINARY], capture_output=True, text=True).stdout.split('\n')
    entries = [line.split() for line in listing if re.match(r'^[0-9a-f]+ [Tt] ', line)]
    demangled = subprocess.run(['c++filt'], input='\n'.join(e[2] for e in entries),
                               capture_output=True, text=True).stdout.split('\n')
    found = {}
    for (address, _, _), name in zip(entries, demangled):
        for label, prefix in PHASES.items():
            if name.startswith(prefix) and '::{lambda' not in name and label not in found:
                found[label] = int(address, 16)
    return found


def file_offset(address):
    text = subprocess.run(['readelf', '-lW', BINARY], capture_output=True, text=True).stdout
    for line in text.splitlines():
        parts = line.split()
        if parts and parts[0] == 'LOAD':
            offset, vaddr, size = int(parts[1], 16), int(parts[2], 16), int(parts[5], 16)
            if vaddr <= address < vaddr + size:
                return address - vaddr + offset
    raise ValueError(hex(address))


def probes(command):
    subprocess.run(['mount', '-t', 'tracefs', 'nodev', TRACING], capture_output=True)
    found = symbols()
    with open(f'{TRACING}/uprobe_events', 'w') as events:
        for label, address in found.items():
            offset = file_offset(address)
            events.write(f'p:orca/{label} {BINARY}:{offset:#x}\n')
            events.write(f'r:orca/{label}__ret {BINARY}:{offset:#x}\n')
    code = subprocess.run([PERF, 'record', '-q', '-e', 'orca:*', '-o', '/job/probes.data', '--', *command]).returncode
    script = subprocess.run([PERF, 'script', '-i', '/job/probes.data', '-F', 'tid,time,event'],
                            capture_output=True, text=True).stdout
    open_calls, result = {}, {}
    for line in script.splitlines():
        match = re.match(r'\s*(\d+)\s+([\d.]+):\s+orca:(\w+?)(__ret)?:', line)
        if not match:
            continue
        tid, stamp, label, ret = int(match[1]), float(match[2]), match[3], match[4]
        entry = result.setdefault(label, {'first': None, 'last': None, 'inclusive': 0., 'calls': 0})
        if ret:
            started = open_calls.get((label, tid), [])
            if started:
                entry['inclusive'] += stamp - started.pop()
            entry['last'] = stamp
        else:
            open_calls.setdefault((label, tid), []).append(stamp)
            entry['calls'] += 1
            entry['first'] = stamp if entry['first'] is None else entry['first']
    json.dump({'exit_code': code, 'phases': result, 'probed': sorted(found)}, open('/job/probes.json', 'w'))
    return code


def sample(command):
    subprocess.run(['sh', '-c', 'echo -1 > /proc/sys/kernel/perf_event_paranoid; echo 0 > /proc/sys/kernel/kptr_restrict'])
    code = subprocess.run([PERF, 'record', '-q', '-F', '199', '--call-graph', 'dwarf,16384',
                           '-o', '/job/perf.data', '--', *command]).returncode
    if code == 0:
        report = []
        for title, arguments in (
                ('SELF (top symbols, % of samples)', ['--no-children', '--percent-limit', '0.3']),
                ('CHILDREN (inclusive)', ['--children', '--percent-limit', '1.5'])):
            out = subprocess.run([PERF, 'report', '-i', '/job/perf.data', *arguments, '--sort', 'symbol',
                                  '--stdio', '-g', 'none', '--no-inline'], capture_output=True, text=True).stdout
            report.append(f'==== {title} ====\n' + '\n'.join(
                line[:220] for line in out.splitlines() if line and not line.startswith('#'))[:60000])
        open('/job/perf-report.txt', 'w').write('\n\n'.join(report))
    return code


def local_io(command):
    """Run on a copy of /job in the container's own filesystem (LOCAL_IO=1): separates code time
    from virtiofs bind-mount I/O. Outputs are copied back to /job afterwards."""
    import os
    import shutil
    shutil.copytree('/job', '/tmp/job')
    os.chdir(os.getcwd().replace('/job', '/tmp/job', 1))
    return [a.replace('/job', '/tmp/job') for a in command]


def copy_back():
    import shutil
    for name in ('plate_1.gcode', 'result.json', '00000.log'):
        try:
            shutil.copyfile(f'/tmp/job/section-01/{name}', f'/job/section-01/{name}')
        except FileNotFoundError:
            pass


if __name__ == '__main__':
    import os
    mode, command = sys.argv[1], sys.argv[sys.argv.index('--') + 1:]
    stage()
    local = os.environ.get('LOCAL_IO') == '1'
    if local:
        command = local_io(command)
    code = probes(command) if mode == 'probes' else sample(command)
    if local:
        copy_back()
    sys.exit(code)
