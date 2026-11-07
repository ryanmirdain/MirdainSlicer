"""Package native Linux runtimes and assemble the web slicer's release manifest."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess


def checksum(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    package = commands.add_parser('package')
    package.add_argument('--arch', choices=('amd64', 'arm64'), required=True)
    package.add_argument('--tag', required=True)
    package.add_argument('--repository', required=True, help='GitHub owner/repository')
    package.add_argument('--output', type=Path, required=True)
    manifest = commands.add_parser('manifest')
    manifest.add_argument('directory', type=Path)
    args = parser.parse_args()

    if args.command == 'package':
        root = Path(__file__).resolve().parents[1]
        revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip()
        version = re.search(r'set\(SoftFever_VERSION "([^"]+)"\)', (root / 'version.inc').read_text())[1]
        source = root / 'build' / f'OrcaSlicer_Linux_V{version}.AppImage'
        filename = f'MirdainSlicer_Linux_Ubuntu2404_{args.arch}_{revision}.AppImage'
        repository = f'https://github.com/{args.repository}'
        release = {
            'repository': repository, 'version': version, 'revision': revision,
            'layer_format': 'orca_nonplanar_layers', 'layer_format_version': 1,
            'assets': {args.arch: {'file': filename, 'sha256': checksum(source),
                                  'url': f'{repository}/releases/download/{args.tag}/{filename}'}},
        }
        args.output.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, args.output / filename)
        (args.output / f'{args.arch}.json').write_text(json.dumps(release, indent=2) + '\n')
        return

    releases = [json.loads((args.directory / f'{arch}.json').read_text()) for arch in ('amd64', 'arm64')]
    release = {key: value for key, value in releases[0].items() if key != 'assets'}
    assets = {}
    for arch, item in zip(('amd64', 'arm64'), releases):
        assert {key: value for key, value in item.items() if key != 'assets'} == release
        asset = item['assets'][arch]
        assert checksum(args.directory / asset['file']) == asset['sha256']
        assets[arch] = asset
    release['assets'] = assets
    (args.directory / 'orca-native-release.json').write_text(json.dumps(release, indent=2) + '\n')


if __name__ == '__main__':
    main()
