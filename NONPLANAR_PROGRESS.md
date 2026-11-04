# Native nonplanar implementation and build

User-selected repository: https://github.com/ryanmirdain/MirdainSlicer
Working branch: `main`; initial checkout: `79985162e4c83a0749dbebeb77dff562aeca45da`.
This provided checkout identifies itself as `2.5.0-dev` in `version.inc`.
The user's latest instruction selects this existing Orca version. It supersedes
the pasted request to start a fresh checkout at the release tag.
Upstream `git ls-remote` independently verified on 2026-11-04:
`refs/tags/v2.4.2 = 8500fcdccaa10b5099ac20d252af3a7c560046f1`.
The local checkout has not been reset/rebased, and no branch was created.
Implementation commit: `4f44528708` (native geometry, planning, CLI/GUI,
metadata, preview and dependency build option). Follow-up source review changes
were committed as `231eac0286` (input validation, translations, normal-spacing
planning and preview transitions). Commit `f729520b42` accounts for shared
deposited footprints once, bounds full bead heights for clearance, connects
flat-ended infill to walls, supplies first/last-path template values, rejects
floating first beads, and fixes the CLI's explicit formatting include.
Commit `e90e3b79f3` fixes the generated AppImage helper's source-directory lookup
so existing application metadata is included during packaging. The first complete
application build was started from this clean source commit and linked successfully.
`9c7631b933` connects native first-layer bounds/area/tool metadata to CLI and GUI
project exports; `9bf5a75cb9` accepts either winding of a convex bed outline.

No tests may be written or run. Only source review, compiling and packaging are
authorized. Compilation establishes type/link compatibility; it does not prove
geometric correctness, machine safety, GUI behavior or physical print quality.
The initial upstream dependency build unexpectedly ran CPython's PGO test
workload via `--enable-optimizations` despite `BUILD_TESTS=OFF`. This was caught
during build-script inspection after that dependency finished. No Orca tests
were run or authored. The build was stopped and `ORCA_PYTHON_PGO` was added,
default OFF, to prevent that implicit workload in subsequent builds.

## Implemented source

- `src/libslic3r/Nonplanar/LayerStack.*`: v1 parser, production validation,
  physical bilinear evaluator, gradient, bounds and spacing.
- `src/libslic3r/Nonplanar/Slicer.*`: direct original-mesh intersection,
  complete-band regions, wall/infill/skin planning, 3D segments and native flow.
- `src/libslic3r/Nonplanar/Print.*`, `Print.*`, `PrintApply.cpp`: gated native
  processing/emission, saved sidecar, placement provenance and standard-mode
  isolation.
- CLI registration and integration in `PrintConfig.cpp` and `OrcaSlicer.cpp`.
- GUI import/standard-mode controls in `Plater` and `MainFrame`.
- GCodeProcessor, LibVGCodeWrapper, PathVertex, ViewerImpl and shaders: native
  segment dimensions and band preview plumbing.
- `Model.hpp`: imported layer/provenance metadata follows internal undo snapshots.
- Schema/examples: `resources/nonplanar/`; design/usage: `docs/HLSD/Nonplanar.md`.

## Build environment and commands

Host is macOS arm64; Docker provides native Linux aarch64. No web-slicer image,
source, experiment or service is used. An independent Ubuntu 24.04 build image
was built by the fork's own `build_linux.sh`. Image identity:
`sha256:58484fe69f930399c7f0b75e03d0d75cba6db598947ac5e57c1e234f4790fe3d`.
The dependency/full/package build
command with the CPython workload explicitly disabled is:

```sh
./build_linux.sh -g -dsi -j 6 -J 1 -r -- -DSLIC3R_CAD=OFF -DBUILD_TESTS=OFF -DORCA_PYTHON_PGO=OFF
```

For an Ubuntu 24.04 host, install build prerequisites using the fork's normal
script, then build dependencies and the application/package without `-g`:

```sh
./build_linux.sh -u
./build_linux.sh -dsi -j 6 -J 1 -r -- -DSLIC3R_CAD=OFF -DBUILD_TESTS=OFF -DORCA_PYTHON_PGO=OFF
```

`-t` must never be supplied. These commands select the host architecture; use
native amd64 and arm64 builders for separate artifacts. Existing dependencies
are in `deps/build/OrcaSlicer_dep/usr/local`, and application output is under
`build`. Resume the same Release build rather than building dependencies again.
This build sets `SLIC3R_CAD=OFF`, so the artifact omits OpenCASCADE/STEP support;
native STL/3MF mesh inputs remain compiled. Other platforms and amd64 have no
compiled artifacts from this run. The AppImage follows this fork's host-library
policy: GTK, WebKitGTK 4.1, OpenGL and platform libraries come from the Linux host.
The geometry modules were individually compiled with GCC 13 in Ubuntu arm64.
All eight changed desktop/OpenGL ES segment and shadow shader stages compiled
with glslang 15.1, and the four stage pairs linked. This was offline shader
compilation, without opening the GUI or executing a slicing job.

The initial 8 GB Docker VM needed serial compilation for the largest translation
units and two jobs for the remaining GUI. At the user's explicit authorization,
Docker Desktop was stopped after the successful full link, configured with
`MemoryMiB=32768`, and restarted. All three previously running containers were
restored. Docker now reports 16 CPUs and 33,600,802,816 usable memory bytes.
The final incremental build uses eight jobs in the same independent builder:

```sh
docker run --rm --name mirdain-orca-native-build --user 501:20 \
  --mount type=bind,src="$PWD",dst=/__w/OrcaSlicer/OrcaSlicer \
  -w /__w/OrcaSlicer/OrcaSlicer \
  orcaslicer-linux-builder:ubuntu-24.04--cmake-4.3.0--4fc88e178d94 \
  cmake --build build --config Release --parallel 8
```

For a native Ubuntu 24.04 machine with 32 GB available to the build, `-J 8`
can replace `-J 1` in the full build commands above. Packaging the already built
application separately uses:

```sh
docker run --rm --name mirdain-orca-native-package --user 501:20 \
  --mount type=bind,src="$PWD",dst=/__w/OrcaSlicer/OrcaSlicer \
  -w /__w/OrcaSlicer/OrcaSlicer/build \
  orcaslicer-linux-builder:ubuntu-24.04--cmake-4.3.0--4fc88e178d94 \
  ./src/build_linux_image.sh -i -R Release
```

On the Ubuntu host itself, the equivalent packaging command from `build/` is
`./src/build_linux_image.sh -i -R Release`. The `501:20` container user above is
this checkout's macOS owner/group; select the checkout owner's UID/GID elsewhere.

Host logs: `/tmp/MirdainSlicer-native-build.log` (dependency/configuration build),
`/tmp/MirdainSlicer-native-heavy-gui.log` (serial large GUI units),
`/tmp/MirdainSlicer-native-app-final.log` (successful full application build), and
`/tmp/MirdainSlicer-native-release.log` (successful final incremental build), and
`/tmp/MirdainSlicer-native-package.log` (successful AppImage packaging).

## Completed build artifacts

The complete Release application, native slicing library, GUI library and Python
module linked successfully, including the final first-layer metadata changes.
The final build used clean commit
`93f1a9faba00c19ce61045315268158e0caf4d2b` on 2026-11-04. Subsequent report-only
commits do not change the source compiled into this artifact.

- AppImage: `build/OrcaSlicer_Linux_V2.5.0-dev.AppImage`, aarch64 Linux,
  207,706,632 bytes. SHA256:
  `ed50da37197ea61ca7f3019593d720c32a39d74fbd86d3d52763edb1ef0b4d21`.
- Checksum file: `build/OrcaSlicer_Linux_V2.5.0-dev.AppImage.sha256`.
- Unpacked application: `build/package/`; standalone ELF:
  `build/src/Release/orca-slicer`.
- The final ELF and the packaged `build/package/bin/orca-slicer` have identical
  SHA256 `f0e43509d5629990adeb4ab93990ac54a5df435fbe94c289616ecb623f899e8c`.
- The packaging dependency audit passed. The schema and both example files are
  included under the package's `resources/nonplanar/` directory. AppImage creation
  completed successfully using appimagetool `4919cc4` (build 368).

Link-time optimization reports an existing TBB `notify_type` enum ODR warning
from dependency sources; this was not changed or behaviorally evaluated.
`appstreamcli` was absent, so optional AppStream metadata validation did not run;
the existing application metadata was included in the AppImage.

## Scope and remaining unverified behavior

The basic native source path and Linux arm64 build/package are implemented.
The detailed surface, original-mesh intersection, complete-band containment,
bead/flow, CLI and GUI conventions are in [the design document](docs/HLSD/Nonplanar.md).
The source contains explicit nonplanar-only guards for unsupported features and
configurations. In particular, use a single tool, one watertight volume per object,
classic walls, consistent line widths and rectilinear fills; supports/bridges,
modifiers, multi-tool strategies, ironing, postprocessing and planar caches are
not implemented for native nonplanar jobs. The original standard pipeline remains
selected for ordinary inputs; its runtime regression behavior has not been checked.
No runtime slices, GUI sessions, regression checks or printer jobs have run.
amd64 has not been built. Nothing has been pushed, deployed or submitted upstream.
No functional correctness or physical printing guarantee follows from this build.
Further runtime/GUI/print validation and other architecture builds are follow-up
work. Continue to obey the prohibition on writing or running tests unless the
user changes that instruction. No required build/package stage remains pending.
