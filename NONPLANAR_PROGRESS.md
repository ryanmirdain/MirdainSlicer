# Native nonplanar implementation continuation

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

Host logs: `/tmp/MirdainSlicer-native-build.log` (dependency/configuration build),
`/tmp/MirdainSlicer-native-heavy-gui.log` (serial large GUI units),
`/tmp/MirdainSlicer-native-app-final.log` (successful full application build), and
`/tmp/MirdainSlicer-native-release.log` (final incremental build).

## In progress / not yet established

The complete Release application, native slicing library, GUI library and Python
module linked successfully. The executable is an aarch64 Linux ELF at
`build/src/Release/orca-slicer`. The latest first-layer metadata changes still
need their incremental compile/relink, followed by AppImage packaging and checksum
recording. Link-time optimization reports an existing TBB `notify_type` enum ODR
warning from dependency sources; this was not changed or behaviorally evaluated.
No runtime slices, GUI sessions, regression checks or printer jobs have run.
amd64 has not been built. Nothing has been pushed, deployed or submitted upstream.
Do not describe a partial or compile-only implementation as behaviorally verified.
