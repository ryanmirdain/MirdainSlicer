# Native nonplanar implementation continuation

User-selected repository: https://github.com/ryanmirdain/MirdainSlicer
Working branch: `main`; initial checkout: `79985162e4c83a0749dbebeb77dff562aeca45da`.
The user's latest instruction selects this existing Orca version. It supersedes
the pasted request to start a fresh checkout at the release tag.
Upstream `git ls-remote` independently verified on 2026-11-04:
`refs/tags/v2.4.2 = 8500fcdccaa10b5099ac20d252af3a7c560046f1`.
The local checkout has not been reset/rebased, and no branch was created.

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
was built by the fork's own `build_linux.sh`. Initial dependency/full/package
build was started with:

```sh
./build_linux.sh -g -dsi -j 6 -J 2 -r -- -DSLIC3R_CAD=OFF -DBUILD_TESTS=OFF -DORCA_PYTHON_PGO=OFF
```

For an Ubuntu 24.04 host, install build prerequisites using the fork's normal
script, then build dependencies and the application/package without `-g`:

```sh
./build_linux.sh -u
./build_linux.sh -dsi -j 6 -J 2 -r -- -DSLIC3R_CAD=OFF -DBUILD_TESTS=OFF -DORCA_PYTHON_PGO=OFF
```

`-t` must never be supplied. These commands select the host architecture; use
native amd64 and arm64 builders for separate artifacts. Existing dependencies
are in `deps/build/OrcaSlicer_dep/usr/local`, and application output is under
`build`. Resume the same Release build rather than building dependencies again.
The geometry modules were individually compiled with GCC 13 in Ubuntu arm64.
All eight changed desktop/OpenGL ES segment and shadow shader stages compiled
with glslang 15.1, and the four stage pairs linked. This was offline shader
compilation, without opening the GUI or executing a slicing job.

The low-memory application build is being resumed directly in the independent
builder image, keeping the same mounted source/build paths:

```sh
docker run --rm --name mirdain-orca-native-build --user 501:20 \
  --mount type=bind,src="$PWD",dst=/__w/OrcaSlicer/OrcaSlicer \
  -w /__w/OrcaSlicer/OrcaSlicer \
  orcaslicer-linux-builder:ubuntu-24.04--cmake-4.3.0--4fc88e178d94 \
  cmake --build build --config Release --parallel 2
```

Host logs: `/tmp/MirdainSlicer-native-build.log` (dependency/configuration build)
and `/tmp/MirdainSlicer-native-compile.log` (current application compile).
The three new native modules compiled in the full Orca build after fixing enum
qualification with precompiled headers and selecting the public processor API.
Source edits after that compile still need the final incremental build.

## In progress / not yet established

The complete application and package build is still running. GCC's parallel
Release build exhausted the Docker VM's approximately 8 GB memory in an upstream
CGAL translation unit. Resume with lower application build concurrency (the
successful dependency builds can be reused). This document must
be updated with its actual result and artifact paths before delivery. Native
emission/startup template compatibility and preview shader integration still
need source review and compilation. Final bead-envelope/flow details and machine
clearance restrictions must be kept consistent between code and design docs.
No runtime slices, GUI sessions, regression checks or printer jobs have run.
amd64 has not been built. Nothing has been pushed, deployed or submitted upstream.
Do not describe a partial or compile-only implementation as behaviorally verified.
