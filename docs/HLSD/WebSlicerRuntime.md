# Native runtime releases for the web slicer

This fork keeps only the manually dispatched native-runtime workflow. Inherited
upstream build matrices, tests, scheduled jobs, and repository-maintenance
workflows are removed so pushes and schedules do not consume Actions minutes.

The `Web slicer native runtime` workflow is manually dispatched from the fork's
Actions page. It builds a Release AppImage on a native Ubuntu 24.04 amd64
runner. ARM64 builds remain local and do not consume Actions minutes. It installs
the build script's system dependencies, caches the architecture-specific native dependencies and compiler
objects, and builds without enabling unit tests. Raw profile JSON remains in the
AppImage because the web slicer reads those profiles directly.

Completed native dependencies are saved before compiling the application, so a
source compilation failure does not discard that work. Compiler-cache progress
is saved after successful or failed application compilation. Cold builds still
compile dependencies; subsequent runs reuse the cached dependency installation while its inputs remain unchanged.

The amd64 build must succeed before publishing. A run publishes a prerelease tagged
`web-slicer-<full-commit>-<run-id>-<attempt>`, containing the amd64 AppImage and
`orca-native-release.json`. The tag targets the source commit; run and attempt
IDs prevent later builds from overwriting earlier artifacts. The manifest records
the source commit, native layer format, and each architecture's filename, URL,
and SHA-256. The workflow summary links to the manifest.

The packaging helper assembles the final manifest from the build artifact and
verifies that its bytes match the recorded digest. The web slicer imports that
manifest with
`python3 scripts/pin_native_orca.py <manifest-path-or-https-url>`. Importing the
manifest replaces the runtime pin with the architectures it contains; it does not
rebuild the slicer or retain assets from another source revision. Local ARM64
builds use their own ARM64 pin and matching staged AppImage.

Changes must be committed and pushed before dispatching the workflow. A remote
build cannot include local working-tree edits. Hosted runner access and release
permissions must be enabled on the fork. No separately provisioned AWS builder
is required.
