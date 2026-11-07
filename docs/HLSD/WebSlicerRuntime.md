# Native runtime releases for the web slicer

The `Web slicer native runtime` workflow is manually dispatched from the fork's
Actions page. It builds Release AppImages on native Ubuntu 24.04 amd64 and arm64
runners using the same source commit. It installs the build script's system
dependencies, caches the architecture-specific native dependencies and compiler
objects, and builds without enabling unit tests. Raw profile JSON remains in the
AppImage because the web slicer reads those profiles directly.

Both builds must succeed before publishing. A run publishes a prerelease tagged
`web-slicer-<full-commit>-<run-id>-<attempt>`, containing two AppImages and
`orca-native-release.json`. The tag targets the source commit; run and attempt
IDs prevent later builds from overwriting earlier artifacts. The manifest records
the source commit, native layer format, and each architecture's filename, URL,
and SHA-256. The workflow summary links to the manifest.

The packaging helper assembles the final manifest from both build artifacts and
verifies that their source metadata agrees and their bytes match the recorded
digests. The web slicer imports that manifest with
`python3 scripts/pin_native_orca.py <manifest-path-or-https-url>`. Importing the
manifest pins both architectures together; it does not rebuild the slicer.

Changes must be committed and pushed before dispatching the workflow. A remote
build cannot include local working-tree edits. Hosted runner access and release
permissions must be enabled on the fork. No separately provisioned AWS builder
is required.
