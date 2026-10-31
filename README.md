# Temporary CRbeam_v0 review reference

This orphan branch contains a byte-for-byte source snapshot of the historical
CRbeam_v0 fork, published at the user's request for independent comparison in
https://github.com/okolo/mcray/pull/16. It is reference material, not a modern
implementation branch or a proposed merge into main.

Original code, comments and script contents are preserved. Source files,
Python/shell/plotting scripts and CMake/Make definitions retain their relative
paths. Generated files, binaries, scientific data tables and IDE metadata are
excluded. SNAPSHOT_MANIFEST.json lists the SHA-256 of every copied file.
All twelve transport-reference hashes in PR #16 match this snapshot exactly.

Do not modify or run the historical build/utility scripts during review.
Use the modern PR's hash-checking tests/fork_compat/transport.py against an
archive of this exact commit, in a separate scratch directory. Normal modern
builds and CI remain independent of this temporary branch. The PR's AGENTS.md
contains the pinned commit and retrieval/reproduction instructions.
