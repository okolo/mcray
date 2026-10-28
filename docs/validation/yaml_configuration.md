# YAML configuration verification

Verified on 2026-10-20 on the `talys` working tree based on
`7fc04ae6ce2e76887287ba85562b721d6b31ddf8`. This includes uncommitted application
configuration changes and the separately approved missing-PD compatibility
work. These are software regression checks, not final physical acceptance.

The implemented contract and complete input example are maintained in
[YAML configuration](../usage/yaml_configuration.md) and
[config.full.yaml](../../examples/config.full.yaml).

## Short regression suite

| Platform | Toolchain | Result | Full suite wall time |
| --- | --- | --- | --- |
| macOS arm64 | Homebrew Clang 22.1.5, gfortran 15, Release | 143/143 passed | 14.70 s |
| Linux amd64, local Docker emulation | Ubuntu 24.04, GCC 13.3, Release | 143/143 passed | 52.19 s |

Both full suites used `OMP_NUM_THREADS=2` and two CTest workers. The nine
configuration cases were also repeated successfully on both platforms after
the final help-text change. Linux includes yaml-cpp and OpenSSL development
packages, as documented in the build prerequisites.

Coverage includes strict schema/type/enum validation; all thirteen EBL
choices; the complete commented example; relative YAML paths; a valid and a
missing additional background; fixed nuclear installation and checksum
rejection; explicit 64-bit and automatic seeds; effective thread/batch values;
the standard SHA-256 known vector; and failure status recording. A missing
SOPHIA table now produces a diagnostic rather than unchecked column access.

The CLI/YAML/replay comparison uses twelve sampled power-law proton histories,
seed 123456, one worker and batch size twelve. Complete raw particle output
is identical between the legacy command, equivalent YAML input and replay of
the saved configuration on each platform. This establishes equivalence for
that controlled workload; it does not promise identical histories across
platforms or changed execution settings.

## Full-process installation controls

Short Fe56 runs were launched from `/tmp` with a single absolute YAML argument,
without a CLI nuclear-directory override. A runtime directory contained the
normal application tables plus a `tables/nuclear` symlink to the verified
photon-enriched GZK bundle. TALYS, nuclear decays, nuclear GZK, Full nuclear
pair production, both nuclear photon sources and the electromagnetic cascade
were enabled together. Missing-PD and GZK residual policies used the approved
CRPropa compatibility modes.

The workload used two primary histories, total injection energy 5.6e19 eV,
source redshift 1e-8, CMB without EBL, zero magnetic field, seed 2036102000 and
two workers. Both macOS and Linux completed and recorded:

- `tables.nuclear_lookup: standard_tables`;
- 582 reachable missing-PD states as structured provenance;
- actual worker count two and automatic batch size twenty;
- `config.yaml` and `run.yaml`, with no `params` file;
- `status: completed` and zero discarded cascades.

A separate native control requested eight workers with `OMP_THREAD_LIMIT=1`.
It completed with effective worker count one and batch size ten in the saved
configuration. The native full-process controls precede the final help-text
and scalar-string validation edits; Linux repeated the full-process control
with the final application build.

One deterministic CMB coefficient control at the frozen initial energy-grid
node also matched the existing native reference under the unchanged 5e-12
relative identity guard. Its exported table SHA-256 was
`e0ae99527eafa0609598a5b6ee65150a77f5356720cac8eb40e40066f42d0bfc`.
This is a focused regression control, not a repeat of the full numerical matrix.

## Provenance and limits

Run records include executable and application-source SHA-256 fingerprints.
The final Linux full-process executable SHA-256 was
`f37bf8679d7b2bace7a67765ca4068eb5579e3c00a8b6aa992aa787bb32f3c30`;
its application-source fingerprint was
`f0f8956447a106ca2e107e157d66a7ee363407e20bb3116e1ce2ce07539f99c4`.
The Linux read-only source mount did not expose the enclosing submodule Git
directory, so that build honestly records revision `unknown`; the base revision
above and fingerprints identify the tested development source. Published
server runs must resolve and record their actual committed revision.

The two-history installation controls do not constrain spectra, rare products,
physical accuracy or performance of long propagation. The frozen physical
acceptance protocol, tolerances and history budget remain unchanged. No final
acceptance ensemble was launched as part of this configuration work.

## Follow-up: separate settings and remove OpenSSL

The user subsequently requested that simulation settings exist only in
`config.yaml`. The implementation now omits `nuclear.settings`, the legacy
mixed settings/metadata text, process switches and full simulation CLI
arguments from `run.yaml`. Configuration SHA-256, execution information,
the nuclear mass source and reachable missing-PD exceptions remain. Source
manifests are recorded as data provenance under `tables.nuclear_manifests`.

OpenSSL was replaced by a standard-C++ SHA-256 file implementation. No
external hashing command is used; yaml-cpp is the only additional dependency
introduced by the configuration interface. The previous build fingerprints
and OpenSSL prerequisites above describe the original checks, not the updated
binary. Current installation requirements are in the usage contract.

The updated suites pass 144/144 checks on macOS (15.02 s) and Linux amd64
(51.22 s). The additional SHA-256 case covers empty input, padding boundaries,
whole blocks and one
million bytes against fixed reference digests. An independent binary-file
control matches Python `hashlib.sha256` for fourteen lengths from zero through
one million bytes, including 55/56, 63/64/65 and 65535/65536/65537 bytes.
Configuration tests verify that duplicate settings are absent, source manifests
remain available, and `config_sha256` still matches the saved input.

The updated two-history full-process Fe56 controls also complete from `/tmp`
on both platforms using standard nuclear tables, retain all 582 missing-PD
exceptions and three
source-data manifests, and contain no copied settings in `run.yaml`.
After adding source exception-table snapshots, all ten configuration cases
were repeated successfully on both platforms. The snapshots preserve the
terminal whitelist and all ten explicit photon exceptions as source data.
Native dynamic-library inspection confirms no OpenSSL dependency. Its
executable SHA-256 is
`edca0843c881169e7284454503801526269c776387d5355261674963672f45cb`.

## Current SHA-256 backend: OpenSSL Crypto

The user subsequently selected OpenSSL Crypto instead of the private SHA-256
algorithm. `Sha256File` is now a streaming EVP wrapper; the compression rounds,
constants and padding implementation were removed. CMake requires OpenSSL
>=1.1.1, component Crypto. It links only `libcrypto`, uses the API shared by
OpenSSL 1.1.1 and 3.x, and retains all existing SHA-256 manifest fields and
values. Simulation settings remain solely in `config.yaml`.

Verification of this final backend on 2026-10-20:

| Platform | Dependencies | Result | CTest wall time |
| --- | --- | --- | --- |
| macOS arm64, Clang 22.1.5 | Boost 1.88, yaml-cpp 0.8, OpenSSL 3.6.4 | 144/144 passed | 14.80 s |
| Ubuntu 22.04 amd64, Docker emulation, GCC 11 | Boost 1.74.0, yaml-cpp 0.7.0, OpenSSL 3.0.2, CMake 3.22.1 | 144/144 passed | 50.95 s |

The Ubuntu image was built from the public Dockerfile with
`--build-arg UBUNTU_VERSION=22.04`, using only stock distribution development
packages. The image SHA-256 is
`37682238eeaf9f01f7eddea17d4e55551d2e317307dfa5c8ea8dd0dc8084baf6`.
The build script no longer requires `cmake --fresh`: it clears generated
CMake cache/configuration state explicitly, supporting the stock CMake 3.22.
The Docker context includes the commented example required by configuration
tests. Header-only Boost.Hash2 and a newer Boost are not required.

Both suites retain the fixed SHA-256 vectors, padding/block boundaries,
million-byte input, manifest mismatch rejection, configuration hash linkage,
CLI/YAML/replay equivalence and settings/provenance separation checks.
Fourteen binary-file fingerprints also match Python `hashlib.sha256` and the
previous implementation exactly. Dependency inspection on both platforms
confirms `libcrypto` linkage without `libssl` linkage.

A two-history native Fe56 control with all agreed nuclear processes, both
photon sources and EM transport also completes, verifies the standard nuclear
bundle, preserves 582 missing-PD exceptions and writes separate configuration
and provenance records. Its executable SHA-256 is
`a0639008a8a2bafb87418345a5624ffe6b1c8b43e35eccdc58b58bce497e662f`.
These remain short software checks, not the final physical acceptance campaign.
