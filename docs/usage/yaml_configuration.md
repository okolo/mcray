# YAML configuration and run records

CRbeam accepts one YAML file, or the existing CLI options:

```sh
CRbeam simulation.yaml
CRbeam --particle proton --z 0.01 --emax 1e19 --seed 123 --output results
```

The complete version 1 schema, enum choices, units, inactive-field semantics
and CLI mappings are documented in
[examples/config.full.yaml](../../examples/config.full.yaml). Its explicit
Fe56 settings enable all agreed processes and are not the legacy defaults.
Prepare and install the required data before using that example.

YAML fills the same command model and process setup as the CLI. Omitted
settings retain legacy defaults, except that nuclear tables now have a
standard location. Required YAML source redshift must be positive; nuclear
input requires explicit a>=2 and z>0. All existing cross-process restrictions
remain, including Full nuclear pair production requiring EM transport and
CRPropa GZK residuals requiring nuclear decays. The physical kernels and their
approved approximations are unchanged.

## Input rules

- One positional `.yaml` or `.yml` file; YAML and simulation CLI arguments
  cannot be mixed. A CLI option value ending in `.yaml` remains an option value.
- `schema_version: 1` or `2` is mandatory. Version 1 retains its original
  magnetic-field behavior. Version 2 adds explicitly selected grid/none/sum
  fields; see [the magnetic-field contract](grid_magnetic_fields.md).
  It also supports the optional [photon-thinning policy](photon_thinning.md)
  and the [truncated vMF source](vmf_source.md).
  Unknown/duplicate keys, invalid enum values,
  unsupported types, nonfinite numbers and invalid ranges are errors.
- Booleans are `true`/`false`. Numeric values must be numbers, not quoted
  numeric strings. Units are part of field names; unit expressions are unsupported.
- Use exactly one monoenergetic or power-law spectrum. Monoenergetic energy
  and power-law endpoints are total particle energies, including for nuclei.
- Relative YAML paths resolve against the YAML file's directory. Existing
  CLI path semantics are preserved. Absolute additional-background paths are
  supported by the shared matrix-background loader.
- Includes, anchors, aliases, merge keys, explicit YAML tags and multiple
  documents are unsupported. There is no environment-variable expansion.
- `simulation.seed` accepts a nonnegative integer fitting unsigned long on
  the platform, or `auto`. The Linux amd64/macOS targets support 64-bit seeds.
  CLI explicit seeds retain their existing unsigned-32-bit input contract.
- `threads` accepts a positive integer or `auto`; it uses the same OpenMP
  execution mechanism as the CLI. Runtime limits determine the effective
  team size, which is then fixed and recorded. `batch_size: auto` resolves to
  ten times that size; trajectory logging instead uses the full primary count.

A seed alone does not promise identical trajectories for different
thread/batch settings, platforms or data. These controls are preserved in
both the effective configuration and run provenance.

## Cosmology

Both schema versions support [cosmology settings](cosmology.md):
`cosmology.hubble_km_s_Mpc` and `cosmology.omega_lambda`. Defaults remain
67.8 km/s/Mpc and 0.692; effective snapshots record them explicitly. The
matter density is derived from flatness. Existing CLI invocations retain the
same cosmology and acquire these explicit values in their saved YAML.

## Nuclear tables

The default nuclear package is **`<standard tables root>/nuclear`**. The root
uses the existing application mechanism: `../share/<executable>/tables`, then
`tables` beside the executable, then the existing working-directory `tables`
fallback. This is also the root for the other application tables. Normal
installed deployments therefore do not depend on the YAML file's location.
There is no `nuclear.data_directory` or separate PD/photon path in YAML.

Prepare a verified immutable bundle using the existing public preparation
tools. For full Fe56 support it must contain the approved GZK mass/decay
compatibility data, raw PD cross sections and corrected decay/photon tables,
including the ten explicit photon exceptions. Install the complete bundle,
not selected individual TSV files. For example, from a shell with appropriate
access to the selected installation:

```sh
# TABLES_ROOT is the existing application's table directory; choose its real path.
# The destination must not already contain another bundle.
cp -a /path/to/verified/photon-enriched-gzk-bundle "$TABLES_ROOT/nuclear"
```

A directory symlink to the complete immutable bundle is also supported.
Do not replace the other table directories, such as `sophia2`, when installing
nuclear data. Standard nuclear packages require a nonempty payload hash map
in `manifest.json`; available nuclear/decay/photon manifest payloads are
verified before loading. Missing required files and checksum mismatches fail
before creating simulation output.

The old `--talys-data` override remains available **only in CLI mode** for
compatibility with existing commands. Its actual directory and hashes are
recorded as `legacy_cli_override` in `run.yaml`. It has no YAML counterpart.
To replay an exported configuration from such a command, install the matching
recorded bundle in the standard location. Replaying any configuration also
requires matching the recorded external table versions; the YAML does not
embed numerical datasets.

## Output records

Both input forms now write **`config.yaml` and `run.yaml` instead of `params`**.
This deliberately changes the result-metadata format; consumers of `params`
must migrate. Particle output files and their schemas are unchanged.

`config.yaml` is a valid simulation input containing resolved settings,
including the actual isotope charge, random seed, worker/batch counts and
absolute selected paths. It preserves legacy clamps such as raising the
power-law injection minimum to the propagation cutoff. Set a new output
directory before replaying it; overwrite is false by default. It contains
no run-only fields or configurable nuclear table path.

`run.yaml` records:

- input mode, executable path, YAML path and SHA-256 when applicable;
- executable hash, build revision, tracked patch hash, compiled application
  source fingerprint, compiler/build type and relevant runtime environment;
- standard table root, installed-file SHA-256 inventory, additional background
  hash and active nuclear package directory/hashes;
- explicit/automatic seed origin, nuclear mass source, source-data manifests
  and structured reachable missing-PD exceptions;
- UTC start/finish timestamps, execution status and discarded-cascade count.

Simulation settings belong only in `config.yaml`. `run.yaml` links to that
file through `config_sha256`; it does not copy nuclear settings, process
switches, legacy mixed settings/metadata text or simulation CLI arguments.
Source-data manifests and explicit terminal/photon exception tables describe
the installed data rather than simulation settings and remain part of data
provenance.

The installed-file inventory is labelled as an inventory, not a claim every
file was used. Directory symlinks are followed with cycle protection; the
current output directory is excluded from its own data inventory. The
executable hash identifies the exact binary.

Records are initially written before background/process initialization, then
updated on completion or failure. A failure before output creation, such as
invalid input or failed nuclear preflight, creates no result directory.
`completed`, `failed` and `interrupted` describe execution, not physical
acceptance. A hard kill can leave `running`; that does not establish that a
process is still alive or that its partial output is valid. Writes use a
same-directory temporary file and rename.

## Build dependencies and checks

The application additionally requires yaml-cpp >=0.6 and OpenSSL Crypto
>=1.1.1. File hashing uses the EVP SHA-256 interface, compatible with OpenSSL
1.1.1 and 3.x, with no private algorithm implementation or external command.
Only `libcrypto` is linked; OpenSSL's TLS library `libssl` is not needed.
The physics library and standalone numerical adapters do not depend on YAML.
For the full CRbeam build, install the development packages alongside the
existing compiler/Fortran/Boost/GSL dependencies:

```sh
# Debian/Ubuntu, including the existing build dependencies
sudo apt-get install build-essential cmake gfortran libgsl-dev \
  libboost-filesystem-dev libboost-system-dev libyaml-cpp-dev libssl-dev

# macOS/Homebrew
brew install yaml-cpp openssl@3
```

Ubuntu 22.04's stock Boost 1.74, yaml-cpp 0.7, OpenSSL 3 and CMake 3.22
packages suffice; Boost.Hash2 and a Boost upgrade are not required. Build with
`bash bin/build.sh Release` from the repository root. The build script supports
CMake versions without `--fresh` and locates Homebrew OpenSSL on macOS.
For direct CMake configuration with a nonstandard OpenSSL installation, pass
`-DOPENSSL_ROOT_DIR=/path/to/openssl`.

Run the short regression suite after building:

```sh
OMP_NUM_THREADS=2 ctest --test-dir build/Release --output-on-failure --parallel 2
```

The Docker image includes these dependencies. GoogleTest checks strict YAML
parsing, all EBL choices, the complete example, seeded CLI/YAML/snapshot
history equivalence, path resolution, fixed nuclear installation, checksum
rejection, actual seed/batch recording, SHA-256 and failed-run metadata.
Long physical acceptance remains a separate server workload.

Source geometry, strict sphere membership and the compatible spatial table
format are described in [the spatial input contract](spatial_source_observer.md).
