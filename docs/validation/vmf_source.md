# Truncated vMF source verification

This is the frozen record of the earlier normalized model, superseded on
2026-10-31 by the user-requested [original weighted contract](vmf_fork_contract.md).
Its point/CDF results must not be presented as verification of that restored
contract. Original cases, criteria and results below remain historical.

The implementation starts from main 5b5399d69e34f6dc3a883ac966bd801edc626f8e
and follows the [historical normalized source contract](https://github.com/okolo/mcray/blob/8cf97128330f93e7f7a8ed25b70a5ce0dbd43056/docs/usage/vmf_source.md).
The new law is opt-in; default pencil and existing cone behavior are preserved.
Gaussian profiles, FFT fields, backgrounds and observer changes are excluded.

## Frozen software controls

The new target unit.vmf.* contains six invariant/point cases and 24
parameterized distribution cases. Three configuration cases check input,
replay and output through the application.

- Compare vmf(kappa=0) and cone directions and the next RNG value exactly
  for openings 0,0.001,10,180 degrees and a rotated axis. Zero opening also
  matches pencil for kappa=0,100 and the maximum finite double.
- Check primary weights, energy, nuclear A/Z/mass, IDs, ancestry pointers,
  time and non-owning custom data. Positions and directions are the intended
  mutable source attributes. Reject invalid/inapplicable constructor inputs.
- Quantile and azimuth point checks use independent Boost.Multiprecision
  cpp_dec_float_100 arithmetic with exp/log and trigonometric functions, not
  the production expm1/log1p path. Nine fixed (kappa,angle) cases cover uniform,
  small/moderate/large kappa*U, a very narrow cone and kappa=1e300. Each uses
  eight seeded pairs (seed 8307). Above product=1000 the reference mass is
  one to its 100-digit precision (exp(-1000)<1e-434). Coordinate tolerances
  are 512*epsilon relative for chord offsets/transverse scale and absolute
  for the longitudinal component; comparisons retain coord_type precision.
- Extreme finite concentration checks exercise the small-product uniform
  limit and the overflow-safe large-product limit, with exactly two draws
  per positive-opening direction. The positive-width underflow rejection is
  checked on platforms where a double input angle can underflow coord_type.
- Fixed ensemble: N=40000 directions per combination of kappa in
  {0,1,100,1000000}, angle in {0.001,10,180} degrees and axis in
  {(0,0,1), normalized(1,-2,3)}. Seed 401715 per case. Validate unit norm,
  support and unchanged weight/energy/ID/position. Recover u from squared
  chord distance to avoid small-angle subtraction cancellation. Transform
  through the analytic CDF; use 20 fixed equal bins with tolerance
  6*sqrt(N*p*(1-p)), p=1/20, and mean tolerance 6*sqrt(1/(12*N)).
  No adaptive reruns, ensemble changes or statistical-tolerance tuning.
- Application checks compare CLI/YAML/config.yaml results byte-for-byte,
  verify required/inapplicable/malformed fields, duplicate keys/options and
  schema v1 exclusion, and compare vmf(kappa=0) with a cone run. A rotated
  massless neutrino ballistic control checks weight, direction metadata and
  X-X0=direction*FlightMpc using named columns and precision-scaled tolerances.
  Zero-angle vmf retains the existing spatial output layout; the complete
  published example is parsed too.

These are short software sampling and integration checks, not astrophysical
acceptance thresholds. They do not reproduce historical-fork trajectories or
establish equivalence of its field-boundary or angular-weight conventions.

An initial emulated Linux point check failed for small kappa*U. A standalone
probe found relative log1p(long double) error about 1.8e-6 at x=-1.52e-14
under the emulator, while expm1 retained precision. The new small-argument
atanh series fixed the same 100-digit point test without changing cases or
tolerances. The original diagnostic and failed run are retained with the logs.

## Results on 2026-10-30

- macOS arm64, Clang 22.1.5 Release: complete suite 230/230 passed;
  the final three interface/example cases also passed separately.
- Emulated Ubuntu 22.04 Linux amd64, GCC/GFortran 11.4.0 Release:
  complete suite 230/230 passed after the small-argument fix.
- The same emulated Linux amd64 source with C/C++/Fortran ASan/UBSan and
  leak detection: complete suite 230/230 passed with no sanitizer reports
  or suppressions. These are local checks, not an assumed hosted CI outcome.

Local artifacts are retained outside Git under build/checks of the feature
worktree. Linux uses the existing image
sha256:37682238eeaf9f01f7eddea17d4e55551d2e317307dfa5c8ea8dd0dc8084baf6,
GCC/GFortran 11.4.0, a 1500-second limit and project-backed temporary paths.
Its sanitizer flags and leak/error settings match the required
[CI jobs](../usage/ci.md). System libraries are not rebuilt with instrumentation.
No TSan or native Linux performance result is asserted.
