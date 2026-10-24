# Corrected ICS/GammaPP secondary spectra with fast invariant selection

This ports the production corrections from
1cc28f381e961f3d4947d08f0149ac72d63c3635 onto talys based on
e7ffd18e49fdec72246106dc8319dbdc092a9d33. It preserves the previous fast
rejection algorithm for selecting the invariant S. No invariant-CDF sampler,
analytic rejection envelope or additional background-table machinery is added.

## Corrected behavior

ICS now samples the logarithm of the recoil-electron energy fraction. The
photon fraction is recovered with -expm1(log_recoil), avoiding cancellation
near both physical endpoints. Its integrated cross section uses the stable
Compton primitive: the short-interval series is selected by the relative width
(y-x)/x, rather than the absolute width y-x. This cross-section correction can
change rates and the resulting S distribution; the algorithm selecting S is
unchanged.

GammaPP uses a rationalized lower lepton-energy bound and stable primitive
differences. Conditional energy inversion resolves both probability tails with
a tolerance proportional to min(u,1-u), removing the former endpoint cutoff
and pileup near equal splitting. Invalid conditional distributions or failed
inversion raise non-retryable errors, preventing selective history resampling.

These inverse CDFs determine secondary ENERGY at an already selected S; they
are separate from the background invariant sampler. All nontrivial numerical
implementations live in EMSecondarySampling.cpp; its header contains only
short declarations. Formulas and inversion rules match the corrected source
branch, with an out-of-line callable interface replacing the header template.
No new numerical/physical approximation is introduced by this source layout.

## Short verification

The six dedicated GoogleTest cases compare the integrated kernels with
independent positive differential quadrature over threshold, intermediate and
high invariants; exercise actual production secondaries; check recoil and
soft-lepton endpoints, equal splitting, probabilities down to 1e-12 at both
limits, energy balance and fatal invalid-CDF handling. They use the same
assertions and seeds as the corrected source branch. Test helpers live in
one .cpp file; the checks require no historical NR headers or server adapter.

Build with the normal project workflow, then run the focused suite:

```bash
ctest --test-dir build/Release -L sampling --output-on-failure
```

The full short suite also checks interactions with the current talys code.
These checks do not establish 1%/10% precision for every marginal-spectrum bin
or replace the separately recorded server studies.

Background.cpp/.h, MathUtils.cpp/.h and TableFunction.cpp/.h are byte-identical
to the talys base above. In particular, fast S selection retains its existing
sampled maximum and fallback behavior; this port does not provide the analytic
envelope guarantee studied on the experimental branch.

Native macOS arm64 Release verification: all 151 CTest cases passed in
31.03 seconds, including all six new secondary-spectrum cases.

Ubuntu 22.04 Linux amd64 under emulation: all 151 CTest cases passed in
103.70 seconds, including all six secondary-spectrum cases. Both checks used
Release builds; they are short regression checks, not new acceptance ensembles.
Byte-level source invariants and verification counts are retained in
[em_secondary_talys_port_20261024.json](data/em_secondary_talys_port_20261024.json).
