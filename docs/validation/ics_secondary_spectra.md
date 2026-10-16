# ICS secondary spectra: formulas, sources and server comparison

Plot review date: 2026-10-17. Data are from the completed single-interaction
server run `em-secondary-20261012T070453-2242518`, physics commit
`cb5b2d9210846c6ba6b46ff2c4bb80ffb47b9751`. Both sampled modes use that
conditional-energy implementation; only the invariant-s sampler is switched
between historical NR and current rejection sampling.

The [fixed-S diagnosis](ics_endpoint_diagnosis.md) identifies both causes of
the residual ICS endpoint error: an invalid Taylor-expansion domain in the
integral and inadequate inverse-CDF endpoint resolution.

The spectra are normalized per **hard ICS interaction**, with secondary photon
energy at least 1e12 eV. Parent electron energies are 1e13, 1e14, 1e15, 1e16,
1e18 and 1e20 eV. The isotropic CMB has T=2.7255 K, z=0, and the frozen photon
range is [0.001,1000] kBT. These are one-interaction distributions, rather than
spectra after cascade propagation.

## Figures

- [Photons: all six parent energies](figures/ics_20261012/ics_photons.png)
  ([vector PDF](figures/ics_20261012/ics_photons.pdf),
  [SVG](figures/ics_20261012/ics_photons.svg)).
- [Recoil electrons: all six parent energies](figures/ics_20261012/ics_recoil_electrons.png)
  ([vector PDF](figures/ics_20261012/ics_recoil_electrons.pdf),
  [SVG](figures/ics_20261012/ics_recoil_electrons.svg)).
- [Photon endpoint at 1e18 and 1e20 eV](figures/ics_20261012/ics_photon_endpoint.png)
  ([vector PDF](figures/ics_20261012/ics_photon_endpoint.pdf),
  [SVG](figures/ics_20261012/ics_photon_endpoint.svg)).
- [Plotted numerical bin data](figures/ics_20261012/ics_spectra.tsv).

Black horizontal segments show theoretical **bin averages**, orange points
historical NR and blue points current sampling. The ordinate is probability
per natural logarithm of secondary energy, dP/dln(E). Each lower panel shows
sampled/theoretical bin probability. The gray bands indicate ratios within
1% and 10%; the original frozen classification determines which target applies
to each region. Error bars use the archived simultaneous confidence intervals;
open downward triangles are upper limits for zero observations. Ratios are
shown where the reference predicts at least 20 events, solely for readability;
this display rule does not alter any test or precision decision.

A logarithmic photon-energy axis compresses the endpoint near E_gamma=E0.
The third figure resolves it using q=1-E_gamma/E0=E_e/E0. It transforms the
**original photon bins**, with the corresponding logarithmic-width Jacobian;
it does not resample or rebin the events.

## Differential law actually used

Use natural units c=hbar=1. Let E0 be the incident electron energy,
s=(p_e+k)^2 the squared invariant collision energy, and define

\[
 a=\frac{m_e^2}{s},\qquad y=\frac{E_e}{E_0},\qquad
 x=\frac{E_\gamma}{E_0}=1-y,\qquad \eta=\frac{E_{\rm cut}}{E_0}.
\]

The soft incident photon energy is neglected in the final energy sum, as in
the production kinematics. The conditional secondary law is the
ultrarelativistic soft-target Klein-Nishina expression, retaining the full
Klein-Nishina dependence rather than restricting to Thomson scattering:

\[
 H(y,a)=\frac1y\left[
 \frac{1+y^2}{2}
 -\frac{2a(y-a)(1-y)}{y(1-a)^2}\right],
 \qquad a\le y\le1-\eta.
\]

This is **equation (12)** in the published-version ELMAG paper,
M. Kachelrieß, S. Ostapchenko and R. Tomàs,
[arXiv:1106.5508v2, section 2.2](https://arxiv.org/pdf/1106.5508v2),
DOI [10.1016/j.cpc.2011.12.025](https://doi.org/10.1016/j.cpc.2011.12.025).
The source prints the distribution proportional to H; the overall factor
cancels when conditioning on s. It is implemented directly in
[secondary_oracle.h](../../tests/sampling/secondary_oracle.h) as
`comptonDensity`, not obtained from the production inverse CDF.

Explicitly, define

\[
 A(s)=\int_a^{1-\eta}H(v,a)\,dv,\qquad
 p_e(y\mid s)=\frac{H(y,a)}{A(s)},\qquad
 p_\gamma(x\mid s)=p_e(1-x\mid s).
\]

There is no extra Jacobian factor for the y-to-x density transformation,
because |dy/dx|=1. A density per physical energy includes the factor 1/E0.
A logarithmic-energy density includes the energy fraction:

\[
 \frac{dP_e}{d\ln E_e}=y\,p_e(y),\qquad
 \frac{dP_\gamma}{d\ln E_\gamma}=x\,p_\gamma(x).
\]

For cross-section normalization, differentiating ELMAG equation (9) gives

\[
 \frac{d\sigma_{\rm ICS}}{dy}
 =\frac{3\sigma_T}{2}\frac{a}{1-a}\,H(y,a),\qquad
 \sigma_{\rm hard}(s)=\frac{3\sigma_T}{2}\frac{a}{1-a}\,A(s).
\]

The oracle uses the normalized form above; the displayed overall factor is
not an additional weight introduced into the sampling comparison.

As a normalization cross-check, a 128-point Gauss integration in log(y)
was compared with equation (9) at s/m_e^2=2,10,1e3,1e8 and
eta=0.1,1e-3,1e-8. The largest relative difference was 7.8e-15.

## Averaging over CMB and the invariant s

The plotted theoretical distribution is not a fixed-s curve. In the
ultrarelativistic limit the invariant weight can be written, up to an
s-independent factor, as

\[
 W_E(s)=(s-m_e^2)\,\sigma_{\rm hard}(s)\,
 I\!\left(\frac{s-m_e^2}{4E_0}\right),\qquad
 I(u)=\int_{\max(u,\epsilon_{\min})}^{\epsilon_{\max}}
 \frac{n_{\rm CMB}(\epsilon)}{\epsilon^2}\,d\epsilon,
\]

\[
 n_{\rm CMB}(\epsilon)=
 \frac{\epsilon^2}{\pi^2(\hbar c)^3}
 \frac1{\exp[\epsilon/(k_BT)]-1}.
\]

This is the invariant-rate representation of **ELMAG equations (6), (8)**.
The full finite-beta kernel used by mcray retains 2E0(1+beta) in the lower
background-integral argument and subtracts the upper-angle term. Both reduce
to the expression above in the regime of this test.

The normalized marginal distributions are

\[
 Z_E=\int W_E(s)\,ds,\qquad
 p_e(y)=\frac1{Z_E}\int W_E(s)\,p_e(y\mid s)\,ds,
 \qquad p_\gamma(x)=p_e(1-x),
\]

with the conditional density set to zero outside its kinematic support and
I(u)=0 for u at or above the maximum background photon energy.
The lower s bound is m_e^2/(1-eta); the upper bound follows from the frozen
maximum background photon energy.

**Shared-input limitation:** the archived oracle obtains W_E(s), including
sigma_hard and the background integral, from the same production kernel as
mcray. The outer quadrature and direct integrations of H are independent of
both tested samplers and the production secondary CDF. Thus the theoretical
curve tests numerical sampling of the specified model; it is not a completely
independent validation of the background convolution or all cross sections.
The saved reference grids use 128/16, 256/24 and 512/32 outer cells/Gauss orders
with kinematic-boundary splitting. Their convergence envelope is included in
the archived report and remains an empirical numerical guard.

## Original physical source

The original unpolarized Compton cross section underlying this differential
law is due to **O. Klein and Y. Nishina**,
*Über die Streuung von Strahlung durch freie Elektronen nach der neuen
relativistischen Quantendynamik von Dirac*, **Zeitschrift für Physik 52,
853–868 (1929)**, DOI
[10.1007/BF01366453](https://link.springer.com/article/10.1007/BF01366453).
That paper is the physical primary source. ELMAG equation (12) is the direct
source for the particular secondary-energy-variable form used by this oracle;
it should not be described as literally equation (12) of the 1929 paper.

## Bin normalization and reproduction

No new simulation was run. For an energy bin [E_i,E_(i+1)], both theory and
sampling are displayed as

\[
 D_i=\frac{P_i}{\ln(E_{i+1}/E_i)},\qquad
 P_i^{\rm MC}=\frac{N_i}{N_{\rm interactions}}.
\]

Consequently the probability integral is the sum of D_i times the logarithmic
bin width, and equals one per secondary species. Using bin integrals avoids
comparing a broad-bin histogram with a theory evaluated only at its midpoint.
For the photon endpoint graph, replace E by q in this bin-width formula; the
ordinate is then dP/dln(q), not dP/dln(E_gamma).

```sh
python tests/sampling/plot_ics_spectra.py \
  --report /path/to/extracted/report.json \
  --reference /path/to/extracted/frozen-reference.json \
  --output-dir /path/to/ics_figures
```

The plotter requires optional Matplotlib and NumPy. All original bins,
counts, reference probabilities and simultaneous intervals remain available
in the saved report. The [server-results review](em_sampling_server_results.md)
contains acceptance decisions and the diagnosis of the shared endpoint defect.
