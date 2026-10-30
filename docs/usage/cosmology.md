# Cosmology settings

Both YAML schema versions accept an optional `cosmology` mapping:

```yaml
cosmology:
  hubble_km_s_Mpc: 67.8
  omega_lambda: 0.692
```

These are the existing defaults, now written explicitly in every effective
`config.yaml`, including runs started through the legacy CLI. Old YAML files
without the mapping retain the same defaults. Either field may be omitted;
the other retains its default. No new CLI flags are introduced.

| Field | Meaning and allowed values |
| --- | --- |
| `hubble_km_s_Mpc` | Present-day Hubble constant H0, in km/s/Mpc; finite and strictly positive. |
| `omega_lambda` | Present-day vacuum-energy density parameter; finite, `0 <= omega_lambda < 1`. |

The existing flat matter-plus-Lambda model is preserved:

```text
Omega_m = 1 - Omega_Lambda
H(z) = H0 * sqrt(Omega_m*(1+z)^3 + Omega_Lambda)
```

There is no independent `omega_m` input, curvature/radiation term or selectable
new cosmological model. `omega_lambda: 0` selects the existing matter-only
limit. Pure de Sitter (`omega_lambda: 1`) has no finite universe age in this
implementation and is rejected. Invalid numbers, quoted numeric strings,
unknown/duplicate keys and values outside these ranges fail before output
creation. Internal-unit conversion and the computed universe age must also
remain finite and positive at the platform's time and rate precision;
otherwise initialization fails rather than producing an invalid run. The
Hubble rate must fit the double-precision redshift-loss interface even when
cosmological time uses a wider type.

To select the historical CRbeam fork's application settings, use:

```yaml
cosmology:
  hubble_km_s_Mpc: 70
  omega_lambda: 0.7
```

This selects its cosmology only. CMB model/evolution, EBL, enabled interactions,
source geometry, stopping and thinning remain independent settings. Matching
these two numbers alone does not establish matching cascade spectra.

Cosmology changes redshift-to-time/distance conversion, redshift losses and
the time available for interactions and deflection. The existing integration
accuracy (100000) and initialization range (source maximum redshift plus 10)
are retained. No changes to the cosmological integration formulas or default
particle/RNG behavior are made.

`config.yaml` contains the resolved two input parameters and can be replayed
in a fresh process. `run.yaml` retains its `config_sha256` link without copying
these settings. In an embedded application the library cosmology is global
and initialized once: a conflicting H0/Lambda request, or a source beyond the
initialized table range, fails before output creation with a request to use
a fresh process. It cannot silently retain another run's cosmology.

See [YAML configuration and provenance](yaml_configuration.md) and
[verification](../validation/yaml_cosmology.md).
