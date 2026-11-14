# Native Linux EeV field ensemble results

Run `20261113T202246Z-94c347cb` completed on 2026-11-13 on native Ubuntu
22.04 amd64 with GCC 11.4.0 and up to 16 worker processes. It used clean
candidate revision `78570988b942555ff93254a920596137bc2feccb` and reference
snapshot `e1fc71f613efe7fecf3210e15339c34814245cd4`. The launcher exited
with code 0 after 1135.25 seconds (18 min 55 s). The final ensemble and
archive verification took 1007.19 seconds (16 min 47 s).

All 253 CTest cases, 28 Python guards and the 18000-history fixture passed.
The final calculation completed 9000000 histories in each of 18 streams
for **162000000 histories** and zero failed histories. Each stream used its
frozen reservation from the [protocol](../../tests/fork_compat/EEV_FIELD_ENSEMBLE.md).
Those final reservations are now spent.

## Statistical result

All 126 difference intervals contain zero in the simultaneous 95% family
of 144 intervals. The run supports compatibility with zero difference for
the frozen observables at 1 EeV and accuracy 40 in this field/source model.

| Observable | Largest absolute difference | Interval half-width |
| --- | ---: | ---: |
| Hit fraction | 0.00015833 | 0.00098086 |
| Arrival indicator mean | 0.00018800 | 0.00098086 |
| Final coordinate mean | 0.04264 kpc | 0.49043 kpc |
| Emitted direction component mean | 0.00019351 | 0.00196171 |
| Final direction component mean | 0.00019351 | 0.00196171 |

The largest absolute differences in all five groups came from photons in the
rotated cone. The largest ratio of absolute difference to interval half-width
was 0.192, for arrival by 0.098 Mpc of c*time. The 18 individual hit fractions
have a half-width of 0.00069357.

The historical implementation recorded 8776 end-time overruns, all in cone
streams. The largest was 9.997 pc of flight parameter, consistent with the
known observer step floor. The modern implementation recorded zero overruns.
These diagnostics use the same stopping semantics described in the protocol.

## Evidence and local verification

The eighteen compressed streams occupy 15.142 GiB. The server verified each
archive by full gzip decompression and raw SHA-256 comparison. A local
read verified all eighteen compressed SHA-256 values and byte counts against
the server summaries. All stream summaries agreed with the report, the
executable and reference input hashes matched, and all 144 intervals were
recomputed from the retained stream sums with exact agreement to the report.

[Machine-readable evidence](eev_field_ensemble_results.json) retains the
statistical report, configuration, provenance and checksums of the source
artifacts. The compressed histories remain in the run's test-data directory.
