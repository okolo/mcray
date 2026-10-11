#!/usr/bin/env python3
"""Audit network PD rates with independent interval quadrature (no acceptance claim).

Photon integrals are analytic for the agreed piecewise density interpolant.
Gauss-Legendre rest-energy quadrature is split at every cross-section node.
Official CRPropa rates and their gamma interpolation remain a separate control.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from numpy.polynomial.legendre import leggauss

MPC_CM = 3.085677581491367e24


def interpolate(x, grid, values):
    j = np.clip(np.searchsorted(grid, x) - 1, 0, len(grid) - 2)
    f = np.log(x / grid[j]) / np.log(grid[j + 1] / grid[j])
    positive = (values[j] > 0) & (values[j + 1] > 0)
    left = np.log(values[j], where=values[j] > 0, out=np.zeros_like(f))
    right = np.log(values[j + 1], where=values[j + 1] > 0, out=np.zeros_like(f))
    result = np.where(positive, np.exp(left * (1-f) + right * f), values[j] * (1-f) + values[j+1] * f)
    return np.where((x >= grid[0]) & (x <= grid[-1]), result, 0.)


class PhotonIntegral:
    """Integral of dn/dE / E^2 from x to the table upper endpoint."""
    def __init__(self, energy, density):
        self.e, self.d = energy, density
        self.width = np.log(energy[1:] / energy[:-1])
        self.positive = (density[:-1] > 0) & (density[1:] > 0)
        self.slope = np.zeros_like(self.width)
        self.slope[self.positive] = np.log(density[1:][self.positive] / density[:-1][self.positive]) / self.width[self.positive]
        self.linear = np.diff(density) / self.width
        self.parts = self.segment(np.arange(len(self.width)), np.zeros_like(self.width), self.width)
        self.tail = np.r_[np.cumsum(self.parts[::-1])[::-1], 0.]

    def segment(self, j, lower, upper):
        # Integrate directly from lower to upper, avoiding subtraction of two
        # nearly equal primitives in steep exponential tails.
        span = upper - lower
        q = self.slope[j] - 1
        factor = np.divide(np.expm1(q * span), q, out=span.copy(), where=np.abs(q) > 1e-12)
        power = self.d[j] / self.e[j] * np.exp(q * lower) * factor
        linear = np.exp(-lower) / self.e[j] * ((self.d[j] + self.linear[j] * lower) * (-np.expm1(-span))
                   + self.linear[j] * (-np.expm1(-span) - span * np.exp(-span)))
        return np.where(self.positive[j], power, linear)

    def __call__(self, x):
        bounded = np.clip(x, self.e[0], self.e[-1])
        j = np.clip(np.searchsorted(self.e, bounded) - 1, 0, len(self.e) - 2)
        lower = np.log(bounded / self.e[j])
        value = self.segment(j, lower, self.width[j]) + self.tail[j+1]
        return np.where(x >= self.e[-1], 0., np.maximum(value, 0.))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--reference-data', type=Path, required=True)
    parser.add_argument('--mcray-rates', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--orders', default='8,16,32')
    args = parser.parse_args()
    orders = [int(v) for v in args.orders.split(',')]
    if len(orders) < 2 or min(orders) < 2 or orders != sorted(set(orders)):
        parser.error('Require at least two increasing quadrature orders >=2')
    # Verify retained sources and official payload before using either.
    manifest = json.loads((args.reference_data/'shared_manifest.json').read_text())
    for name, expected in manifest['payload_sha256'].items():
        if hashlib.sha256((args.reference_data/name).read_bytes()).hexdigest() != expected:
            raise ValueError('Reference payload checksum mismatch: '+name)
    sys.path.insert(0, str(args.reference_data/'sources'))
    import interactionRate
    from units import eV
    energy, density = np.loadtxt(args.reference_data/'shared_photons.tsv').T
    if np.any(~np.isfinite(energy)) or np.any(~np.isfinite(density)) or np.any(energy <= 0) or np.any(np.diff(energy) <= 0) or np.any(density < 0):
        raise ValueError('Invalid photon field')
    photons = PhotonIntegral(energy, density)
    class Field:
        def getDensity(self, x, z=0):
            return interpolate(np.asarray(x)/eV, energy, density)*1e6/eV
    tables = {}
    for filename, channels in [('total.tsv', False), ('channels.tsv', True)]:
        for line in (args.data/filename).read_text().splitlines():
            if not line or line.startswith('#'): continue
            r = line.split(); a,z = map(int,r[:2]); offset = 3 if channels else 2
            tables[a,z,r[2] if channels else 'full'] = np.array(list(map(float,r[offset:])))
    grids = {label: np.loadtxt(args.data/f'grid_{label}.tsv')*1e6 for label in ('talys','external')}
    official_tables = {}
    for line in (args.reference_data/'Photodisintegration/rate_Shared.txt').read_text().splitlines():
        if not line or line.startswith('#'): continue
        r=line.split(); z,n=map(int,r[:2]); official_tables[z+n,z,'full']=np.array(list(map(float,r[2:])))
    for line in (args.reference_data/'exclusive_rates.tsv').read_text().splitlines():
        if not line or line.startswith('#'): continue
        r=line.split(); a,z=map(int,r[:2]); official_tables[a,z,r[2]]=np.array(list(map(float,r[3:])))
    quadratures = {}
    for label, grid in grids.items():
        for order in orders:
            nodes, weights = leggauss(order)
            lo, hi = np.log(grid[:-1]), np.log(grid[1:])
            logx = ((hi+lo)[:,None] + (hi-lo)[:,None]*nodes)/2
            quadratures[label,order] = (np.exp(logx).ravel(), ((hi-lo)[:,None]/2*weights).ravel())
    requested = {}
    failed = []
    for line in args.mcray_rates.read_text().splitlines():
        if line.startswith('# FAILED'): failed.append(line)
        if not line or line.startswith('#'): continue
        a,z,g,code,rate = line.split()
        if code == 'sum': continue
        requested.setdefault((int(a),int(z),code), []).append((float(g),float(rate)))
    rows = []
    for (a,z,code), points in sorted(requested.items()):
        xs = tables.get((a,z,code))
        if xs is None:
            if any(rate != 0 for _,rate in points): raise ValueError('Nonzero missing cross section')
            for gamma,rate in points:
                rows.append(dict(a=a,z=z,channel=code,gamma=gamma,mcray=rate,terminal=True))
            continue
        label = 'external' if a < 12 else 'talys'; grid = grids[label]
        gammas = np.array([g for g,_ in points])
        official = interactionRate.calc_rate_eps(grid*eV, xs*1e-31, gammas, Field())
        values = {}
        for order in orders:
            x,w = quadratures[label,order]
            sigma = interpolate(x,grid,xs)
            values[order] = np.array([np.sum(w*x*x*sigma*1e-27*photons(x/(2*g)))*MPC_CM/(2*g*g) for g in gammas])
        for i,(g,rate) in enumerate(points):
            fine = float(values[orders[-1]][i]); previous = float(values[orders[-2]][i])
            rows.append(dict(a=a,z=z,channel=code,gamma=g,mcray=rate,official_direct=float(official[i]),
                official_gamma_interpolated=float(np.interp(np.log10(g),np.linspace(6,14,201),official_tables[a,z,code])),
                independent=fine, quadrature={str(o):float(v[i]) for o,v in values.items()},
                refinement_absolute=abs(fine-previous),
                mcray_relative=rate/fine-1 if fine else None,
                official_relative=float(official[i])/fine-1 if fine else None))
    summary = {}
    for label, selected in [('full',[r for r in rows if r['channel']=='full']), ('exclusive',[r for r in rows if r['channel']!='full'])]:
        resolved = [r for r in selected if r.get('independent',0) > .1]
        summary[label] = dict(rows=len(selected), resolved_above_0_1_per_mpc=len(resolved))
        for metric in ('mcray_relative','official_relative'):
            summary[label][metric+'_worst'] = max(resolved,key=lambda r:abs(r[metric])) if resolved else None
        summary[label]['max_refinement_relative'] = max((r['refinement_absolute']/r['independent'] for r in resolved),default=None)
    inputs = [Path(__file__), args.mcray_rates, args.data/'total.tsv', args.data/'channels.tsv', args.data/'grid_talys.tsv', args.data/'grid_external.tsv', args.reference_data/'shared_manifest.json']
    result = dict(status='diagnostic only; sparse gamma coverage; no acceptance verdict', summary=summary, failed_cells=failed,
        orders=orders, photon_integral='analytic piecewise; rest quadrature split at cross-section nodes',
        dependencies=dict(python=sys.version,numpy=np.__version__),
        inputs_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}, rows=rows)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    main()
