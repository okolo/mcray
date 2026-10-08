#!/usr/bin/env python3
"""Prepare a labeled, one-channel mcray benchmark bundle; never production input."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools'))
from prepare_nuclear_decay_data import verify
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--A', type=int, required=True)
parser.add_argument('--Z', type=int, required=True)
parser.add_argument('--channel', required=True)
args = parser.parse_args()
verify(args.data)
rows, selected = [], 0
for line in (args.data / 'decays.tsv').read_text().splitlines():
    if not line or line.startswith('#'): continue
    a, z, code, tau = line.split()
    if (int(a), int(z)) == (args.A, args.Z):
        if code != args.channel: continue
        selected += 1
    rows.append(line)
if selected != 1: raise ValueError('Exactly one existing channel must be selected')
shutil.copytree(args.data, args.output)
(args.output / 'decays.tsv').write_text('# CONDITIONAL BENCHMARK ONLY: not unconditional rates\n' + '\n'.join(rows) + '\n')
with (args.output / 'decay_provenance.txt').open('a') as out:
    out.write(f'CONDITIONAL-BENCHMARK-ONLY=A{args.A}Z{args.Z};channel={args.channel};not-production-data\n')
manifest = json.loads((args.output / 'decay_manifest.json').read_text())
manifest['conditional_benchmark_only'] = {'A': args.A, 'Z': args.Z, 'channel': args.channel,
    'origin_manifest_sha256': hashlib.sha256((args.data / 'decay_manifest.json').read_bytes()).hexdigest()}
for name in ('decays.tsv', 'decay_provenance.txt'):
    manifest['payload_sha256'][name] = hashlib.sha256((args.output / name).read_bytes()).hexdigest()
(args.output / 'decay_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
verify(args.output)
