#!/usr/bin/env python3
"""Topology-only importance guides; never substitute a reference yield."""
from collections import defaultdict, deque
import json
from pathlib import Path
from run_transport_photon_pilots import sha


def products(a, z, code, process):
    if process == 'P':
        result = []
        for count, (da, dz) in zip(map(int, code), [(1, 0), (1, 1), (2, 1), (3, 1), (3, 2), (4, 2)]):
            result.extend([(da, dz)]*count)
            a -= da*count
            z -= dz*count
        result.append((a, z))
    else:
        z += int(code[0])-int(code[1])
        result = []
        for count, (da, dz) in zip(map(int, code[2:]), [(4, 2), (1, 1), (1, 0)]):
            for _ in range(count):
                result.append((da, dz)); a -= da; z -= dz
                if a > 1 and z in [0, a]:
                    result.extend([(1, int(z > 0))]*a)
                    return result
        result.append((a, z))
    if any(a < 1 or z < 0 or z > a for a, z in result):
        raise ValueError('Invalid topology channel')
    return result


def channel_inventory(data):
    result = []
    for process, name in [('P', 'channels.tsv'), ('D', 'decays.tsv')]:
        for line in (data/name).read_text().splitlines():
            if not line or line.startswith('#'):
                continue
            a, z, code = line.split()[:3]
            a, z = int(a), int(z)
            if len(code) != (6 if process == 'P' else 5) or not code.isdigit():
                raise ValueError('Invalid channel code')
            result.append((process, a, z, code, products(a, z, code, process)))
    if len({r[:4] for r in result}) != len(result):
        raise ValueError('Duplicate topology channel')
    return result


def scores(inventory, target):
    if target not in ['minus', 'plus']:
        raise ValueError('Unknown guide target')
    digit = int(target == 'plus')
    reverse, seeds = defaultdict(set), set()
    for process, a, z, code, children in inventory:
        for child in children:
            reverse[child].add((a, z))
        if process == 'D' and int(code[digit]):
            seeds.add((a, z))
    distance = dict.fromkeys(seeds, 0)
    queue = deque(sorted(seeds))
    while queue:
        child = queue.popleft()
        for parent in sorted(reverse[child]):
            if parent not in distance:
                distance[parent] = distance[child]+1
                queue.append(parent)
    rows = []
    for process, a, z, code, children in inventory:
        score = 1. if process == 'D' and int(code[digit]) else max(
            [10.**(-distance[c]) for c in children if c in distance], default=0.)
        rows.append((process, a, z, code, score))
    if (56, 26) not in distance:
        raise ValueError('Target unreachable from Fe56 topology')
    return rows


def prepare(data, output):
    # Caller additionally verifies the frozen bundle manifests and payloads.
    inventory = channel_inventory(data)
    output.mkdir(parents=True, exist_ok=False)
    hashes = {}
    for target in ['minus', 'plus']:
        path = output/(target+'.tsv')
        path.write_text('# Topology guide only: process A Z code score; no physical rates modified\n'+
                        ''.join(f'{p} {a} {z} {c} {s:.17g}\n' for p, a, z, c, s in scores(inventory, target)))
        hashes[path.name] = sha(path)
    manifest = dict(schema='rare-topology-guide-v1', channel_count=len(inventory),
        input_sha256={name: sha(data/name) for name in ['channels.tsv', 'decays.tsv']},
        payload_sha256=hashes, guide_distance_factor=10.,
        note='Shortest channel count to a beta parent; same topology guide for M and C. No reference yield used.')
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return manifest
