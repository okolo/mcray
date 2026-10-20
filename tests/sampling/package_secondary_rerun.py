#!/usr/bin/env python3
"""Package the committed correction and portable server tools without a push."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True,help='New package directory')
    a=p.parse_args()
    tools=Path(__file__).resolve().parent
    repo=tools.parents[1]
    def git(*args):
        return subprocess.check_output(['git','-C',str(repo),*args],text=True).strip()
    if git('diff','--name-only','HEAD','--','src/lib','tests/sampling','tests/CMakeLists.txt'):
        raise SystemExit('Commit the production and sampling changes before packaging')
    branch=git('branch','--show-current')
    if not branch:raise SystemExit('A named correction branch is required')
    root=a.output.resolve();root.mkdir(parents=True,exist_ok=False)
    shutil.copytree(tools,root/'sampling',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    git('bundle','create',str(root/'mcray.bundle'),'refs/heads/'+branch)
    manifest={'commit':git('rev-parse','HEAD'),'branch':branch,
        'historical_s_sampler':'d05723f449c1eecee4f00af392014df4c898ebe9',
        'alpha':.01,'core_relative_target':.01,'tail_relative_target':.1,
        'maximum_hours':24,'scope':'Analytic bounded rejection versus frozen rejection and CDF, with corrected conditional spectra; global alpha=1% across both comparisons.',
        'current_algorithm':'analytic_bounded_log_rejection',
        'default_baseline':'both', 'alpha_allocation':{'rejection':.005,'cdf':.005}, 'cdf_control_commit':'fa562d2', 'rejection_commit':'1cc28f381e961f3d4947d08f0149ac72d63c3635'}
    (root/'package-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (root/'README.txt').write_text('Frozen EM secondary-spectrum rerun\n\n'
        'On Ubuntu 22.04 / Python 3.10 (no Python upgrade required):\n'
        '  bash sampling/run_fixed_server.sh --data-dir /home3/ok/mcray-test-data --hours 24 --jobs 4 --raw-gib 20\n\n'
        'The wrapper checks hashes and installs pinned Python requirements in .venv.\n'
        'The 24-hour experiment budget begins when the Python runner starts.\n'
        'System dependencies: git, cmake, make, gcc, g++, gfortran, libgsl-dev,\n'
        'libboost-filesystem-dev, libboost-system-dev, python3-venv, tar, gzip, CA certificates.\n'
        'NR is downloaded at its frozen revision unless --nr /path/to/nr is supplied.\n'
        'For infrastructure verification only: --smoke --hours 0.25 --raw-gib 0.01.\n'
        'Smoke exit code 3 denotes inconclusive and is expected.\n'
        'All original bins, seeds, inference criteria and endpoint checks are retained.\n'
        'Raw events, reference convergence, exact commit, binaries and hashes are archived.\n'
        'Default: bounded rejection against BOTH frozen rejection and CDF controls.\n'
        'Two comparison families receive alpha=0.5% each, preserving global alpha=1%.\nSingle controls: --baseline rejection, --baseline cdf or --baseline nr.\n'
        'See sampling/BOUNDED_REJECTION_SERVER.md for tail quantiles and timing controls.\n')
    entries=[]
    for path in sorted(root.rglob('*')):
        if path.is_file():entries.append(hashlib.sha256(path.read_bytes()).hexdigest()+'  '+str(path.relative_to(root)))
    (root/'SHA256SUMS').write_text('\n'.join(entries)+'\n')
    archive=root.with_suffix('.tar.gz')
    with tarfile.open(archive,'w:gz') as out:out.add(root,arcname=root.name)
    print(archive)
    print('commit='+manifest['commit'])
    print('sha256='+hashlib.sha256(archive.read_bytes()).hexdigest())

if __name__=='__main__':main()
