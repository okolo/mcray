#!/usr/bin/env python3
"""Bounded external Linux amd64 experiment; never start a long local run."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import sys
import time
import traceback
from datetime import datetime, timezone

OLD='d05723f449c1eecee4f00af392014df4c898ebe9'
REJECTION='1cc28f381e961f3d4947d08f0149ac72d63c3635'
NR_HEADERS=('nr3.h','stepper.h','odeint.h','stepperdopr5.h')
EDGES=[0,1e-12,1e-10,1e-8,1e-6,1e-5,1e-4,1e-3,.01,.03,.1,.2,.35,.5,.65,.8,.9,.97,.99,.999,.9999,1-1e-5,1-1e-6,1-1e-8,1-1e-10,1-1e-12,1]
CASES=[('ICS',e) for e in (1e13,1e14,1e15,1e16,1e18,1e20)]+[('PP',e) for e in (1e14,1e15,1e16,1e18,1e20)]

def save(path,value):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    temporary.replace(path)

def seed(label):
    return 1+int.from_bytes(hashlib.sha256(('mcray-secondary-v1:'+label).encode()).digest()[:4],'big')%4294967294

def run(cmd,log,deadline,cwd=None):
    remaining=deadline-time.monotonic()
    if remaining<=0:
        raise TimeoutError('Experiment deadline reached')
    with log.open('ab') as out:
        out.write(('\nCOMMAND '+json.dumps(list(map(str,cmd)))+'\n').encode());out.flush()
        p=subprocess.Popen(list(map(str,cmd)),stdout=out,stderr=subprocess.STDOUT,cwd=cwd,start_new_session=True)
        try:
            rc=p.wait(timeout=remaining)
        except BaseException as error:
            try:os.killpg(p.pid,signal.SIGTERM)
            except ProcessLookupError:pass
            try:p.wait(timeout=min(2,max(.1,deadline-time.monotonic())))
            except subprocess.TimeoutExpired:pass
            # The parent may exit while a descendant ignores TERM.
            try:os.killpg(p.pid,signal.SIGKILL)
            except ProcessLookupError:pass
            p.wait()
            if isinstance(error,(subprocess.TimeoutExpired,TimeoutError)):
                raise TimeoutError('Deadline terminated process group: '+str(cmd[0])) from error
            raise
    if rc:
        raise RuntimeError('Command failed ('+str(rc)+'): '+str(cmd[0])+'; see '+str(log))

NR_REPOSITORY='https://github.com/blackstonep/Numerical-Recipes.git'

def prepare_nr(local, root, tools, deadline, repository=NR_REPOSITORY):
    """Freeze and validate historical headers before compiling either sampler."""
    expected=json.loads((tools/'nr_reference_hashes.json').read_text())
    nr=root/'nr-private';nr.mkdir()
    provenance={'historical_commit':expected['commit'],'source':'local' if local else repository,
                'expected_headers':expected['headers']}
    if local is not None:provenance['local_directory']=str(local.resolve())
    save(root/'nr-provenance.json',provenance)
    if local is None:
        upstream=root/'nr-upstream'
        run(['git','clone','--no-checkout',repository,upstream],root/'logs/nr-fetch.log',deadline)
        run(['git','-C',upstream,'checkout','--detach',expected['commit']],root/'logs/nr-fetch.log',deadline)
        actual=subprocess.check_output(['git','-C',str(upstream),'rev-parse','HEAD'],text=True).strip()
        if actual!=expected['commit']:raise RuntimeError('Downloaded NR revision mismatch')
        provenance['resolved_commit']=actual
        local=upstream
    else:
        provenance['local_directory']=str(local.resolve())
    hashes={}
    for name in NR_HEADERS:
        shutil.copy2(local/name,nr/name)
        hashes[name]=hashlib.sha256((nr/name).read_bytes()).hexdigest()
    provenance['headers']=hashes
    save(root/'nr-provenance.json',provenance)
    for name in NR_HEADERS:
        if hashes[name]!=expected['headers'][name]:
            raise RuntimeError('Historical NR header hash mismatch: '+name)
    return nr

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--server',action='store_true',required=True,help='Explicit external server acknowledgement')
    p.add_argument('--repo',required=True,help='Remote URL or server-side Git repository')
    p.add_argument('--ref',default='talys')
    p.add_argument('--baseline',choices=('rejection','nr'),default='rejection',
                   help='Frozen rejection control (1cc28f3), or historical NR control')
    p.add_argument('--nr',type=Path,help='Local historical NR headers; otherwise download the pinned NR repository')
    p.add_argument('--data-dir',type=Path,required=True)
    p.add_argument('--hours',type=float,default=24)
    p.add_argument('--jobs',type=int,default=4,help='Build threads; sampling is sequential')
    p.add_argument('--raw-gib',type=float,default=20,help='Maximum planned uncompressed event bytes')
    p.add_argument('--smoke',action='store_true',help='Two cases, small grids and 1000 events; no acceptance claim')
    a=p.parse_args()
    if platform.system()!='Linux' or platform.machine() not in ('x86_64','amd64'):
        p.error('Run on an external native Linux amd64 server')
    if not 0<a.hours<=24 or a.jobs<1 or not 0<a.raw_gib<=1000:
        p.error('Require 0 < hours <= 24, positive jobs and raw-gib')
    started=time.monotonic()
    shutdown_guard=min(5,a.hours*36)
    deadline=started+a.hours*3600-shutdown_guard
    def alarm_handler(signum, frame):
        raise TimeoutError('Total experiment deadline reached')
    signal.signal(signal.SIGALRM,alarm_handler)
    signal.setitimer(signal.ITIMER_REAL,a.hours*3600-shutdown_guard)
    reserve=min(3600,a.hours*360);work_deadline=deadline-reserve
    a.data_dir.mkdir(parents=True,exist_ok=True)
    root=a.data_dir.resolve()/('em-secondary-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S')+'-'+str(os.getpid()))
    root.mkdir();logs=root/'logs';logs.mkdir();results=root/'results';results.mkdir()
    status={'status':'running','started_utc':datetime.now(timezone.utc).isoformat(),'arguments':{k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()},'maximum_seconds':a.hours*3600,'work_seconds':work_deadline-started}
    save(root/'status.json',status)
    exitcode=1;source=root/'source';tools=Path(__file__).resolve().parent
    try:
        # Preflight dependencies before fetching/building.
        for name in ('git','cmake','ctest','gcc','g++','gfortran','tar','gzip'):
            if not shutil.which(name):raise RuntimeError('Missing tool: '+name)
        import scipy, numpy
        import secondary_statistics as stats
        nr=prepare_nr(a.nr,root,tools,work_deadline)
        run(['git','clone','--no-checkout',a.repo,source],logs/'fetch.log',work_deadline)
        run(['git','-C',source,'fetch','origin',a.ref],logs/'fetch.log',work_deadline)
        run(['git','-C',source,'checkout','--detach','FETCH_HEAD'],logs/'fetch.log',work_deadline)
        commit=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()
        # Resolve historical code from this same repository, never a moving old branch.
        run(['git','-C',source,'cat-file','-e',OLD+':src/lib/MathUtils.cpp'],logs/'fetch.log',work_deadline)
        run(['git','-C',source,'archive','--format=tar','--output='+str(root/'historical-source.tar'),OLD],logs/'fetch.log',work_deadline)
        status['resolved_commit']=commit;status['historical_commit']=OLD;save(root/'status.json',status)
        status['current_algorithm']='full_support_numerical_log_cdf'
        status['rejection_control_commit']=REJECTION
        status['baseline_commit']=REJECTION if a.baseline=='rejection' else OLD
        labels={'old':'Frozen rejection (1cc28f3)' if a.baseline=='rejection' else 'Historical NR (d05723f)',
                'current':'Numerical log-CDF inversion','old_minus_current':'Baseline minus CDF'}
        modes={'old':'rejection' if a.baseline=='rejection' else 'old','current':'current'}
        save(root/'sampler-labels.json',labels)
        run(['git','-C',source,'archive','--format=tar','--output='+str(root/'frozen-source.tar'),commit],logs/'fetch.log',work_deadline)
        run(['git','-C',source,'update-ref','refs/reproducibility/current',commit],logs/'fetch.log',work_deadline)
        run(['git','-C',source,'update-ref','refs/reproducibility/historical',OLD],logs/'fetch.log',work_deadline)
        run(['git','-C',source,'bundle','create',root/'frozen-history.bundle','refs/reproducibility/current','refs/reproducibility/historical'],logs/'fetch.log',work_deadline)
        overlay=root/'adapter';shutil.copytree(tools,overlay,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        shutil.copytree(overlay,source/'tests/sampling',dirs_exist_ok=True)
        cmake=source/'tests/CMakeLists.txt'
        with cmake.open('a') as out:out.write('\ninclude("${CMAKE_CURRENT_SOURCE_DIR}/sampling/secondary_targets.cmake")\n')
        build=root/'build'
        run(['cmake','-S',source/'src/app/crbeam','-B',build,'-DCMAKE_BUILD_TYPE=Release','-DCMAKE_POLICY_VERSION_MINIMUM=3.5','-DMCRAY_NR_REFERENCE_DIR='+str(nr)],logs/'configure.log',work_deadline)
        run(['cmake','--build',build,'--parallel',a.jobs],logs/'build.log',work_deadline)
        run(['ctest','--test-dir',build,'--output-on-failure'],logs/'ctest.log',work_deadline)
        reproducibility=root/'reproducibility';reproducibility.mkdir()
        shutil.copy2(build/'tests/em_secondary_probe',reproducibility/'em_secondary_probe')
        shutil.copy2(build/'CMakeCache.txt',reproducibility/'CMakeCache.txt')
        shutil.copytree(build/'tests/sampling_reference',reproducibility/'generated')
        run([sys.executable,'-m','unittest','discover','-s',overlay,'-p','test_secondary*.py'],logs/'python-tests.log',work_deadline)
        run(['bash','-c','uname -a; gcc --version; g++ --version; gfortran --version; cmake --version; ldd "$1"; df -h "$2"','versions',build/'tests/em_secondary_probe',root],logs/'versions.log',work_deadline)
        run([sys.executable,'-m','pip','freeze','--all'],logs/'pip-freeze.log',work_deadline)
        if shutil.which('dpkg-query'):
            run(['dpkg-query','-W'],logs/'system-packages.log',work_deadline)
        elif shutil.which('rpm'):
            run(['rpm','-qa'],logs/'system-packages.log',work_deadline)
        if Path('/etc/os-release').is_file():
            shutil.copy2('/etc/os-release',root/'os-release')
        save(root/'python-environment.json',{'python':sys.version,'numpy':numpy.__version__,'scipy':scipy.__version__})
        binary=build/'tests/em_secondary_probe';edges=root/'edges.txt';edges.write_text('\n'.join(format(e,'.17g') for e in EDGES)+'\n')
        cases=CASES if not a.smoke else [('ICS',1e15),('PP',1e15)]
        # Independent quantile checks and three repeated full-cost timing
        # blocks per sampler. These use separate processes/seeds and do not
        # change the original spectrum family, bins, pilots or main seeds.
        import invariant_diagnostics
        diagnostics=root/'invariant-diagnostics';diagnostics.mkdir()
        diagnostic_paths=[]
        for process,energy in cases:
            key=process+'-'+format(energy,'.0e');path=diagnostics/(key+'.tsv')
            run([binary,'diagnostic',process,format(energy,'.17g'),'1e12','0.001',edges,path],
                logs/(key+'-invariant.log'),work_deadline,cwd=source/'bin')
            diagnostic_paths.append(path)
        controls=invariant_diagnostics.summarize(diagnostic_paths)
        save(root/'invariant-diagnostics.json',controls)
        if controls['status']!='within_control':
            raise RuntimeError('Invariant quantile controls unresolved; see invariant-diagnostics.json')
        references={};pilots={};case_args={};counts={};completed={};plan={}
        for process,energy in cases:
            key=process+'-'+format(energy,'.0e');directory=results/key;directory.mkdir()
            args=[process,format(energy,'.17g'),'1e12','0.001',edges];case_args[key]=args
            paths=[]
            for cells,order in ([(8,8),(16,12),(32,16)] if a.smoke else [(128,16),(256,24),(512,32)]):
                path=directory/('oracle-'+str(cells)+'.tsv');paths.append(path)
                run([binary,'oracle',*args,cells,order,path],logs/(key+'-oracle.log'),work_deadline,cwd=source/'bin')
            references[key]=stats.freeze_reference(paths)
            pilots[key]={}
            for mode in ('old','current'):
                begin=time.monotonic();n=100 if a.smoke else 1000
                run([binary,'sample',*args,modes[mode],n,seed(key+mode+'pilot'),directory/(mode+'-pilot.events'),directory/(mode+'-pilot.tsv')],logs/(key+'-pilot.log'),work_deadline,cwd=source/'bin')
                pilots[key][mode]=(time.monotonic()-begin)/n
        save(root/'frozen-reference.json',references)
        save(root/'pilot-timing.json',pilots)
        available=work_deadline-time.monotonic()
        if available<=0:raise TimeoutError('No time remains after references/pilots')
        # Allocate equal wall time per case using BOTH measured sampler costs.
        # Samples, seeds and chunk boundaries are frozen before the first main event.
        disk_limit=int(a.raw_gib*1024**3/(2*len(cases)*26))
        if shutil.disk_usage(root).free < a.raw_gib*1024**3*2+2*1024**3:
            raise RuntimeError('Insufficient disk for raw data, archive and build; lower --raw-gib')
        for key in references:
            n=1000 if a.smoke else min(disk_limit,int(.65*available/len(cases)/sum(pilots[key].values())))
            if n<1000:raise RuntimeError('Insufficient runtime or disk budget for the frozen experiment')
            chunks=[]
            for i,offset in enumerate(range(0,n,100000)):
                chunks.append({'index':i,'N':min(100000,n-offset),'seeds':{mode:seed(key+mode+str(i)) for mode in ('old','current')}})
            plan[key]={'N':n,'chunks':chunks};counts[key]={'old':[],'current':[]};completed[key]={'old':0,'current':0}
        all_seeds=[c['seeds'][m] for x in plan.values() for c in x['chunks'] for m in ('old','current')]
        all_seeds += [seed(k+m+'pilot') for k in plan for m in ('old','current')]
        if len(all_seeds)!=len(set(all_seeds)):raise RuntimeError('Seed collision; frozen design invalid')
        from scipy.stats import norm
        family=sum(len(ref)*9 for ref in references.values())
        z=float(norm.isf(.01/(2*family)))
        forecast={k:{oid:{'expected_multiplicity_counts':[plan[k]['N']*p for p in r['p']],
                         'approximate_simultaneous_relative_halfwidth':None if r['mean']==0 else z*(max(0.,r['p'][1]+4*r['p'][2]-r['mean']**2)/plan[k]['N'])**.5/r['mean'],
                         'target':r['target']} for oid,r in ref.items()} for k,ref in references.items()}
        save(root/'expected-statistics.json',{'diagnostic_only':'normal approximation for runtime planning; exact intervals determine results','cases':forecast})
        save(root/'frozen-plan.json',{'alpha':.01,'core_relative_target':.01,'tail_relative_target':.1,'cases':plan,'background':'CMB T=2.7255 K, z=0, k/T=[0.001,1000], grid=0.01','Emin_eV':1e12,'sampler_tolerance':.001,'smoke':a.smoke,'sampler_labels':labels,'baseline_commit':status['baseline_commit']})
        for key in plan:
            for c in plan[key]['chunks']:
                for mode in ('old','current'):
                    stem=results/key/(mode+'-'+str(c['index']))
                    run([binary,'sample',*case_args[key],modes[mode],c['N'],c['seeds'][mode],str(stem)+'.events',str(stem)+'.tsv'],logs/(key+'-'+mode+'.log'),work_deadline,cwd=source/'bin')
                    if Path(str(stem)+'.events').stat().st_size != 8+26*c['N']:
                        raise RuntimeError('Raw event count/record size mismatch')
                    counts[key][mode].append(Path(str(stem)+'.tsv'));completed[key][mode]+=c['N']
                    save(root/'progress.json',completed)
        aggregate={k:{m:stats.read_counts(paths) for m,paths in modes.items()} for k,modes in counts.items()}
        save(root/'diagnostics.json',{k:{m:[path.read_text().splitlines()[0] for path in paths] for m,paths in modes.items()} for k,modes in counts.items()})
        # Analysis in a subprocess is also bounded by the total deadline.
        save(root/'analysis-input.json',{'references':references,'counts':aggregate,'planned_n':{k:v['N'] for k,v in plan.items()},'complete':not a.smoke,'sampler_labels':labels})
        run([sys.executable,overlay/'analyze_secondary.py',root/'analysis-input.json',root/'report.json'],logs/'analysis.log',deadline-reserve*.3)
        report=json.loads((root/'report.json').read_text());status['status']=report['status'];status['smoke']=a.smoke
        exitcode={'within_target':0,'deviation':2,'inconclusive':3}[report['status']]
        status['completed_utc']=datetime.now(timezone.utc).isoformat()
    except TimeoutError as e:
        status.update(status='deadline',error=str(e));exitcode=124
    except Exception as e:
        status.update(status='execution_error',error=str(e),traceback=traceback.format_exc());exitcode=1
    finally:
        status['elapsed_before_archive_seconds']=time.monotonic()-started;status['exit_code']=exitcode
        save(root/'status.json',status)
        # Retain the checkout, binaries, raw events, partial chunks and logs even
        # when analysis or compression fails. Include licensed headers privately.
        try:
            manifest={}
            for path in root.rglob('*'):
                if not path.is_file() or '.git' in path.parts or 'build' in path.relative_to(root).parts:continue
                h=hashlib.sha256()
                with path.open('rb') as f:
                    for block in iter(lambda:f.read(1024*1024),b''):
                        h.update(block)
                        if time.monotonic()>=deadline:raise TimeoutError('Deadline during checksums')
                manifest[str(path.relative_to(root))]=h.hexdigest()
                if time.monotonic()>=deadline:raise TimeoutError('Deadline during checksums')
            save(root/'sha256.json',manifest)
            archive=root.with_suffix('.tar.gz')
            run(['tar','--exclude=./source/.git','--exclude=./nr-upstream/.git','--exclude=./build','-I','gzip -1','-cf',archive,'-C',root,'.'],root.with_suffix('.archive.log'),deadline)
            print(str(archive))
        except Exception as e:
            status['archive_error']=str(e);status['status']='archive_incomplete';exitcode=124 if isinstance(e,TimeoutError) else 1
        status['total_elapsed_seconds']=time.monotonic()-started;status['exit_code']=exitcode
        signal.setitimer(signal.ITIMER_REAL,0)
        save(root/'final-status.json',status)
        save(root.with_suffix('.status.json'),status)
        print(json.dumps(status));print('Retained directory: '+str(root))
    return exitcode

if __name__=='__main__':sys.exit(main())
