"""Frozen-family exact inference; one Bernoulli trial per interaction/category."""
import math
from scipy.stats import beta, binom, fisher_exact

ALPHA = .01

def holm(values, alpha=ALPHA):
    order = sorted(range(len(values)), key=values.__getitem__)
    adjusted = [1.] * len(values)
    previous = 0.
    for rank, i in enumerate(order):
        previous = max(previous, min(1., (len(values)-rank)*values[i]))
        adjusted[i] = previous
    return adjusted

def interval(k, n, alpha):
    return (0. if k == 0 else float(beta.ppf(alpha/2, k, n-k+1)),
            1. if k == n else float(beta.isf(alpha/2, k+1, n-k)))

def composite_p(k, n, lo, hi):
    """Conservative two-sided equal-tail test of a numerical null interval."""
    if lo <= k/n <= hi:
        return 1.
    return min(1., 2*float(binom.cdf(k, n, lo) if k/n < lo else binom.sf(k-1, n, hi)))

def read_reference(path):
    rows = {}
    for line in path.read_text().splitlines():
        if line.startswith('#'):
            continue
        t = line.split()
        rows[t[0]] = dict(id=t[0], view=t[1], kind=t[2], lo=float(t[3]), hi=float(t[4]),
                          p=list(map(float,t[5:8])), possible=list(map(int,t[8:11])))
    return rows

def freeze_reference(paths):
    levels = [read_reference(p) for p in paths]
    if any(set(x) != set(levels[0]) for x in levels):
        raise ValueError('Reference observable families differ')
    frozen = {}
    for key, row in levels[-1].items():
        differences = [[abs(levels[j][key]['p'][c]-levels[j-1][key]['p'][c])
                        for c in range(3)] for j in range(1,len(levels))]
        error = [max(4*max(d[c] for d in differences), 2e-14) for c in range(3)]
        row = dict(row, error=error)
        allowed=[c for c in range(3) if row['possible'][c]]
        if len(allowed)==1:
            row['p']=[float(c==allowed[0]) for c in range(3)]
            row['error']=[0.,0.,0.]
        row['mean'] = row['p'][1]+2*row['p'][2]
        row['mean_error'] = row['error'][1]+2*row['error'][2]
        # The bin family is classified using the reference only, before MC.
        row['target'] = .1 if row['kind']=='tail' else .01
        frozen[key] = row
    # Classify bins wholly outside central 98% of each view as tails.
    for view in {r['view'] for r in frozen.values()}:
        bins = sorted((r for r in frozen.values() if r['view']==view and r['kind']=='bin'), key=lambda r:r['lo'])
        total = sum(r['mean'] for r in bins)
        cdf = 0.
        for r in bins:
            end = cdf+r['mean']
            if total and (end <= .01*total or cdf >= .99*total):
                r['target'] = .1
            cdf = end
    return frozen

def read_counts(paths):
    counts = {}
    for path in paths:
        for line in path.read_text().splitlines():
            if line.startswith('#'):
                continue
            t=line.split(); a=list(map(int,t[1:]))
            if len(a)!=3 or min(a)<0:
                raise ValueError('Invalid event multiplicity counts')
            counts[t[0]] = [x+y for x,y in zip(counts.get(t[0],[0,0,0]),a)]
    return counts

def analyze(references, counts, planned_n, complete):
    # Keep even structurally zero cells in the frozen family. No data-dependent
    # hypothesis removal, merging, selection, or stopping is allowed.
    family = sum(len(r)*9 for r in references.values())  # 2 oracle + 1 old/new, 3 categories
    if not family:
        raise ValueError('Empty frozen hypothesis family')
    for case,ref in references.items():
        for mode in ('old','current'):
            if set(counts[case][mode])!=set(ref):
                raise ValueError('Counts changed the frozen observable family')
            if any(sum(counts[case][mode][key])!=planned_n[case] for key in ref):
                complete=False
    ci_alpha = ALPHA / family
    tests, spectra = [], []
    for case, ref in references.items():
        old, new = counts[case]['old'], counts[case]['current']
        for key, r in ref.items():
            if key not in old or key not in new:
                raise ValueError('Missing frozen observable')
            ns = {mode:sum(counts[case][mode][key]) for mode in ('old','current')}
            if min(ns.values())<=0:
                raise ValueError('Empty sample')
            for c in range(3):
                lo=max(0.,r['p'][c]-r['error'][c]); hi=min(1.,r['p'][c]+r['error'][c])
                for mode in ('old','current'):
                    k=counts[case][mode][key][c]; n=ns[mode]
                    p=composite_p(k,n,lo,hi)
                    tests.append(dict(case=case,id=key,category=c,comparison=mode+'_oracle',p=p))
                a,b=old[key][c],new[key][c]
                p=float(fisher_exact([[a,ns['old']-a],[b,ns['current']-b]]).pvalue)
                tests.append(dict(case=case,id=key,category=c,comparison='old_current',p=p))
            cis={}
            for mode in ('old','current'):
                n=ns[mode]; cc=counts[case][mode][key]
                bounds=[interval(cc[c],n,ci_alpha) for c in (1,2)]
                cis[mode]=(bounds[0][0]+2*bounds[1][0],bounds[0][1]+2*bounds[1][1])
                mu=r['mean']; err=r['mean_error']
                numeric_ok = (mu==0 and not r['possible'][1] and not r['possible'][2]) or (mu>0 and err <= .05*r['target']*mu)
                # Numerical refinement is an empirical guard, not a proof bound.
                if mu>0:
                    var=r['p'][1]+4*r['p'][2]-mu*mu
                    numeric_ok &= err <= .1*math.sqrt(max(0.,var)/n)
                lower,upper=cis[mode]
                effect = None if mu==0 else [(lower-mu-err)/mu,(upper-mu+err)/mu]
                within=effect is not None and max(abs(x) for x in effect)<=r['target']
                material=effect is not None and (effect[0]>r['target'] or effect[1]<-r['target'])
                forbidden=mu==0 and not r['possible'][1] and not r['possible'][2]
                if forbidden:
                    within=(cc[1]+cc[2]==0);material=not within
                state='inconclusive'
                if complete and n==planned_n[case] and numeric_ok:
                    state='deviation' if material else ('within_target' if within else 'inconclusive')
                spectra.append(dict(case=case,id=key,sampler=mode,N=n,counts=cc,mean=(cc[1]+2*cc[2])/n,
                    reference=mu,reference_error=err,mean_ci=[lower,upper],relative_effect_ci=effect,
                    target=r['target'],numerical_guard=bool(numeric_ok),status=state))
            mu=r['mean'];err=r['mean_error']
            difference=[cis['old'][0]-cis['current'][1],cis['old'][1]-cis['current'][0]]
            effect=None if mu==0 else [x/mu for x in difference]
            state='inconclusive'
            if complete and all(n==planned_n[case] for n in ns.values()):
                if mu>0 and err <= .05*r['target']*mu:
                    state='within_target' if max(abs(x) for x in effect)<=r['target'] else ('deviation' if effect[0]>r['target'] or effect[1]<-r['target'] else 'inconclusive')
                elif mu==0 and not r['possible'][1] and not r['possible'][2]:
                    state='within_target' if old[key][1]+old[key][2]+new[key][1]+new[key][2]==0 else 'deviation'
            spectra.append(dict(case=case,id=key,sampler='old_minus_current',N=ns,
                mean_ci=difference,relative_effect_ci=effect,reference=mu,target=r['target'],status=state))
    if len(tests)!=family:
        raise ValueError('Changed hypothesis family')
    adjusted=holm([t['p'] for t in tests])
    for t,p in zip(tests,adjusted):
        t['holm_p']=p;t['significant']=bool(complete and p<=ALPHA)
    states=[r['status'] for r in spectra if 'status' in r]
    status='deviation' if 'deviation' in states else ('within_target' if states and all(s=='within_target' for s in states) else 'inconclusive')
    return dict(alpha=ALPHA,family_size=family,complete=complete,status=status,tests=tests,spectra=spectra,
                interval_method='simultaneous Bonferroni Clopper-Pearson, summed multiplicity bounds',
                significance_method='Holm, exact binomial composite oracle / Fisher old-current')
