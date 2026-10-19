#!/usr/bin/env python3
"""Reanalyze frozen inputs and write machine-readable and readable reports."""
import csv
import json
from collections import Counter
from pathlib import Path
import sys
from secondary_statistics import analyze
if __name__=='__main__':
    data=json.loads(Path(sys.argv[1]).read_text())
    labels=data.pop('sampler_labels',{'old':'Historical NR','current':'Current sampler','old_minus_current':'Old minus current'})
    output=Path(sys.argv[2]);report=analyze(**data)
    report['sampler_labels']=labels
    report['significant_comparisons']=dict(Counter(t['comparison'] for t in report['tests'] if t['significant']))
    report['precision_states']=dict(Counter(r['status'] for r in report['spectra'] if 'status' in r))
    output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    with output.with_suffix('.tsv').open('w') as f:
        columns=['case','id','sampler','N','counts','mean','reference','reference_error','mean_ci','relative_effect_ci','target','numerical_guard','status']
        writer=csv.DictWriter(f,fieldnames=columns,delimiter='\t',extrasaction='ignore');writer.writeheader();writer.writerows(report['spectra'])
    lines=['# Secondary sampling comparison','',
        'Result: **'+report['status']+'**. Complete frozen acceptance ensemble: '+str(report['complete'])+'.','',
        'Global significance level: 1%; hypotheses: '+str(report['family_size'])+'.',
        'Precision targets: 1% in core bins, 10% in tails. All interactions are independent trials; paired products remain correlated.',
        'Reference convergence is an empirical numerical guard, not a certified error bound.','',
        '| Case | Sampler | Within target | Material deviations | Inconclusive | Significant category tests |',
        '| --- | --- | ---: | ---: | ---: | ---: |']
    for case in data['references']:
        for mode in ('old','current','old_minus_current'):
            counts=Counter(r['status'] for r in report['spectra'] if r['case']==case and r['sampler']==mode)
            significant=sum(t['significant'] for t in report['tests'] if t['case']==case and t['comparison']==('old_current' if mode=='old_minus_current' else mode+'_oracle'))
            lines.append(f"| {case} | {labels[mode]} | {counts['within_target']} | {counts['deviation']} | {counts['inconclusive']} | {significant} |")
    lines+=['','See report.tsv for every frozen bin/tail and simultaneous confidence interval; report.json includes all Holm-adjusted tests and direct old/current comparisons.',
            'Zero observed tail events retain a positive exact upper confidence limit. Incomplete runs and smoke runs cannot establish equivalence.','']
    output.with_suffix('.md').write_text('\n'.join(lines))
