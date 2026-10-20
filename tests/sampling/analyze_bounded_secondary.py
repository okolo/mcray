#!/usr/bin/env python3
"""Joint comparison against two controls; Bonferroni allocation across families."""
import json
from collections import Counter
from pathlib import Path
import sys
from analyze_secondary import write_report


def analyze_run(data, output):
    output=Path(output);controls=data.get('controls',['old'])
    if not controls or len(set(controls))!=len(controls) or set(controls)-{'old','cdf'}:
        raise ValueError('Invalid frozen control family')
    reports={}
    labels=data['sampler_labels']
    for control in controls:
        selected={k:data[k] for k in ('references','planned_n','complete')}
        selected['counts']={case:{'old':values[control],'current':values['current']} for case,values in data['counts'].items()}
        selected['sampler_labels']={'old':labels[control],'current':labels['current'],'old_minus_current':labels[control]+' minus '+labels['current']}
        selected['alpha']=.01/len(controls)
        reports[control]=write_report(selected,output.with_name(output.stem+'-'+control+'.json'))
    tests=[];spectra=[]
    for control,report in reports.items():
        for row in report['tests']:
            row=dict(row,comparison_family=control)
            row['holm_p']=min(1.,len(controls)*row['holm_p'])
            tests.append(row)
        spectra.extend(dict(row,comparison_family=control) for row in report['spectra'])
    states=[r['status'] for r in reports.values()]
    status='deviation' if 'deviation' in states else 'within_target' if all(s=='within_target' for s in states) else 'inconclusive'
    result=dict(alpha=.01,alpha_allocation={c:.01/len(controls) for c in controls},
                family_size=sum(r['family_size'] for r in reports.values()),complete=data['complete'],
                status=status,tests=tests,spectra=spectra,
                precision_states=dict(Counter(r['status'] for r in spectra)),
                significant_comparisons=dict(Counter(r['comparison_family']+':'+r['comparison'] for r in tests if r['significant'])),
                pair_report_files={c:output.stem+'-'+c+'.json' for c in controls},
                significance_method='Bonferroni across control families; Holm and exact tests inside each family',
                interval_method='Bonferroni across families and original simultaneous exact intervals within each')
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    lines=['# Bounded secondary sampling comparison','',f'Result: **{status}**. Global significance level: 1%.',
           'Precision targets remain 1% in core bins and 10% in tails.','',
           'The bounded sample is reused for both comparisons; the multiplicity correction does not require independence.','']
    for control,report in reports.items():
        lines.append(f"- {labels[control]}: {report['status']}; allocated alpha={report['alpha']:.8g}; [full bin/tail report]({output.stem}-{control}.md).")
    output.with_suffix('.md').write_text('\n'.join(lines)+'\n')
    return result


if __name__=='__main__':
    analyze_run(json.loads(Path(sys.argv[1]).read_text()),sys.argv[2])
