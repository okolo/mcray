#!/usr/bin/env python3
"""Extract the historical sampler locally; never redistribute Numerical Recipes."""
import argparse
import hashlib
import json
import pathlib
import subprocess

OLD = 'd05723f449c1eecee4f00af392014df4c898ebe9'
p = argparse.ArgumentParser()
p.add_argument('--repo', type=pathlib.Path, required=True)
p.add_argument('--nr', type=pathlib.Path, required=True)
p.add_argument('--output', type=pathlib.Path, required=True)
a = p.parse_args()
a.output.mkdir(parents=True, exist_ok=True)
old = subprocess.check_output(['git', '-C', str(a.repo), 'show', OLD + ':src/lib/MathUtils.cpp'], text=True)
start = old.index('class MathUtilsLogODE')
end = old.index('bool MathUtils::SampleLogDistributionBoost', start)
reference = old[start:end]
wrapper_start = reference.index('bool MathUtils::SampleLogDistribution(')
wrapper_end = reference.index('bool MathUtils::SampleLogDistributionNR(')
reference = reference[:wrapper_start] + reference[wrapper_end:]
reference = reference.replace('MathUtils::SampleLogDistributionNR', 'HistoricalSampleLogDistribution')
logger_start = reference.index('\tSafePtr<DebugFunctionX<double> > func;')
logger_end = reference.index('\n\tdouble distrXmin', logger_start)
reference = reference[:logger_start] + reference[logger_end:]
headers = ['nr3.h', 'stepper.h', 'odeint.h', 'stepperdopr5.h']
expected = json.loads((pathlib.Path(__file__).parent / 'nr_reference_hashes.json').read_text())
actual_nr = subprocess.check_output(['git', '-C', str(a.repo), 'rev-parse', OLD + ':src/external/nr'], text=True).strip()
assert actual_nr == expected['commit'], 'Historical NR gitlink mismatch'
for name in headers:
    assert hashlib.sha256((a.nr/name).read_bytes()).hexdigest() == expected['headers'][name], 'Wrong historical NR header: ' + name

prefix = '#include "MathUtils.h"\n' + ''.join('#include "' + h + '"\n' for h in headers)
(a.output / 'historical_sampler.cpp').write_text(prefix + 'using namespace Utils;\n' + reference)
background = (a.repo / 'src/lib/Background.cpp').read_text()
needle = 'math.SampleLogDistributionBounded(*kern, aRand, aS, totalRate, Smin, Smax, fEpsRel)'
assert background.count(needle) == 1, 'Sampling call changed; review adapter'
background = background.replace(needle, 'SamplingProbe(*kern, aRand, aS, totalRate, Smin, Smax, fEpsRel)')
background = '#include "MathUtils.h"\nextern bool SamplingProbe(const Utils::Function&, mcray::Randomizer&, double&, double&, double, double, double);\n' + background
(a.output / 'probe_background.cpp').write_text(background)
# Current production code is linked unchanged. Legacy rejection is a frozen
# validation-only control with the same corrected conditional secondary laws.
LEGACY = '1cc28f381e961f3d4947d08f0149ac72d63c3635'
math = (a.repo / 'src/lib/MathUtils.cpp').read_text()
CDF = 'fa562d207b58af6d1357220e3809b44e0e2aa1ea'
assert (a.repo/'src/lib/NumericalLogCDF.h').read_bytes() == subprocess.check_output(
    ['git','-C',str(a.repo),'show',CDF+':src/lib/NumericalLogCDF.h']), 'CDF control header changed'
frozen_math = subprocess.check_output(['git','-C',str(a.repo),'show',CDF+':src/lib/MathUtils.cpp'],text=True)
def cdf_body(source):
    start=source.index('bool MathUtils::SampleLogDistributionCDF(')
    return source[start:source.index('\n}',start)+2]
assert cdf_body(math)==cdf_body(frozen_math), 'CDF control function changed'
(a.output / 'audited_math.cpp').write_text(math)
legacy = subprocess.check_output(['git','-C',str(a.repo),'show',LEGACY+':src/lib/MathUtils.cpp'],text=True)
start = legacy.index('class LogscaleDistr :')
end = legacy.index('template<typename X>',start)
body = legacy[start:end]
body = body.replace('bool MathUtils::SampleLogDistribution(', 'bool LegacyRejectionSampleLogDistribution(')
body = body.replace('return SampleLogDistribution(aDistrib,', 'return LegacyRejectionSampleLogDistribution(aDistrib,')
body = body.replace('aOutputIntegral = Integration_qag (', 'thread_local MathUtils math;\n    aOutputIntegral = math.Integration_qag (')
def replace_once(old, new):
    global body
    assert body.count(old)==1, 'Frozen rejection control changed'
    body=body.replace(old,new)
replace_once('double f = pdf(x);', 'const double previous_max = max_val;\n        ++secondary::audit.proposals;\n        double f = pdf(x);\n        if (f > previous_max) {\n            ++secondary::audit.envelope;\n            secondary::audit.largestEnvelopeRatio = std::max(secondary::audit.largestEnvelopeRatio, f/previous_max);\n        }')
replace_once('if (max_recurse <= 0){', 'if (max_recurse <= 0){\n        ++secondary::audit.argmax;')
replace_once('return LegacyRejectionSampleLogDistribution(aDistrib,', '++secondary::audit.narrowing;\n    return LegacyRejectionSampleLogDistribution(aDistrib,')
(a.output/'legacy_rejection.cpp').write_text('#include "MathUtils.h"\n#include "secondary_oracle.h"\nusing namespace Utils;\n'+body)
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
manifest = {'old_commit': OLD, 'legacy_rejection_commit': LEGACY, 'current_commit': subprocess.check_output(['git', '-C', str(a.repo), 'rev-parse', 'HEAD'], text=True).strip(),
            'nr_files': {h: sha(a.nr / h) for h in headers},
            'historical_nr_gitlink': subprocess.check_output(['git', '-C', str(a.repo), 'rev-parse', OLD + ':src/external/nr'], text=True).strip(),
            'sources': {str(path.relative_to(a.repo)): sha(path) for path in [a.repo / name for name in subprocess.check_output(['git', '-C', str(a.repo), 'ls-files', 'src/lib'], text=True).splitlines() if (a.repo / name).is_file()] + [a.repo / 'tests/sampling/probe.cpp', a.repo / 'tests/sampling/prepare_reference.py']},
            'generated': {path.name: sha(path) for path in a.output.glob('*.cpp')}}
(a.output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
