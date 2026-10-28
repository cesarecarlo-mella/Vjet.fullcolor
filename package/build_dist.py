#!/usr/bin/env python3
"""assemble the shippable package dist/Zgqq_evaluator-<VERSION>/ and dist/Zgqq_evaluator-<VERSION>.tar.gz from the
repository: the user files of package/, the evaluator of iter int /generators_deq (code, data, C++ source) and the
amplitudes of Zgqq_finite_remainders.  Run make_reference.py first when the data change."""
import hashlib, os, shutil, tarfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
G = os.path.join(ROOT, 'iter int ', 'generators_deq')
Z = os.path.join(ROOT, 'Zgqq_finite_remainders')
VERSION = open(os.path.join(HERE, 'VERSION')).read().strip()
NAME = 'Zgqq_evaluator-' + VERSION
DIST = os.path.join(ROOT, 'dist')
DST = os.path.join(DIST, NAME)

FILES = {                                            # destination -> source
    'README.md': os.path.join(HERE, 'README.md'),
    'Makefile': os.path.join(HERE, 'Makefile'),
    'VERSION': os.path.join(HERE, 'VERSION'),
    'zgqq_eval.py': os.path.join(HERE, 'zgqq_eval.py'),
    'tests/run_tests.py': os.path.join(HERE, 'tests', 'run_tests.py'),
    'tests/lc_reference.json': os.path.join(HERE, 'tests', 'lc_reference.json'),
    'tests/reference_values.json': os.path.join(HERE, 'tests', 'reference_values.json'),
    'core/gen_eval.py': os.path.join(G, 'gen_eval.py'),
    'core/remainder_eval.py': os.path.join(G, 'remainder_eval.py'),
    'core/out/meta.json': os.path.join(G, 'out', 'meta.json'),
    'core/out/generators_bs.bin': os.path.join(G, 'out', 'generators_bs.bin'),
    'core/out/start_y.head': os.path.join(G, 'out', 'start_y.head'),
    'core/cpp/tridents_bs.cpp': os.path.join(G, 'cpp', 'tridents_bs.cpp'),
    'core/cpp/third_party.tar.gz': os.path.join(G, 'cpp', 'third_party.tar.gz'),
    'amplitudes/Alphabet20.wl': os.path.join(Z, 'Alphabet20.wl'),
    'amplitudes/Generators.wl': os.path.join(Z, 'Generators.wl'),
    'amplitudes/RationalFunctionBasis.wl': os.path.join(Z, 'RationalFunctionBasis.wl'),
    'amplitudes/TranscendentalMultipliers.wl': os.path.join(Z, 'TranscendentalMultipliers.wl'),
}
for n in ('alpha1', 'alpha2', 'alpha3', 'gamma1', 'gamma2', 'gamma3'):
    FILES['amplitudes/FiniteRemainders/%s.wl' % n] = os.path.join(Z, 'FiniteRemainders', n + '.wl')

if os.path.isdir(DST):
    shutil.rmtree(DST)
for dst, src in FILES.items():
    os.makedirs(os.path.dirname(os.path.join(DST, dst)), exist_ok=True)
    shutil.copy2(src, os.path.join(DST, dst))
with open(os.path.join(DST, 'MANIFEST.sha256'), 'w') as f:
    for dst in sorted(FILES):
        f.write('%s  %s\n' % (hashlib.sha256(open(os.path.join(DST, dst), 'rb').read()).hexdigest(), dst))
tgz = os.path.join(DIST, NAME + '.tar.gz')
with tarfile.open(tgz, 'w:gz') as t:
    t.add(DST, arcname=NAME)
print('%s: %d files, %.1f MB; %s: %.1f MB' % (DST, len(FILES) + 1, sum(os.path.getsize(os.path.join(DST, d)) for d in FILES) / 1e6,
                                              tgz, os.path.getsize(tgz) / 1e6))
