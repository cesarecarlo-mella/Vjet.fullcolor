#!/usr/bin/env python3
"""installation test of Zgqq_evaluator (make test):
  1. decay: the Nc^3 layer at 45 points against exact values: the published planar three-loop results (GPL form,
     >= 15 digits), independent of this evaluator;
  2. decay: all 84 remainders at 12 points against the reference values of this release (same data, possibly another
     compiler or machine), within the error estimates;
  3. production regions x, y, z: all 84 remainders at 3 points each against the reference values of this release;
  4. production region y (q qbar-type channel qbar g -> V qbar): the layers Nc^3, Nc^2 Nf, Nc Nf^2 against the
     published planar production results (lmL.m, lpL.m, GiNaC at 25 digits), independent of this evaluator;
  5. the command-line interface.
run_tests.py --quad (make test-quad): the quad precision option instead:
  6. three decay points and one point in each production region in quad: the Nc^3 layer against the exact planar
     values (decay), all 84 against the double evaluation within its error estimate, and the quad error estimate
     itself (two paths) far below double precision."""
import json, os, subprocess, sys, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
import zgqq_eval as Z    # noqa: E402

ok = True
QUAD = '--quad' in sys.argv
cplx = lambda a: np.array([[[complex(*v) for v in lay] for lay in co] for co in a])

if QUAD:
    lc = json.load(open(os.path.join(HERE, 'lc_reference.json')))
    sel = [0, len(lc['points']) // 2, len(lc['points']) - 1]
    pts = [tuple(lc['points'][k]) for k in sel]
    E = np.array([[complex(*v) for v in lc['values'][k]] for k in sel])
    rq, rd = Z.evaluate(pts, precision='quad'), Z.evaluate(pts)
    l = Z.LAYERS.index('Nc^3*Nf^0*NfZ^0')
    rel_lc = (np.abs(rq.values[:, :, l] - E) / np.abs(E)).max()
    ins, tot, qrel = [], 0, []
    for reg, P in [('decay', pts)] + [(A, [(0.2, 0.5)]) for A in ('x', 'y', 'z')]:
        q, d = (rq, rd) if reg == 'decay' else (Z.evaluate(P, precision='quad', region=reg), Z.evaluate(P, region=reg))
        ins.append(np.abs(q.values - d.values) <= d.errors + 1e-13 * np.abs(q.values))
        qrel.append(q.errors / np.maximum(np.abs(q.values), 1e-300))
    inside = np.concatenate([a.ravel() for a in ins])
    qrel = np.median(np.concatenate([a.ravel() for a in qrel]))
    t6 = rel_lc < 1e-12 and inside.mean() > 0.99 and qrel < 1e-20
    print('6. quad, 3 decay points + 1 point in each production region: Nc^3 layer vs exact planar values %.1e '
          '(limited by the >= 15 digits of the reference); double within its error estimate for %d of %d values; '
          'median relative quad error estimate %.0e  %s'
          % (rel_lc, inside.sum(), inside.size, qrel, 'ok' if t6 else 'FAILED'))
    print('PASSED' if t6 else 'FAILED')
    sys.exit(0 if t6 else 1)

# 1. leading colour vs exact
lc = json.load(open(os.path.join(HERE, 'lc_reference.json')))
pts = [tuple(p) for p in lc['points']]
E = np.array([[complex(*v) for v in row] for row in lc['values']])
res = Z.evaluate(pts)
l = Z.LAYERS.index('Nc^3*Nf^0*NfZ^0')
V, Er = res.values[:, :, l], res.errors[:, :, l]
rel = np.abs(V - E) / np.abs(E)
inside = np.abs(V - E) <= Er + 1e-13 * np.abs(E)
t1 = rel.max() < 1e-8 and inside.mean() > 0.95
ok &= t1
print('1. decay: Nc^3 layer vs exact planar values, %d points x 6 coefficients: worst relative difference %.1e, '
      'difference within the error estimate for %d of %d values  %s'
      % (len(pts), rel.max(), inside.sum(), inside.size, 'ok' if t1 else 'FAILED'))

# 2. all 84 vs the release reference
ref = json.load(open(os.path.join(HERE, 'reference_values.json')))
assert ref['coefficients'] == Z.COEFFICIENTS and ref['layers'] == Z.LAYERS
pts = [tuple(p) for p in ref['points']]
Vr, Er = cplx(ref['values']), np.array(ref['errors'])
res = Z.evaluate(pts)
diff = np.abs(res.values - Vr)
inside = diff <= res.errors + Er + 1e-13 * np.abs(Vr)
t2 = inside.all()
ok &= t2
print('2. decay: all 84 remainders vs the release reference, %d points: worst relative difference %.1e, '
      'within the error estimates for %d of %d values  %s'
      % (len(pts), (diff / np.maximum(np.abs(Vr), 1)).max(), inside.sum(), inside.size, 'ok' if t2 else 'FAILED'))

# 3. production regions vs the release reference
ref = json.load(open(os.path.join(HERE, 'region_reference_values.json')))
worst, n_in, n_all = 0, 0, 0
for A, d in ref['regions'].items():
    pts = [tuple(p) for p in d['points']]
    Vr, Er = cplx(d['values']), np.array(d['errors'])
    res = Z.evaluate(pts, region=A)
    diff = np.abs(res.values - Vr)
    inside = diff <= res.errors + Er + 1e-13 * np.abs(Vr)
    worst = max(worst, (diff / np.maximum(np.abs(Vr), 1)).max()); n_in += inside.sum(); n_all += inside.size
t3 = n_in == n_all
ok &= t3
print('3. production regions x, y, z: all 84 remainders vs the release reference, 3 points each: worst relative '
      'difference %.1e, within the error estimates for %d of %d values  %s' % (worst, n_in, n_all, 'ok' if t3 else 'FAILED'))

# 4. region y vs the published planar production results
pub = json.load(open(os.path.join(HERE, 'region_lc_reference.json')))
pts = [tuple(p) for p in pub['points']]
res = Z.evaluate(pts, region='y')
rel = []
for e in pub['values']:
    v = res.values[e['point'], Z.COEFFICIENTS.index(e['coefficient']), Z.LAYERS.index(e['layer'])]
    rel.append(abs(v - complex(*e['value'])) / abs(complex(*e['value'])))
t4 = max(rel) < 1e-8
ok &= t4
print('4. region y vs the published planar production results (layers Nc^3, Nc^2 Nf, Nc Nf^2), %d points: worst '
      'relative difference %.1e  %s' % (len(pts), max(rel), 'ok' if t4 else 'FAILED'))

# 5. command line
with tempfile.TemporaryDirectory() as d:
    out = os.path.join(d, 'r.csv')
    r = subprocess.run([sys.executable, os.path.join(HERE, '..', 'zgqq_eval.py'), '0.3', '0.3', '--out', out],
                       capture_output=True, text=True)
    lines = open(out).read().split('\n') if r.returncode == 0 else []
    r2 = subprocess.run([sys.executable, os.path.join(HERE, '..', 'zgqq_eval.py'), '--physical', '--no-error', '--',
                         '-0.75', '-0.75'], capture_output=True, text=True)
    t5 = r.returncode == 0 and len([x for x in lines if x]) == 1 + 84 and r2.returncode == 0 and 'region z' in r2.stdout
    ok &= t5
    print('5. command line (csv output, physical point in a production region)  %s' % ('ok' if t5 else 'FAILED ' + r.stderr[-300:] + r2.stderr[-300:]))

print('PASSED' if ok else 'FAILED')
sys.exit(0 if ok else 1)
