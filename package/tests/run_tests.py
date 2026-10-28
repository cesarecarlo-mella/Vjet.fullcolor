#!/usr/bin/env python3
"""installation test of Zgqq_evaluator (make test):
  1. the Nc^3 layer at 40 points against exact values: the published planar three-loop results (GPL form, >= 15 digits),
     independent of this evaluator;
  2. all 84 remainders at 12 points against the reference values of this release (same data, possibly another
     compiler or machine), within the error estimates;
  3. the command-line interface."""
import json, os, subprocess, sys, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
import zgqq_eval as Z    # noqa: E402

ok = True

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
print('1. Nc^3 layer vs exact planar values, %d points x 6 coefficients: worst relative difference %.1e, '
      'difference within the error estimate for %d of %d values  %s'
      % (len(pts), rel.max(), inside.sum(), inside.size, 'ok' if t1 else 'FAILED'))

# 2. all 84 vs the release reference
ref = json.load(open(os.path.join(HERE, 'reference_values.json')))
assert ref['coefficients'] == Z.COEFFICIENTS and ref['layers'] == Z.LAYERS
pts = [tuple(p) for p in ref['points']]
Vr = np.array([[[complex(*v) for v in lay] for lay in co] for co in ref['values']])
Er = np.array(ref['errors'])
res = Z.evaluate(pts)
diff = np.abs(res.values - Vr)
inside = diff <= res.errors + Er + 1e-13 * np.abs(Vr)
t2 = inside.all()
ok &= t2
print('2. all 84 remainders vs the release reference, %d points: worst relative difference %.1e, '
      'within the error estimates for %d of %d values  %s'
      % (len(pts), (diff / np.maximum(np.abs(Vr), 1)).max(), inside.sum(), inside.size, 'ok' if t2 else 'FAILED'))

# 3. command line
with tempfile.TemporaryDirectory() as d:
    out = os.path.join(d, 'r.csv')
    r = subprocess.run([sys.executable, os.path.join(HERE, '..', 'zgqq_eval.py'), '0.3', '0.3', '--out', out],
                       capture_output=True, text=True)
    lines = open(out).read().split('\n') if r.returncode == 0 else []
    t3 = r.returncode == 0 and len([x for x in lines if x]) == 1 + 84
    ok &= t3
    print('3. command line (csv output)  %s' % ('ok' if t3 else 'FAILED ' + r.stderr[-300:]))

print('PASSED' if ok else 'FAILED')
sys.exit(0 if ok else 1)
