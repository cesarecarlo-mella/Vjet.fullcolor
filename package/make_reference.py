#!/usr/bin/env python3
"""write the reference data of the package tests (run in the repository, before build_dist.py):
  tests/lc_reference.json      the Nc^3 layer (HA1-6) at 40 of the 406 grid points of
                               reproduceLC/HA_series_Nc3_vs_exact_digits: exact values of the published planar results
                               (GPL form, >= 15 digits), independent of this evaluator
  tests/reference_values.json  all 84 remainders and their error estimates at 12 points (including points that need
                               another path deformation or the average across a spurious curve), from this release"""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
sys.path.insert(0, os.path.join(ROOT, 'iter int ', 'generators_deq'))
import remainder_eval    # noqa: E402

grid = json.load(open(os.path.join(ROOT, 'reproduceLC', 'HA_series_Nc3_vs_exact_digits', 'results', 'digits_order18.json')))['rows']
keep = [r for i, r in enumerate(grid) if i % 10 == 0]
for xy in ((1 / 30, 1 / 30), (1 / 30, 28 / 30), (28 / 30, 1 / 30), (7 / 30, 19 / 30)):    # corners, a path switch
    keep += [r for r in grid if abs(r['x'] - xy[0]) < 1e-12 and abs(r['y'] - xy[1]) < 1e-12]
json.dump(dict(source='reproduceLC/HA_series_Nc3_vs_exact_digits (Nc^3 layer of the published planar decay.m, '
                      'alpha[i,3], gamma[i,3], q2 = 1); values [Re, Im] for alpha1..3, gamma1..3',
               points=[[r['x'], r['y']] for r in keep],
               values=[[[float(a), float(b)] for a, b in r['exact']] for r in keep]),
          open(os.path.join(HERE, 'tests', 'lc_reference.json'), 'w'), indent=1)
print('lc_reference.json: %d points' % len(keep))

pts = [(0.3, 0.3), (0.25, 0.35), (0.1, 0.7), (0.6, 0.2), (0.05, 0.05), (0.45, 0.5), (0.02, 0.5), (0.8, 0.1),
       (7 / 30, 19 / 30), (0.36, 0.24), (0.21, 0.7), (0.16, 0.6)]
R, info, err = remainder_eval.evaluate_remainders(pts, 1e-12, error=True)
names = list(remainder_eval.NAMES)
layers = [k for n, k in remainder_eval.KEYS if n == names[0]]
col = np.array([[remainder_eval.KEYS.index((n, k)) for k in layers] for n in names])
V, E = R[:, col], err[:, col]
json.dump(dict(coefficients=names, layers=layers, points=[list(p) for p in pts],
               values=[[[[V[p, i, l].real, V[p, i, l].imag] for l in range(len(layers))] for i in range(6)] for p in range(len(pts))],
               errors=[[[float(E[p, i, l]) for l in range(len(layers))] for i in range(6)] for p in range(len(pts))],
               notes=['averaged across %s' % i['averaged']['letter'] if i.get('averaged') else
                      ('path deformation %s' % (i['delta'],) if i['delta'] != (0.1, 0.2) else '') for i in info]),
          open(os.path.join(HERE, 'tests', 'reference_values.json'), 'w'), indent=1)
print('reference_values.json: %d points' % len(pts), [n for n in ['averaged' if i.get('averaged') else '' for i in info]])
