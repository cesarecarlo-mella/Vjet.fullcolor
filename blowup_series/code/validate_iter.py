#!/usr/bin/env python3
"""the ITER family (iterated-integral system) in blow-up variables against the numerical evaluator: the series of the
functions (exact_<tag>_w<w>.m, weights 1..6, plus sol0) are combined into the 579 generators (generators.wl) and
evaluated at points close to the corner of each chart; gen_eval gives the same generators by numerical integration.
usage: validate_iter.py <ITER dist dir> <ITER inputs dir> <generators.wl> [charts, default 1,3,5,2]"""
import math, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', '..', 'iter int ', 'generators_deq'))
import mparse as mp      # noqa: E402
import gen_eval          # noqa: E402

dist, inputs, genfile = sys.argv[1:4]
charts = [int(c) for c in (sys.argv[4] if len(sys.argv) > 4 else '1,3,5,2').split(',')]
labels = re.findall(r'"([^"]+)"', open(os.path.join(inputs, 'labels.m')).read())
idx = {l: i for i, l in enumerate(labels)}
gens = []
for line in open(genfile):
    m = re.match(r'(T\[\d+\]) -> (.*?),?\s*$', line.strip())
    if m:
        gens.append((m.group(1), [(float(eval(c.replace('/', '/'))), idx[l]) for c, l in re.findall(r'\((-?[\d/]+)\)\*([A-Z]\[[\d,]+\])', m.group(2))]))
sol0 = [float(eval(s.strip())) for s in mp.split_top_list(open(os.path.join(inputs, 'sol0.m')).read())]
CONST = {'Pi': math.pi, 'Zeta[3]': 1.2020569031595942854, 'Zeta[5]': 1.0369277551433699263, 'I': 1j}
# charts: (BIG name, map (t, B) -> (x, y))
MAP = {1: ('y', lambda t, B: (t * B * B, B)), 2: ('y', lambda t, B: (1 - B - t * B * B, B)),
       3: ('z', lambda t, B: (t * B * B, 1 - B - t * B * B)), 5: ('z', lambda t, B: (1 - B - t * B * B, t * B * B))}
for c in charts:
    big, f = MAP[c]
    ser = [mp.parse_file(os.path.join(dist, 'exact_c%d_w%d.m' % (c, w))) for w in range(1, 7)]
    for tv, Bv in ((0.01, 0.04), (0.1, 0.1)):
        val = {'t': tv, big: Bv, 'Log[t]': math.log(tv), 'Log[%s]' % big: math.log(Bv)}
        val.update(CONST)
        F = np.array(sol0, complex)
        for S in ser:
            for i, poly in enumerate(S):
                s = 0j
                for mono, coef in poly.items():
                    term = float(coef)
                    for sym, e2 in mono:
                        term = term * val[sym] ** (e2 / 2)
                    s += term
                F[i] += s
        T = np.array([sum(a * F[i] for a, i in comb) for _, comb in gens])
        P = f(tv, Bv)
        N, info = gen_eval.evaluate([P], 1e-12)
        names = [n for n, _ in gens]
        order = [gen_eval.NAMES.index(n) for n in names]
        d = np.abs(T - N[0][order]) / np.maximum(np.abs(N[0][order]), 1)
        print('chart %d (t = %g, %s = %g -> x = %.3g, y = %.3g): series vs numerical, 579 generators: max rel diff %.1e (%s), median %.1e'
              % (c, tv, big, Bv, P[0], P[1], d.max(), names[int(d.argmax())], np.median(d)), flush=True)
