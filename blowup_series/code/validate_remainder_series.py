#!/usr/bin/env python3
"""the remainder series (assemble_remainders.py output: <root>/<coefficient>/<layer>/series<c>.m) against the numerical
evaluator, at points near each corner.  The series is evaluated in its native variables; the numerical values come
with their error estimate.  Also printed: the size of the highest-order terms (a rough truncation estimate).
usage: validate_remainder_series.py <root> [corners c1,c2,c3] [distance 0.05]"""
import math, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', '..', 'iter int ', 'generators_deq'))
import mparse as mp            # noqa: E402
import remainder_eval as R     # noqa: E402

root = sys.argv[1]
corners = (sys.argv[2] if len(sys.argv) > 2 else 'c1,c2,c3').split(',')
d0 = float(sys.argv[3]) if len(sys.argv) > 3 else 0.05
VARS = {'c1': ('x', 'y'), 'c2': ('y', 'z'), 'c3': ('x', 'z')}
CONST = {'Pi': math.pi, 'Zeta[3]': 1.2020569031595942854, 'Zeta[5]': 1.0369277551433699263, 'I': 1j}
LNAME = lambda key: '_'.join('%s%s' % (n, v) for n, v in zip(('Nc', 'Nf', 'NfZ'), re.findall(r'-?\d+', key)))


def point(c, s, b):
    """(x, y) from the native small/big variables of the corner"""
    return {'c1': (s, b), 'c2': (1 - s - b, s), 'c3': (s, 1 - s - b)}[c]


for c in corners:
    sm, bg = VARS[c]
    pts = [(d0, d0), (d0 / 2, d0), (d0, d0 / 2)]
    P = [point(c, s, b) for s, b in pts]
    Rv, info, err = R.evaluate_remainders(P, 1e-12, error=True)
    worst = (0, None)
    nlay = 0
    for j, (name, key) in enumerate(R.KEYS):
        fn = os.path.join(root, name, LNAME(key), 'series%s.m' % c)
        if not os.path.isfile(fn):
            continue
        poly = mp.parse_file(fn)
        N = int(re.search(r'rectangle of order (\d+)', open(fn).readline()).group(1))
        nlay += 1
        for ip, (s, b) in enumerate(pts):
            val = {sm: s, bg: b, 'Log[%s]' % sm: math.log(s), 'Log[%s]' % bg: math.log(b)}
            val.update(CONST)
            tot = top = 0j
            for mono, q in poly.items():
                v = float(q)
                deg = 0
                for sym, e2 in mono:
                    v = v * val[sym] ** (e2 / 2)
                    if sym in (sm, bg):
                        deg = max(deg, e2 / 2)
                tot += v
                if deg >= N:
                    top += v
            ref = Rv[ip, j]
            rel = abs(tot - ref) / abs(ref)
            bound = (err[ip, j] + abs(top)) / abs(ref)
            if rel / max(bound, 1e-300) > worst[0]:
                worst = (rel / max(bound, 1e-300), (name, key, pts[ip], rel, err[ip, j] / abs(ref), abs(top) / abs(ref)))
    w = worst[1]
    print('%s: %d layers x %d points near the corner (%s, %s ~ %g): worst |series - numerical| / (numerical error + order-N terms) = %.2f'
          % (c, nlay, len(pts), sm, bg, d0, worst[0]))
    if w:
        print('    at %s %s, (%s, %s) = %s: relative difference %.1e, numerical error %.1e, order-N terms %.1e' % (w[0], w[1], sm, bg, w[2], w[3], w[4], w[5]))
