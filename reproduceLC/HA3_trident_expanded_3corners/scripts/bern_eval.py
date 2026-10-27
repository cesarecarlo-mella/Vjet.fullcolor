#!/usr/bin/env python3
"""Bernoulli re-expansion, consistency check against HAser (c1) and digit evaluation on the grid
for the fully expanded HA3 series of expand_corner.py.

Bernoulli variables: bp = -log(1-p), bq = -log(1-q) for the corner variables (p, q)
(c1: bx, by;  c2: by, bz;  c3: bx, bz);  p = 1 - e^-bp.  Log[p], Log[q] are kept as they are.
Truncation: total degree in (bp, bq) <= the order of the native series.

usage: python3 bern_eval.py <corner> <dir> <exact.json> [HAser.m]
"""
import json
import math
import pickle
import re
import sys
import time
from fractions import Fraction

import mpmath as mp

corner, DIR, EXJ = sys.argv[1], sys.argv[2], sys.argv[3]
HASER = sys.argv[4] if len(sys.argv) > 4 else None
mp.mp.dps = 40
t0 = time.time()
d = pickle.load(open('%s/HA3exp_%s.pkl' % (DIR, corner), 'rb'))
H, N, VN = d['H'], int(d['order']), d['vars']
H = {tuple(int(e) for e in st): {(int(a), int(b)): c for (a, b), c in cs.items()} for st, cs in H.items()}


def mul(a, b, n):
    out = {}
    for i, x in a.items():
        for j, y in b.items():
            if i + j <= n:
                out[i + j] = out.get(i + j, 0) + x * y
    return out


def bpow(a, n):
    """(1 - e^-b)^a as a Laurent series in b up to b^n (a may be negative)"""
    e = {k: Fraction((-1) ** (k + 1), math.factorial(k)) for k in range(1, n + 8)}   # 1-e^-b = sum e_k b^k
    w = {k - 1: c for k, c in e.items()}                                                # (1-e^-b)/b
    if a >= 0:
        r = {0: Fraction(1)}
        for _ in range(a):
            r = mul(r, w, n + 8)
        return {k + a: c for k, c in r.items() if k + a <= n}
    inv = {0: Fraction(1)}                                                                # (b/(1-e^-b))
    for k in range(1, n + 7):
        inv[k] = -sum(w[j] * inv[k - j] for j in range(1, k + 1))
    r = {0: Fraction(1)}
    for _ in range(-a):
        r = mul(r, inv, n + 6)
    return {k + a: c for k, c in r.items() if k + a <= n}


cache = {}


def bp(a, n):
    if (a, n) not in cache:
        cache[(a, n)] = bpow(a, n)
    return cache[(a, n)]


Bn = {}
for st, cs in H.items():
    acc = {}
    for (a, b), c in cs.items():
        A = bp(a, N - b)
        for i, ca in A.items():
            Bq = bp(b, N - i)
            for j, cb in Bq.items():
                if i + j <= N:
                    acc[(i, j)] = acc.get((i, j), 0) + c * ca * cb
    Bn[st] = {m: c for m, c in acc.items() if c != 0}
print('%s Bernoulli: %d coefficients (%.0fs)' % (corner, sum(len(v) for v in Bn.values()), time.time() - t0), flush=True)
pickle.dump({'corner': corner, 'vars': ['b' + v for v in VN], 'order': N, 'H': Bn}, open('%s/HA3bern_%s.pkl' % (DIR, corner), 'wb'))

BV = ('b' + VN[0], 'b' + VN[1])


def term(st, a, b, c):
    f = []
    for sym, e in ((BV[0], a), (BV[1], b), ('Log[%s]' % VN[0], st[0]), ('Log[%s]' % VN[1], st[1]),
                   ('Pi', st[2]), ('Zeta[3]', st[3]), ('Zeta[5]', st[4]), ('I', st[5])):
        if e:
            f.append(sym if e == 1 else '%s^%d' % (sym, e))
    return '*'.join([str(c)] + f)


name = 'HA3bern["C3F0Z0", 3, "%s"]' % corner
with open('%s/C3F0Z0_HA3_bernoulli_%s.m' % (DIR, corner), 'w') as fo:
    fo.write('(* %s: Bernoulli form of HA3exp["C3F0Z0", 3, "%s"]:\n'
             '   %s = -Log[1-%s], %s = -Log[1-%s]  (i.e. %s = 1 - Exp[-%s], ...); Log[%s], Log[%s] kept.\n'
             '   Truncation: total degree in (%s, %s) <= %d. *)\n'
             % (name, corner, BV[0], VN[0], BV[1], VN[1], VN[0], BV[0], VN[0], VN[1], BV[0], BV[1], N))
    fo.write('%s = (0' % name)
    for st in sorted(Bn):
        for (a, b) in sorted(Bn[st], key=lambda m: (m[0] + m[1], m)):
            fo.write('\n + %s' % term(st, a, b, Bn[st][(a, b)]))
    fo.write(');\n')

# ---------------------------------------------------------------- c1: compare with HAser (from GPLs)
if HASER:
    b = re.sub(r'\(\*.*?\*\)', '', open(HASER).read(), flags=re.S).split('=', 1)[1].strip()[1:].rstrip()[:-2]
    ser = {}
    for tm in b.split('\n + '):
        st = [0] * 6; a = bb = 0; c = Fraction(1)
        for f in tm.strip().split('*'):
            if f[0] == '-' and not f[1].isdigit():
                c = -c; f = f[1:]
            base, _, e = f.partition('^'); e = int(e) if e else 1
            if base == 'x': a = e
            elif base == 'y': bb = e
            elif base == 'Log[x]': st[0] = e
            elif base == 'Log[y]': st[1] = e
            elif base == 'Pi': st[2] = e
            elif base == 'Zeta[3]': st[3] = e
            elif base == 'Zeta[5]': st[4] = e
            elif base == 'I': st[5] = 1
            else: c *= Fraction(f)
        if a + bb <= N:
            k = (tuple(st), (a, bb))
            ser[k] = ser.get(k, 0) + c
    mine = {(st, m): c for st, cs in H.items() for m, c in cs.items() if c != 0}
    ser = {k: v for k, v in ser.items() if v != 0}
    diff = [k for k in set(ser) | set(mine) if ser.get(k, 0) != mine.get(k, 0)]
    print('c1 check vs HAser (GPL series) restricted to total degree <= %d: %d coefficients, %d differ'
          % (N, len(ser), len(diff)), flush=True)

# ---------------------------------------------------------------- numerics on the grid
CON = {}


def const(st):
    if st not in CON:
        CON[st] = mp.pi ** st[2] * mp.zeta(3) ** st[3] * mp.zeta(5) ** st[4] * (mp.mpc(0, 1) ** st[5])
    return CON[st]


def evaluate(S, p, q, lp, lq):
    tot = mp.mpc(0)
    for st, cs in S.items():
        s = mp.mpf(0)
        for (a, b), c in cs.items():
            s += mp.mpf(c.numerator) / c.denominator * p ** a * q ** b
        tot += s * lp ** st[0] * lq ** st[1] * const(st)
    return tot


rows = json.load(open(EXJ))['rows']
out = []
for r in rows:
    x = mp.mpf(Fraction(r['x']).limit_denominator(1000).numerator) / Fraction(r['x']).limit_denominator(1000).denominator
    y = mp.mpf(Fraction(r['y']).limit_denominator(1000).numerator) / Fraction(r['y']).limit_denominator(1000).denominator
    z = 1 - x - y
    p, q = {'c1': (x, y), 'c2': (y, z), 'c3': (x, z)}[corner]
    lp, lq = mp.log(p), mp.log(q)
    nat = evaluate(H, p, q, lp, lq)
    ber = evaluate(Bn, -mp.log(1 - p), -mp.log(1 - q), lp, lq)
    out.append({'x': r['x'], 'y': r['y'], 'exact': r['exact'][2],
                'native': [mp.nstr(nat.real, 25), mp.nstr(nat.imag, 25)],
                'bernoulli': [mp.nstr(ber.real, 25), mp.nstr(ber.imag, 25)]})
json.dump({'corner': corner, 'order': N, 'rows': out}, open('%s/eval_%s.json' % (DIR, corner), 'w'))
print('%s evaluated on %d points (%.0fs)' % (corner, len(out), time.time() - t0))
