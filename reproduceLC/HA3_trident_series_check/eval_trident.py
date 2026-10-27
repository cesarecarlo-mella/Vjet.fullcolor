#!/usr/bin/env python3
"""HA3 (C3F0Z0) with every Trident[k] replaced by its series around corner cK, truncated at
order N in both corner variables (u^a v^b with a, b <= N), evaluated on the 406-point grid and
compared with the exact HA3 (our GPL form, gpl_eval.py, ~1e-16).

usage: python3 eval_trident.py <corner> <N> <exact.json> <out.json>
"""
import json
import sys
import time
from fractions import Fraction

import numpy as np
import mpmath as mp

from wparse import parse, peval

LD = np.longdouble
mp.mp.dps = 40
corner, N, EXJ, OUT = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
CHUNK = 4


def ld(q):
    """Fraction or int -> longdouble with ~25 significant digits"""
    q = Fraction(q)
    if q == 0:
        return LD(0)
    s = mp.nstr(mp.mpf(q.numerator) / q.denominator, 25)
    return LD(s)


t0 = time.time()
d = np.load('parsed_%s.npz' % corner)
keys = [tuple(k) for k in d['keys']]
ex = d['ex'].astype(np.int64)
coef = np.load('parsed_%s_coef.npy' % corner)
keep = (ex[:, 1] <= N) & (ex[:, 2] <= N)
ex = ex[keep]; coef = coef[keep]
# fold the constants  Pi^p Zeta3^z3 Zeta5^z5  into the coefficient
PI, Z3, Z5 = (LD(mp.nstr(c, 25)) for c in (mp.pi, mp.zeta(3), mp.zeta(5)))
coef = coef * PI ** ex[:, 5].astype(LD) * Z3 ** ex[:, 6].astype(LD) * Z5 ** ex[:, 7].astype(LD)
print('coefficients ready: %d terms (%.0fs)' % (len(ex), time.time() - t0), flush=True)

fid = ex[:, 0]
nf = len(keys)
starts = np.flatnonzero(np.r_[True, fid[1:] != fid[:-1]])
present = fid[starts]
A_, B_, LA, LB = ex[:, 1], ex[:, 2], ex[:, 3], ex[:, 4]

# prefactors
wl = parse('C3F0Z0_HA3.wl')
pos = {k: i for i, k in enumerate(keys)}
for part in wl.values():
    for k, _, _ in part:
        assert k in pos, k

rows = json.load(open(EXJ))['rows']
if len(sys.argv) > 5: rows = [rows[int(i)] for i in sys.argv[5].split(',')]
out = []
VAR = {'c1': lambda x, y: (x, y), 'c2': lambda x, y: (y, 1 - x - y), 'c3': lambda x, y: (x, 1 - x - y)}
for c0 in range(0, len(rows), CHUNK):
    chunk = rows[c0:c0 + CHUNK]
    P = len(chunk)
    uv = []
    for r in chunk:
        xq = Fraction(r['x']).limit_denominator(1000); yq = Fraction(r['y']).limit_denominator(1000)
        r['xq'], r['yq'] = xq, yq
        uv.append(VAR[corner](xq, yq))
    U = np.array([ld(a) for a, _ in uv]); V = np.array([ld(b) for _, b in uv])
    LU = np.array([LD(mp.nstr(mp.log(mp.mpf(a.numerator) / a.denominator), 25)) for a, _ in uv])
    LV = np.array([LD(mp.nstr(mp.log(mp.mpf(b.numerator) / b.denominator), 25)) for _, b in uv])
    PU = U[None, :] ** np.arange(N + 1, dtype=LD)[:, None]
    PV = V[None, :] ** np.arange(N + 1, dtype=LD)[:, None]
    QU = LU[None, :] ** np.arange(7, dtype=LD)[:, None]
    QV = LV[None, :] ** np.arange(7, dtype=LD)[:, None]
    val = coef[:, None] * PU[A_] * PV[B_] * QU[LA] * QV[LB]
    T = np.zeros((nf, P), dtype=LD)
    T[present] = np.add.reduceat(val, starts, axis=0)
    del val
    for p, r in enumerate(chunk):
        t, u = r['yq'], 1 - r['xq'] - r['yq']
        res = {}
        for part, lst in wl.items():
            s = LD(0); s_abs = LD(0)
            for k, nm, dn in lst:
                R = ld(peval(nm, t, u) / peval(dn, t, u))
                term = R * T[pos[k], p]
                s += term; s_abs += abs(term)
            res[part] = (s, s_abs)
        re_, im_ = res['re'][0], mp.pi * mp.mpf(str(res['ipi'][0]))
        e_re, e_im = r['exact'][2]
        out.append({'x': r['x'], 'y': r['y'], 'series': [str(re_), mp.nstr(im_, 25)],
                    'exact': [e_re, e_im],
                    'cancel': [float(res['re'][1] / max(abs(res['re'][0]), LD(1e-300))),
                               float(res['ipi'][1] / max(abs(res['ipi'][0]), LD(1e-300)))],
                    # high-precision copies of the difference
                    'diff': [mp.nstr(mp.mpf(str(re_)) - mp.mpf(e_re), 8), mp.nstr(im_ - mp.mpf(e_im), 8)]})
    print('%d/%d points (%.0fs)' % (c0 + P, len(rows), time.time() - t0), flush=True)

json.dump({'corner': corner, 'N': N, 'rows': out}, open(OUT, 'w'))
print('done', time.time() - t0)
