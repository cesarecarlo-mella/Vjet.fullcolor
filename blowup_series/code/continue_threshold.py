#!/usr/bin/env python3
"""
continue_threshold.py -- series at the threshold corner of each continued region (the old CC1): the decay corner's
series, continued.  No differential equation is solved (see transport/threshold.py for the regions and the convention).

    rec_c<c>_w<w>.pkl  ->  exact_cont_c<c>_w<w>.m        c = 1, 4: z > 1    c = 2, 5: x > 1    c = 3, 6: y > 1

In the decay chart SMALL = t BIG^2, with SMALL and BIG both negative in the continued region:
    Log[t] -> Log[t] - I Pi ,  Log[BIG] -> Log[BIG] + I Pi ,  t^e -> E^(-I Pi e) t^e ,  BIG^m -> E^(I Pi m) BIG^m
where in the output t and BIG stand for |t| and |BIG|, i.e. SMALL = -t BIG^2 and BIG_phys = -BIG.  This is the CC
convention Log b -> log|b| + I Pi for the two small variables, Sqrt[x y z] -> -Sqrt|x y z|.

    python3 continue_threshold.py --workdir families/H/dist_22_66 [-c 1,2,3,4,5,6] [-w 1,2,3,4,5,6] [--check X,Y]

--check X,Y compares, in physical variables and in the window SMALL^X BIG^Y, with the continuation of the physical
series done directly in (SMALL, BIG) (Log -> log|.| + I Pi, half-integer powers -> I^n |.|^(n/2)).
"""
import argparse
import math
import os
import pickle
from types import SimpleNamespace

import charts
import mparse as mp
import phys
from gmpy2 import mpq

ap = argparse.ArgumentParser()
ap.add_argument('--workdir', required=True)
ap.add_argument('-c', default='1,2,3,4,5,6')
ap.add_argument('-w', default='1,2,3,4,5,6')
ap.add_argument('--check', default=None)
a = ap.parse_args()

IPI = {(('I', 2), ('Pi', 2)): mpq(1)}
MIPI = {(('I', 2), ('Pi', 2)): mpq(-1)}


def ipow_unit(n):
    """I^n as a Poly"""
    n %= 4
    return {(): mpq(1)} if n == 0 else {(('I', 2),): mpq(1)} if n == 1 else {(): mpq(-1)} if n == 2 else {(('I', 2),): mpq(-1)}


def binom_shift(atom, n, shift):
    """(atom + shift)^n with shift a Poly (I Pi or -I Pi)"""
    out, sp = {}, {(): mpq(1)}
    pw = [sp]
    for _ in range(n):
        pw.append(mp.pmul(pw[-1], shift))
    for j in range(n + 1):
        mono = {((atom, 2 * j),) if j else (): mpq(math.comb(n, j))}
        out = mp.padd(out, mp.pmul(mono, pw[n - j]))
    return out


def continue_slots(slots, c, ncomp):
    tv, bv = 't', charts.big(c)
    lt, lb = 'Log[t]', 'Log[%s]' % bv
    comps = [dict() for _ in range(ncomp)]
    cache = {}
    for (tag, j, it, iy), q in slots.items():
        pt, py, la, lbb, ck = tag
        e2t, e2y = 2 * (it - 1) + pt, 2 * (iy - 1) + py
        k = (e2t, e2y, la, lbb)
        if k not in cache:
            mono = {tuple(sorted(([(tv, e2t)] if e2t else []) + ([(bv, e2y)] if e2y else []))): mpq(1)}
            p = mp.pmul(mono, ipow_unit(e2y - e2t))                 # E^(-I Pi e2t/2) E^(I Pi e2y/2)
            p = mp.pmul(p, binom_shift(lt, la, MIPI))
            cache[k] = mp.pmul(p, binom_shift(lb, lbb, IPI))
        comps[j] = mp.padd(comps[j], mp.pmul({ck: q}, cache[k]))
    return comps


def to_phys(poly, c, X, Y):
    """continued chart Poly (t, BIG) -> physical (SMALL, BIG) Poly in the window: t^e BIG^m -> SMALL^e BIG^(m-2e),
    Log[t] -> Log[SMALL] - 2 Log[BIG] (names SMALL, BIG as in charts.small_big)"""
    sm, bg = charts.small_big(c)
    out = {}
    for key, q in poly.items():
        d = dict(key)
        e2t, e2y = d.pop('t', 0), d.pop(bg, 0)
        la, lbb = d.pop('Log[t]', 0) // 2, d.pop('Log[%s]' % bg, 0) // 2
        ey = e2y - 2 * e2t
        if e2t > 2 * X or ey > 2 * Y:
            continue
        rest = list(d.items()) + ([(sm, e2t)] if e2t else []) + ([(bg, ey)] if ey else [])
        for kk in range(la + 1):
            key2 = rest + ([('Log[%s]' % sm, 2 * (la - kk))] if la - kk else []) + ([('Log[%s]' % bg, 2 * (lbb + kk))] if lbb + kk else [])
            out = mp.padd(out, {tuple(sorted(key2)): q * math.comb(la, kk) * (-2) ** kk})
    return out


def direct(pd, c, ncomp):
    """physical slots {(j, ex2, ey2, lx, ly, ck): q} continued directly in (SMALL, BIG)"""
    sm, bg = charts.small_big(c)
    comps = [dict() for _ in range(ncomp)]
    for (j, ex2, ey2, lx, ly, ck), q in pd.items():
        mono = {tuple(sorted(([(sm, ex2)] if ex2 else []) + ([(bg, ey2)] if ey2 else []))): mpq(1)}
        p = mp.pmul(mp.pmul(mono, ipow_unit(ex2 + ey2)), binom_shift('Log[%s]' % sm, lx, IPI))
        p = mp.pmul(p, binom_shift('Log[%s]' % bg, ly, IPI))
        comps[j] = mp.padd(comps[j], mp.pmul({ck: q}, p))
    return comps


for c in [int(x) for x in a.c.split(',')]:
    for w in [int(x) for x in a.w.split(',')]:
        src = os.path.join(a.workdir, 'rec_c%d_w%d.pkl' % (c, w))
        slots = pickle.load(open(src, 'rb'))
        ncomp = 1 + max(j for (_, j, _, _) in slots)
        exact = mp.parse_file(os.path.join(a.workdir, 'exact_c%d_w%d.m' % (c, w)))
        ncomp = max(ncomp, len(exact))
        comps = continue_slots(slots, c, ncomp)
        dst = os.path.join(a.workdir, 'exact_cont_c%d_w%d.m' % (c, w))
        with open(dst + '.tmp', 'w') as f:
            f.write("{" + ",\n ".join(mp.fmt_poly(p) for p in comps) + "}\n")
        os.replace(dst + '.tmp', dst)
        msg = 'c%d w%d: wrote %s' % (c, w, os.path.basename(dst))
        if a.check:
            X, Y = [int(v) for v in a.check.split(',')]
            ref = direct(phys.slots_to_phys(SimpleNamespace(c=c), slots, X, Y), c, ncomp)
            bad = sum(1 for p, r in zip(comps, ref) if mp.padd(to_phys(p, c, X, Y), r, -1))
            msg += '; check in the window %d,%d: %d of %d components differ' % (X, Y, bad, ncomp)
        print(msg, flush=True)
