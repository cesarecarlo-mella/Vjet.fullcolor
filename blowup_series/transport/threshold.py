#!/usr/bin/env python3
"""
threshold.py -- constants of the threshold corner of each continued region (the old CC1).

Region A > 1 (the other two variables negative) touches the decay triangle at the decay corner where the
other two variables vanish, and its threshold corner is that decay corner, continued:

    z > 1  <->  c1 (x, y -> 0)      x > 1  <->  c2 (y, z -> 0)      y > 1  <->  c3 (x, z -> 0)

Continuation (the CC convention): Log b -> log|b| + I Pi for both small variables, Sqrt[x y z] -> -Sqrt|x y z|.
In the decay chart SMALL = t BIG^2 (SMALL and BIG both negative there) this is

    Log[t] -> Log|t| - I Pi ,   Log[BIG] -> Log|BIG| + I Pi ,   t^e -> E^(-I Pi e) |t|^e ,   BIG^m -> E^(I Pi m) |BIG|^m

and the threshold constant (every log of the new region -> 0) is the x^0 y^0 coefficient of the decay series,
a polynomial in Log[t], Log[BIG], at Log[t] -> -I Pi, Log[BIG] -> I Pi.

  low     the x^0 y^0 coefficient of every master from a decay run (rec_c<c>_w<w>.pkl, any order)
          -> families/<F>/inputs/low_c<c>_w<w>.m        (c = 1..3, w = 1..6)
  consts  -> families/<F>/inputs/bdy_R<A>_th_w<w>.m     (w = 0..6; w = 0 is sol0)
          --check-H DIR compares H, region z, with the old bdy_CC1_w<w>.m of transport.boundaries/H/dist.

    python3 threshold.py low --runs <dir with <F>/rec_c*_w*.pkl>
    python3 threshold.py consts [--check-H DIR]
"""
import argparse
import os
import pickle
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'code'))
import charts
import mparse as mp
from gmpy2 import mpq

REGION = {1: 'Z', 2: 'X', 3: 'Y'}
FAMILIES = 'A B1 B2 C D E1 E2 F1 F2 G1 G2 G3 G4 H I1 I2 I3 L3 M3 N3 P3 R2'.split()
WMAX = 6
MINUS_IPI = {(('I', 2), ('Pi', 2)): mpq(-1)}
PLUS_IPI = {(('I', 2), ('Pi', 2)): mpq(1)}


def write_list(path, comps, header):
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        f.write('(* %s *)\n{' % header + ',\n '.join(mp.fmt_poly(c) for c in comps) + '}\n')
    os.replace(tmp, path)


def ncomp(F, fam):
    return len(mp.parse_file(os.path.join(fam, F, 'inputs', 'sol0.m')))


def low(a):
    for F in FAMILIES:
        n = ncomp(F, a.families)
        for c in (1, 2, 3):
            lt, lb = 'Log[t]', 'Log[%s]' % charts.big(c)
            for w in range(1, WMAX + 1):
                d = pickle.load(open(os.path.join(a.runs, F, 'rec_c%d_w%d.pkl' % (c, w)), 'rb'))
                comps = [dict() for _ in range(n)]
                for (tag, j, it, iy), q in d.items():
                    pt, py, la, lbb, ck = tag
                    if it == 1 and iy == 1 and pt == 0 and py == 0:
                        key = list(ck) + ([(lt, 2 * la)] if la else []) + ([(lb, 2 * lbb)] if lbb else [])
                        comps[j][tuple(sorted(key))] = q
                write_list(os.path.join(a.families, F, 'inputs', 'low_c%d_w%d.m' % (c, w)), comps,
                           'family %s, decay corner c%d, weight %d: x^0 y^0 coefficient in Log[t], %s' % (F, c, w, lb))
        print(F, 'low_c{1,2,3}_w{1..%d}.m' % WMAX, flush=True)


def ipow(p, n):
    out = {(): mpq(1)}
    for _ in range(n):
        out = mp.pmul(out, p)
    return out


def threshold_constant(poly, big):
    lt, lb = 'Log[t]', 'Log[%s]' % big
    out = {}
    for key, q in poly.items():
        d = dict(key)
        la, lbb = d.pop(lt, 0) // 2, d.pop(lb, 0) // 2
        assert all(at in ('I', 'Pi', 'Zeta[3]', 'Zeta[5]') for at in d), key
        term = mp.pmul(mp.pmul({tuple(sorted(d.items())): q}, ipow(MINUS_IPI, la)), ipow(PLUS_IPI, lbb))
        out = mp.padd(out, term)
    return out


def consts(a):
    for F in FAMILIES:
        inp = os.path.join(a.families, F, 'inputs')
        sol0 = open(os.path.join(inp, 'sol0.m')).read()
        for c in (1, 2, 3):
            A = REGION[c]
            open(os.path.join(inp, 'bdy_R%s_th_w0.m' % A), 'w').write(sol0)
            for w in range(1, WMAX + 1):
                lowv = mp.parse_file(os.path.join(inp, 'low_c%d_w%d.m' % (c, w)))
                write_list(os.path.join(inp, 'bdy_R%s_th_w%d.m' % (A, w)),
                           [threshold_constant(p, charts.big(c)) for p in lowv],
                           'family %s, region %s > 1, threshold corner (= c%d continued), weight %d' % (F, A.lower(), c, w))
        print(F, 'bdy_R{Z,X,Y}_th_w{0..%d}.m' % WMAX, flush=True)
    if a.check_H:
        inp = os.path.join(a.families, 'H', 'inputs')
        for w in range(0, WMAX + 1):
            mine = mp.parse_file(os.path.join(inp, 'bdy_RZ_th_w%d.m' % w))
            ref = mp.parse_file(os.path.join(a.check_H, 'bdy_CC1_w%d.m' % w))
            bad = sum(1 for x, y in zip(mine, ref) if mp.padd(x, y, -1))
            print('H region z, weight %d: %d of %d components differ from bdy_CC1_w%d.m' % (w, bad, len(ref), w))


ap = argparse.ArgumentParser()
ap.add_argument('step', choices=['low', 'consts'])
ap.add_argument('--runs')
ap.add_argument('--families', default=os.path.join(HERE, '..', 'families'))
ap.add_argument('--check-H', dest='check_H')
a = ap.parse_args()
low(a) if a.step == 'low' else consts(a)
