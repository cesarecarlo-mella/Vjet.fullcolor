#!/usr/bin/env python3
"""
trident_rotation.py -- series of the 1318 trident functions in blow-up
variables, from the family series of run_all_families.py.

At an expansion point c (1: x,y -> 0, chart x = t y^2;  2: y,z -> 0, z = t y^2;
3: x,z -> 0, x = t z^2) the crossed family element BASISLC[F, crossing, i] is
master i of family F, at a permutation of (x, y, z).  In blow-up variables it
is F's series in one of the six charts (1-3 standard, 4-6 mirror), with the
SAME (t, BIG) -- only the name of BIG changes -- so it is read verbatim from
families/<F>/<dist>/rec_c<chart>_w<w>.pkl.  With the 1318 rows of
trident/rotation_rows.json (select_trident_rows.py):

    Trident = m_sel^-1 . (BASISLC_sel / phase)

computed modulo primes < 2^22 (float64 BLAS, exact), then CRT + rational
reconstruction (blowup_array.reconstruct), verified on 2 extra primes; primes
are added until every coefficient verifies.  The 7 algebraic tridents are then
divided by I (real, phase-normalised convention of the published series and of
the amplitude coefficients), and the result must be real.

Output (per point c and weight w, in --out):
    trident_c<c>_w<w>.pkl   {(tag, i, it, iy): mpq}, i = trident index - 1
    trident_c<c>_w<w>.m     list of the 1318 series in t, BIG, Log[t], Log[BIG]
    phys_c<c>_w<w>_ord<N>.m with --phys N (checks only)

    python3 trident_rotation.py --dist dist_22_66 --out ../trident/dist_22_66 -n 64
"""
import argparse
import json
import multiprocessing as mproc
import os
import sys
import time
from types import SimpleNamespace

import numpy as np
from gmpy2 import mpq, mpz

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import blowup_array as ba
import mparse as mp
import phys
import pipe_common as pc

ap = argparse.ArgumentParser()
ap.add_argument('--families', default=os.path.join(HERE, '..', 'families'))
ap.add_argument('--dist', required=True, help='name of the run directory inside each family, e.g. dist_22_66')
ap.add_argument('--rows', default=os.path.join(HERE, '..', 'trident', 'rotation_rows.json'))
ap.add_argument('--out', required=True)
ap.add_argument('-c', '--corners', default='1,2,3')
ap.add_argument('-w', '--weights', default='0,1,2,3,4,5,6')
ap.add_argument('-n', '--jobs', type=int, default=os.cpu_count())
ap.add_argument('--nprimes', type=int, default=12, help='initial number of primes (more are added if needed)')
ap.add_argument('--phys', type=int, default=None, help='also write phys_c<c>_w<w>_ord<N>.m')
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
NTRI = 1318

# ---------------------------------------------------------------- crossings
# BASISLC[F, crossing, i](x, y, z) = F_i(X, Y, Z) with (X, Y, Z) = these
# permutations of (x, y, z); x = s12, y = s13, z = s23, and xij swaps legs i, j
PERM = {"": "xyz", "x12": "xzy", "x13": "zyx", "x23": "yxz", "x123": "zxy", "x132": "yzx"}
POINT = {1: ('x', 'y'), 2: ('z', 'y'), 3: ('x', 'z')}          # (SMALL, BIG) of the chart
CHART = {('X', 'Y'): 1, ('Z', 'Y'): 2, ('X', 'Z'): 3,          # family charts, standard
         ('Y', 'X'): 4, ('Y', 'Z'): 5, ('Z', 'X'): 6}          # and mirror


def chart_of(c, crossing):
    """family chart that gives BASISLC[F, crossing, i] at point c."""
    perm = PERM[crossing]
    name = lambda v: 'XYZ'[perm.index(v)]
    sm, bg = POINT[c]
    return CHART[(name(sm), name(bg))]


# ---------------------------------------------------------------- inputs
rows = json.load(open(a.rows))['rows']
assert len(rows) == NTRI
orders = set()
for F in sorted({r['family'] for r in rows}):
    with open(os.path.join(a.families, F, a.dist, '.order')) as f:
        orders.add(tuple(int(x) for x in f.read().split()))
if len(orders) != 1:
    sys.exit("families were run at different orders: %s" % orders)
ORDt, ORDy = orders.pop()
Nt, Ny = ORDt + 2, ORDy + 2
with open(os.path.join(a.out, '.order'), 'w') as f:
    f.write("%d %d\n" % (ORDt, ORDy))
print("trident rotation: order %d/%d, %d rows, points %s, weights %s" % (ORDt, ORDy, NTRI, a.corners, a.weights), flush=True)


def phase_div(series, phase):
    """series / phase, phase in {1, -1, I, -I};  series = {(tag, it, iy): q}."""
    if phase == '1':
        return series
    out = {}
    for (tag, it, iy), q in series.items():
        pt, py, la, lb, ck = tag
        if phase == '-1':
            ck2, f = ck, -1
        else:                                     # 1/I = -I, 1/(-I) = I
            ck2, s = ba.key_mul(ck, ba.I_KEY)
            f = -s if phase == 'I' else s
        out[((pt, py, la, lb, ck2), it, iy)] = f * q
    return out


def load_selected(c, w):
    """the 1318 selected crossed family elements at point c, weight w (divided by their phase)."""
    need = {}
    for k, r in enumerate(rows):
        need.setdefault((r['family'], chart_of(c, r['crossing'])), []).append(k)
    sel = [None] * NTRI
    for (F, ch), ks in sorted(need.items()):
        if w == 0:
            v = mp.parse_file(os.path.join(a.families, F, 'inputs', 'sol0.m'))
            comps = {}
            for j, poly in enumerate(v):
                for key, q in poly.items():
                    comps.setdefault(j, {})[((0, 0, 0, 0, key), 1, 1)] = q
        else:
            slots = pc.load_pickle(os.path.join(a.families, F, a.dist, 'rec_c%d_w%d.pkl' % (ch, w)))
            want = {rows[k]['index'] - 1 for k in ks}
            comps = {}
            for (tag, j, it, iy), q in slots.items():
                if j in want:
                    comps.setdefault(j, {})[(tag, it, iy)] = q
            del slots
        for k in ks:
            sel[k] = phase_div(comps.get(rows[k]['index'] - 1, {}), rows[k]['phase'])
    return sel


def minv_mod(p):
    """m_sel^-1 mod p (balanced, float64), or None if p is bad for m_sel."""
    M = np.zeros((NTRI, NTRI), dtype=np.int64)
    try:
        for k, r in enumerate(rows):
            for j, v in r['entries']:
                n, _, d = v.partition('/')
                M[k, j - 1] = ba.rmod(mpq(int(n), int(d or 1)), p)
    except ZeroDivisionError:
        return None
    try:
        from flint import nmod_mat
        Mi = nmod_mat(M.tolist(), p).inv()
        Mi = np.array([[int(Mi[i, j]) for j in range(NTRI)] for i in range(NTRI)], dtype=np.int64)
    except ImportError:                           # plain Gauss-Jordan mod p
        A = np.concatenate([M % p, np.eye(NTRI, dtype=np.int64)], axis=1)
        for col in range(NTRI):
            piv = next((r for r in range(col, NTRI) if A[r, col]), None)
            if piv is None:
                return None
            A[[col, piv]] = A[[piv, col]]
            A[col] = A[col] * pow(int(A[col, col]), -1, p) % p
            f = A[:, col].copy(); f[col] = 0
            A -= np.outer(f, A[col]) % p
            A %= p
        Mi = A[:, NTRI:]
    except (ZeroDivisionError, ValueError):
        return None
    Mi = np.where(Mi > p // 2, Mi - p, Mi)
    return Mi.astype(np.float64)


_G = {}          # fork-shared: entries of the selected series by tag, inverses per prime


def _minv(k):
    return k, minv_mod(pc.PRIME_TABLE[k])


def tag_series(tag):
    """reconstruct the trident series of ONE tag over the primes in _G['ks']."""
    ents = _G['bytag'][tag]
    rr = np.array([e[0] for e in ents]); ii = np.array([e[1] for e in ents]); yy = np.array([e[2] for e in ents])
    results, primes = [], []
    for k in _G['ks']:
        p = pc.PRIME_TABLE[k]
        B = np.zeros((NTRI, Nt, Ny), np.int64)
        try:
            B[rr, ii, yy] = [ba.rmod(e[3], p) for e in ents]
        except ZeroDivisionError:
            continue                                  # p divides a denominator of the families
        B = B.reshape(NTRI, Nt * Ny)
        B = np.where(B > p // 2, B - p, B).astype(np.float64)
        T = _G['minv'][k] @ B                         # |terms| < 1318 * 2^42 < 2^53: exact
        results.append({tag: np.remainder(T, p).astype(np.int64).reshape(NTRI, Nt, Ny)})
        primes.append(p)
    out, bad, _ = ba.reconstruct(results, primes, jobs=1)
    return tag, out, bad


# The 7 algebraic tridents (prefactor Sqrt[u v]) contain an odd number of the
# odd letters l7, l8 -- logs of a pure phase -- so m^-1 . BASISLC gives them as
# I * (real series).  The published series and the amplitude coefficients
# (trident_coefficients/*.wl) use the real, phase-normalised functions:
# Trident = (m^-1 . BASISLC) / I for these 7.  Checked exactly against the
# published series at all three points.
PHASE_NORMALISED = {1, 2, 3, 15, 16, 17, 1283}


def normalise(slots):
    out = {}
    for (tag, i, it, iy), q in slots.items():
        if i + 1 in PHASE_NORMALISED:
            pt, py, la, lb, ck = tag
            ck2, s = ba.key_mul(ck, ba.I_KEY)          # q * ck / I = -q * ck * I
            tag, q = (pt, py, la, lb, ck2), -s * q
        out[(tag, i, it, iy)] = q
    return out


def write_series(c, slots, path):
    C = SimpleNamespace(c=c, ncomp=NTRI, tv='t', yv=ba.BIG_VAR[c], ltv='Log[t]', lyv='Log[%s]' % ba.BIG_VAR[c])
    pc.write_exact(C, slots, path)
    return C


T0 = time.time()
allok = True
_G['minv'], next_k = {}, 0


def add_primes(n):
    """make sure inverses exist for n more usable primes."""
    global next_k
    want = len(_G['minv']) + n
    while len(_G['minv']) < want:
        todo = list(range(next_k, next_k + want - len(_G['minv'])))
        next_k += len(todo)
        with mproc.get_context('fork').Pool(min(a.jobs, len(todo))) as pool:
            for k, Mi in pool.imap_unordered(_minv, todo):
                if Mi is None:
                    print("  prime k=%d skipped (m_sel singular or a denominator mod p)" % k)
                else:
                    _G['minv'][k] = Mi


add_primes(a.nprimes + ba.NV)
for c in [int(s) for s in a.corners.split(',')]:
    for w in [int(s) for s in a.weights.split(',')]:
        t0 = time.time()
        bytag = {}
        for r, ser in enumerate(load_selected(c, w)):
            for (tag, it, iy), q in ser.items():
                bytag.setdefault(tag, []).append((r, it, iy, q))
        _G['bytag'] = bytag
        tl = time.time() - t0
        out, todo, nprimes = {}, sorted(bytag), len(_G['minv'])
        while todo:
            _G['ks'] = sorted(_G['minv'])[:nprimes]
            failed = []
            with mproc.get_context('fork').Pool(a.jobs) as pool:
                for tag, o, bad in pool.imap_unordered(tag_series, todo):
                    out.update(o)
                    if bad:
                        failed.append(tag)
            if failed:
                print("  c%d w%d: %d tags unverified with %d primes -> adding primes"
                      % (c, w, len(failed), nprimes - ba.NV), flush=True)
                nprimes += max(4, a.jobs // 4)
                if nprimes > len(_G['minv']):
                    add_primes(nprimes - len(_G['minv']))
            todo = failed
        ks = _G['ks']
        out = normalise(out)
        imag = sum(1 for (tag, i, it, iy) in out if any(at == 'I' for at, _ in tag[4]))
        if imag:
            allok = False
        pc.save_pickle(out, os.path.join(a.out, 'trident_c%d_w%d.pkl' % (c, w)))
        C = write_series(c, out, os.path.join(a.out, 'trident_c%d_w%d.m' % (c, w)))
        if a.phys is not None:
            pd = phys.slots_to_phys(SimpleNamespace(c=c), out, ORDt, a.phys)
            phys.write_phys(C, pd, os.path.join(a.out, 'phys_c%d_w%d_ord%d.m' % (c, w, a.phys)))
        print("c%d w%d: %d coefficients, %d primes (+%d verification)%s  (load %.1fs, total %.1fs)"
              % (c, w, len(out), len(ks) - ba.NV, ba.NV,
                 "  !! %d IMAGINARY coefficients" % imag if imag else "", tl, time.time() - t0), flush=True)
        _G.pop('bytag', None)
print("done in %.1f min%s" % ((time.time() - T0) / 60, "" if allok else "  -- WITH IMAGINARY PARTS (check phases/crossings)"))
sys.exit(0 if allok else 1)
