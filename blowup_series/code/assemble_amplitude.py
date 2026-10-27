#!/usr/bin/env python3
"""
assemble_amplitude.py -- series of helicity coefficients
    A = sum_k R_k(t, u) Trident[i_k, w_k]  +  I Pi sum_k R'_k(t, u) Trident[i_k, w_k]
(t = y, u = z, x = 1 - t - u; files in trident_coefficients/) at the three points.

The rational coefficients have spurious poles at the points (x+y at x,y -> 0, ...), so they are
expanded in the SAME blow-up chart as the trident series (trident_rotation.py) and multiplied
there; only the result goes back to physical variables:

    point 1: x = T y^2, BIG = y      point 2: z = T y^2, BIG = y      point 3: x = T z^2, BIG = z

  * coefficients: R = T^a B^b Ntil(T,B)/Dtil(T,B) with Dtil(0,0) != 0 (asserted: the chart
    resolves every denominator), expanded modulo primes < 2^22;
  * contraction: for every trident tag, one BLAS product over all chart shifts of R; the trident
    residues are cached per prime (--cache) and reused for every coefficient file;
  * CRT + rational reconstruction, verified on 2 extra primes (more primes are added if needed);
  * check: every term with a power of a variable below what the PHYSICAL denominators allow must
    have cancelled (spurious poles);
  * output: the series in the physical variables of the point, to order N in both, and in the
    Bernoulli variables v -> 1 - Exp[-s_v] (sx, sy, sz; LogX, LogY, LogZ placeholders as in
    05_bernoulli.wl), as plain Mathematica expressions:

        <out>/<name>_c<c>_xyz_ord<N>.m        <out>/<name>_c<c>_ber_ord<N>.m

A final order N needs the tridents to (N + q, 3N + p), q and p the largest poles of R in T and B.

    python3 assemble_amplitude.py ../../trident_coefficients/C3F0Z0_HA3.wl \
        --tridents ../trident/dist_22_66 --out ../../amplitude_series -n 32
"""
import argparse
import math
import multiprocessing as mproc
import os
import re
import sys
import time
from types import SimpleNamespace

import numpy as np
import sympy as sp
from gmpy2 import mpq, mpz

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import blowup_array as ba
import mparse as mp
import phys
import pipe_common as pc

ap = argparse.ArgumentParser()
ap.add_argument('files', nargs='+', help='coefficient files sum R(t,u) Trident[i, w]')
ap.add_argument('--tridents', required=True, help='output dir of trident_rotation.py')
ap.add_argument('--out', required=True)
ap.add_argument('--cache', default=None, help='trident residues per prime (default <tridents>/residues)')
ap.add_argument('-c', '--corners', default='1,2,3')
ap.add_argument('--order', type=int, default=None, help='final order N (default: the largest the tridents allow)')
ap.add_argument('-n', '--jobs', type=int, default=os.cpu_count())
ap.add_argument('--nprimes', type=int, default=16, help='initial number of primes (more are added if needed)')
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
CACHE = a.cache or os.path.join(a.tridents, 'residues')
os.makedirs(CACHE, exist_ok=True)
with open(os.path.join(a.tridents, '.order')) as f:
    ORDt, ORDy = (int(x) for x in f.read().split())
NtT, NyT = ORDt + 2, ORDy + 2                      # trident grids: n = -1..ORD
NTRI = 1318
PI_KEY = (('Pi', 2),)

t, u, T, B, S, G = sp.symbols('t u T B S G')
# (t, u) = (y, z) in the blow-up chart of each point ...
CHART = {1: {t: B, u: 1 - B - T * B**2},          # x = T B^2, y = B
         2: {t: B, u: T * B**2},                  # z = T B^2, y = B
         3: {t: 1 - B - T * B**2, u: B}}          # x = T B^2, z = B
# ... and in its physical variables (S = SMALL, G = BIG), to tell physical from spurious poles
PHYS = {1: {t: G, u: 1 - S - G}, 2: {t: G, u: S}, 3: {t: 1 - S - G, u: G}}
VARS = {1: ('x', 'y'), 2: ('z', 'y'), 3: ('x', 'z')}             # (SMALL, BIG), as phys.MAP


# ================================================================ coefficients
TERM = re.compile(r'^\(\((.*)\)/\((.*)\)\)\*Trident\[(\d+), *(\d+)\](?: \+)?$')


def parse(path):
    """-> [(block, i, w, N, D)], block 're' or 'ipi' (inside I*Pi*( ... )), N, D sympy Polys in t, u."""
    out, block = [], 're'
    for k, line in enumerate(open(path)):
        line = line.strip()
        if not line:
            continue
        if line == 'I*Pi*(':
            block = 'ipi'; continue
        if line == ')':
            block = 're'; continue
        m = TERM.match(line)
        if not m:
            sys.exit("%s:%d: unexpected line %r" % (path, k + 1, line[:120]))
        N, D = (sp.Poly(sp.sympify(s.replace('^', '**')), t, u, domain='QQ') for s in m.group(1, 2))
        out.append((block, int(m.group(3)), int(m.group(4)), N, D))
    return out


def to_dict(expr, x1, x2):
    P = sp.Poly(sp.expand(expr), x1, x2)
    return {mono: mpq(int(c.p), int(c.q)) for mono, c in zip(P.monoms(), P.coeffs())}


def split(d):
    """{(i,j): q} -> (a, b, {(i-a, j-b): q}) with the smallest powers pulled out."""
    a_ = min(i for i, _ in d); b_ = min(j for _, j in d)
    return a_, b_, {(i - a_, j - b_): q for (i, j), q in d.items()}


def chart_terms(terms, c):
    """[(block, i, w, shift_e, shift_m, Ntil, Dtil_key)] and the distinct Dtil."""
    dcache, dtil, out = {}, {}, []
    for block, i, w, N, D in terms:
        key = D.as_expr()
        if key not in dcache:
            ad, bd, Dt = split(to_dict(D.as_expr().subs(CHART[c], simultaneous=True), T, B))
            if Dt.get((0, 0), 0) == 0:
                sys.exit("point %d: the chart does not resolve the denominator %s" % (c, D.as_expr()))
            dk = tuple(sorted(Dt.items()))
            dtil[dk] = Dt
            dcache[key] = (ad, bd, dk)
        ad, bd, dk = dcache[key]
        an, bn, Nt = split(to_dict(N.as_expr().subs(CHART[c], simultaneous=True), T, B))
        out.append((block, i, w, an - ad, bn - bd, Nt, dk))
    return out, dtil


def physical_bound(terms, c):
    """lowest powers of (SMALL, BIG) the PHYSICAL denominators allow (spurious factors excluded)."""
    lo = [0, 0]
    seen = {}
    for _, _, _, _, D in terms:
        key = D.as_expr()
        if key in seen:
            continue
        ps = [0, 0]
        for f, mult in sp.factor_list(D.as_expr())[1]:
            d = to_dict(f.subs(PHYS[c], simultaneous=True), S, G)
            if d.get((0, 0), 0) != 0:
                continue                                  # a unit at the point
            if len(d) == 1:                               # a monomial: physical pole
                (i, j), = d
                ps[0] -= i * mult; ps[1] -= j * mult
            # otherwise it vanishes at the point without being a monomial: spurious
        seen[key] = ps
        lo = [min(lo[0], ps[0]), min(lo[1], ps[1])]
    return lo


# ================================================================ residues
def inv_series(Dt, p, ne, nm):
    """1/Dtil mod p on the grid [0, ne) x [0, nm)."""
    d0 = pow(int(ba.rmod(Dt[(0, 0)], p)), -1, p)
    terms = [(i, j, int(ba.rmod(q, p))) for (i, j), q in Dt.items() if (i, j) != (0, 0) and i < ne and j < nm]
    C = np.zeros((ne, nm), np.int64)
    for i in range(ne):
        for j in range(nm):
            s = 1 if (i, j) == (0, 0) else 0
            for k, l, q in terms:
                if k <= i and l <= j:
                    s -= q * int(C[i - k, j - l])
            C[i, j] = s * d0 % p
    return C


def cache_path(c, w, k):
    return os.path.join(CACHE, 'T_c%d_w%d_k%03d.pkl' % (c, w, k))


def fill_cache(job):
    """load the trident series of (c, w) ONCE and write the residues for every prime in ks."""
    c, w, ks = job
    slots = pc.load_pickle(os.path.join(a.tridents, 'trident_c%d_w%d.pkl' % (c, w)))
    bytag = {}
    for (tag, i, it, iy), q in slots.items():
        bytag.setdefault(tag, []).append((i, it, iy, q))
    del slots
    for k in ks:                                  # one prime at a time: bounded memory
        p = pc.PRIME_TABLE[k]
        out = {}
        for tag, ents in bytag.items():
            arr = out[tag] = np.zeros((NTRI, NtT, NyT), np.int64)
            for i, it, iy, q in ents:
                arr[i, it, iy] = ba.rmod(q, p)
        pc.save_residues(out, cache_path(c, w, k))
    return c, w, len(ks)


def prepare_cache(c, weights, ks):
    jobs = [(c, w, [k for k in ks if not os.path.isfile(cache_path(c, w, k))]) for w in weights]
    jobs = [j for j in jobs if j[2]]
    if not jobs:
        return
    t0 = time.time()
    with mproc.get_context('fork').Pool(min(a.jobs, len(jobs))) as pool:
        for c_, w, n in pool.imap_unordered(fill_cache, jobs):
            pass
    print("    trident residues cached for weights %s, %d primes  (%.1fs)"
          % (sorted(j[1] for j in jobs), len(ks), time.time() - t0), flush=True)


def trident_residues(c, w, k):
    """{tag: int64 (NTRI, NtT, NyT)} of the tridents at point c, weight w, prime index k (cached)."""
    return pc.load_residues(cache_path(c, w, k))


_G = {}


def contract(k):
    """A mod prime k: {tagA: int64 (1, NA, MA)}, or None if the prime is unusable."""
    p = pc.PRIME_TABLE[k]
    q_, p_, NA, MA = _G['q'], _G['p'], _G['NA'], _G['MA']
    ne, nm = NA - 1, MA - 1                          # M grid: (e, m) = (ke - q, km - p)
    try:
        inv = {dk: inv_series(Dt, p, ne, nm) for dk, Dt in _G['dtil'].items()}
    except (ZeroDivisionError, ValueError):
        return k, None
    A = {}
    for (block, w), items in _G['groups'].items():
        TR = trident_residues(_G['c'], w, k)
        MS = np.zeros((len(items), ne, nm), np.int64)
        for r, (i, se, sm, Nt, dk) in enumerate(items):
            e0, m0 = se + q_, sm + p_                # grid position of T^se B^sm
            Dinv = inv[dk]
            for (ni, nj), qv in Nt.items():
                oe, om = e0 + ni, m0 + nj
                if oe >= ne or om >= nm:
                    continue
                MS[r, oe:, om:] = (MS[r, oe:, om:] + int(ba.rmod(qv, p)) * Dinv[:ne - oe, :nm - om]) % p
        MS = MS.reshape(len(items), ne * nm)
        live = np.flatnonzero(MS.any(axis=0))       # chart shifts that occur
        MSb = np.where(MS[:, live] > p // 2, MS[:, live] - p, MS[:, live]).astype(np.float64).T
        idx = np.array([it[0] - 1 for it in items])
        for tag, arr in TR.items():
            sub = arr[idx].reshape(len(items), NtT * NyT)
            if not sub.any():
                continue
            subb = np.where(sub > p // 2, sub - p, sub).astype(np.float64)
            R = np.remainder(MSb @ subb, p).astype(np.int64)          # exact: < 1318 * 2^42
            pt, py, la, lb, ck = tag
            sign = 1
            if block == 'ipi':                                       # times I * Pi
                ck, s1 = ba.key_mul(ck, ba.I_KEY)
                ck, s2 = ba.key_mul(ck, PI_KEY)
                sign = s1 * s2
            tA = (pt, py, la, lb, ck)
            acc = A.get(tA)
            if acc is None:
                acc = A[tA] = np.zeros((NA, MA), np.int64)
            for r, sh in enumerate(live):
                ke, km = divmod(int(sh), nm)                         # A index = T index + (ke, km)
                he, hm = min(NtT, NA - ke), min(NyT, MA - km)
                if he <= 0 or hm <= 0:
                    continue
                blk = R[r].reshape(NtT, NyT)[:he, :hm]
                acc[ke:ke + he, km:km + hm] += blk if sign == 1 else p - blk
                acc[ke:ke + he, km:km + hm] %= p
    return k, {tg: v[None] for tg, v in A.items() if v.any()}


# ================================================================ output
def gpow_series(alpha, n):
    """(1 - Exp[-s])/s to the power alpha (mpq), coefficients of s^0..s^n."""
    g = [mpq((-1) ** j, math.factorial(j + 1)) for j in range(n + 1)]
    f = [mpq(1)]
    for m in range(1, n + 1):
        s = mpq(0)
        for j in range(1, m + 1):
            s += ((alpha + 1) * j - m) * g[j] * f[m - j]
        f.append(s / m)
    return f


def bernoulli(physd, sm, bg, N):
    """phys terms {(ex2, ey2, lx, ly, ck): q} -> Bernoulli poly {key: q}, order N in each s."""
    cache, out = {}, {}
    lg = {'x': 'LogX', 'y': 'LogY', 'z': 'LogZ'}
    for (ex2, ey2, lx, ly, ck), q in physd.items():
        ser = []
        for e2 in (ex2, ey2):
            if e2 not in cache:
                cache[e2] = gpow_series(mpq(e2, 2), max(0, 2 * N - e2) // 2 + 1)
            ser.append(cache[e2])
        for kx, cx in enumerate(ser[0]):
            if ex2 + 2 * kx > 2 * N:
                break
            if not cx:
                continue
            for ky, cy in enumerate(ser[1]):
                if ey2 + 2 * ky > 2 * N:
                    break
                if not cy:
                    continue
                key = list(ck)
                if ex2 + 2 * kx: key.append(('s' + sm, ex2 + 2 * kx))
                if ey2 + 2 * ky: key.append(('s' + bg, ey2 + 2 * ky))
                if lx: key.append((lg[sm], 2 * lx))
                if ly: key.append((lg[bg], 2 * ly))
                key = tuple(sorted(key))
                v = out.get(key, 0) + q * cx * cy
                if v:
                    out[key] = v
                else:
                    out.pop(key, None)
    return out


def write_expr(poly, path, header):
    def order(item):
        d = dict(item[0])
        return tuple(d.get(v, 0) for v in ('x', 'y', 'z', 'sx', 'sy', 'sz')) + (item[0],)
    with open(path, 'w') as f:
        f.write("(* " + header.replace("*)", "* )") + " *)\n")
        f.write(mp.fmt_poly(poly, order=order) + "\n")


# ================================================================ main
for path in a.files:
    name = os.path.splitext(os.path.basename(path))[0]
    t0 = time.time()
    terms = parse(path)
    print("=== %s: %d terms (%d in I*Pi), %d distinct denominators  (%.1fs)"
          % (name, len(terms), sum(1 for x in terms if x[0] == 'ipi'),
             len({x[4].as_expr() for x in terms}), time.time() - t0), flush=True)
    for c in [int(s) for s in a.corners.split(',')]:
        t1 = time.time()
        ct, dtil = chart_terms(terms, c)
        q_ = max(0, -min(x[3] for x in ct)); p_ = max(0, -min(x[4] for x in ct))
        At, Am = ORDt - q_, ORDy - p_
        N = min(At, Am // 3)
        if a.order is not None:
            if a.order > N:
                sys.exit("order %d needs tridents to %d/%d (have %d/%d)" % (a.order, a.order + q_, 3 * a.order + p_, ORDt, ORDy))
            N = a.order
        lo = physical_bound(terms, c)
        groups = {}
        for block, i, w, se, sm_, Nt, dk in ct:
            groups.setdefault((block, w), []).append((i, se, sm_, Nt, dk))
        _G.update(c=c, q=q_, p=p_, NA=At + q_ + 2, MA=Am + p_ + 2, dtil=dtil, groups=groups)
        print("  point %d %s: poles T^-%d B^-%d -> final order N = %d; physical poles allow powers >= %s"
              % (c, VARS[c], q_, p_, N, lo), flush=True)
        res, ks, kn = {}, [], 0
        target = a.nprimes + ba.NV
        while True:
            todo = list(range(kn, kn + target - len(ks))); kn += len(todo)
            prepare_cache(c, sorted({w for _, w in groups}), todo)
            with mproc.get_context('fork').Pool(min(a.jobs, len(todo))) as pool:
                for k, r in pool.imap_unordered(contract, todo):
                    if r is None:
                        print("    prime k=%d skipped" % k)
                    else:
                        res[k] = r; ks.append(k)
            if len(ks) < target:
                continue
            ks.sort()
            out, bad, _ = ba.reconstruct([res[k] for k in ks], [pc.PRIME_TABLE[k] for k in ks], jobs=a.jobs)
            if not bad:
                break
            target += max(4, a.jobs // 2)
            print("    %d coefficients unverified with %d primes -> %d primes" % (bad, len(ks) - ba.NV, target - ba.NV), flush=True)
        # chart slots -> physical variables (rectangle of order N in both)
        slots = {(tag, 0, ia - q_, ja - p_): v for (tag, _, ia, ja), v in out.items()}
        pd = phys.slots_to_phys(SimpleNamespace(c=c), slots, N, N)
        pd = {k[1:]: v for k, v in pd.items()}
        sm, bg = VARS[c]
        amin = mpq(min((k[0] for k in pd), default=0), 2); bmin = mpq(min((k[1] for k in pd), default=0), 2)
        below = sum(1 for k in pd if k[0] < 2 * lo[0] or k[1] < 2 * lo[1])
        status = "spurious poles cancelled" if not below else "!! %d TERMS BELOW THE PHYSICAL POLES" % below
        hdr = ("%s at %s,%s -> 0, physical variables, rectangle of order %d in %s and %s (lowest powers %s^%s, %s^%s); "
               "constants Pi, Zeta[3], Zeta[5], I; from tridents at blow-up order %d/%d (assemble_amplitude.py)"
               % (name, sm, bg, N, sm, bg, sm, amin, bg, bmin, ORDt, ORDy))
        poly = {}
        for (ex2, ey2, lx, ly, ck), v in pd.items():
            key = list(ck)
            if ex2: key.append((sm, ex2))
            if ey2: key.append((bg, ey2))
            if lx: key.append(('Log[%s]' % sm, 2 * lx))
            if ly: key.append(('Log[%s]' % bg, 2 * ly))
            poly[tuple(sorted(key))] = v
        write_expr(poly, os.path.join(a.out, '%s_c%d_xyz_ord%d.m' % (name, c, N)), hdr)
        bp = bernoulli(pd, sm, bg, N)
        write_expr(bp, os.path.join(a.out, '%s_c%d_ber_ord%d.m' % (name, c, N)),
                   "%s at %s,%s -> 0, Bernoulli variables %s -> 1 - Exp[-s%s], %s -> 1 - Exp[-s%s], order %d in s%s "
                   "and s%s; Log[%s] -> Log%s, Log[%s] -> Log%s are left as placeholders (substitute at the numeric point)"
                   % (name, sm, bg, sm, sm, bg, bg, N, sm, bg, sm, sm.upper(), bg, bg.upper()))
        print("  point %d: %d chart coefficients, %d primes (+%d) | phys: %d terms, lowest powers %s^%s %s^%s, %s | "
              "Bernoulli: %d terms  (%.1fs)"
              % (c, len(out), len(ks) - ba.NV, ba.NV, len(pd), sm, amin, bg, bmin, status, len(bp), time.time() - t1),
              flush=True)
