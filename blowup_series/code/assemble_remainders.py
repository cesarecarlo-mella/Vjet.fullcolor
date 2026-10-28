#!/usr/bin/env python3
"""
assemble_remainders.py -- series of the 84 finite remainders of Zgqq_finite_remainders (alpha1-3, gamma1-3, 14 colour
layers each) at the three decay corners, from the series of the iterated-integral system (family ITER).

Every remainder is  sum_rows q * R_r(t, u) * T_g * c_m  (t = y, u = z; generators.wl: T_g = sum_i c_gi F_i, the functions
of the ITER family; c_m a constant, or I/r with r = +Sqrt[x y z]).  The rational functions have spurious poles at the
corners, so everything is multiplied in the blow-up chart of the corner and only the result goes back:

    corner c1 (x, y -> 0): chart 1, x = T y^2, BIG = y        variables (x, y)
    corner c2 (y, z -> 0): chart 5, y = T z^2, BIG = z        variables (y, z)   (the chart of the ITER constants at c2)
    corner c3 (x, z -> 0): chart 3, x = T z^2, BIG = z        variables (x, z)

In every chart the third variable is U = 1 - B - T B^2, a unit, and r = T^(1/2) B^(3/2) U^(1/2).  R_r = T^a B^b N/D with
D(0,0) != 0: the quadratic letters that the charts 1 and 3 do not resolve, (1 - T + ...), are units at the corner too
(their pole at T = 1 must cancel, and does: no term below the physical poles survives).

Per prime: the function series (cached per weight), for each layer and constant one BLAS product V . F (V: rational
function x function coefficients), the shifts of the numerators and monomials, and one convolution with 1/D (times
U^(-1/2) for the I/r rows) per distinct denominator.  CRT + rational reconstruction verified on 2 extra primes (more
are added per layer if needed).  Output: <out>/<coefficient>/<layer>/series<c>.m in the variables of the corner (order N
in both) and seriesBer<c>.m in the Bernoulli variables v -> 1 - Exp[-s_v] (LogX, LogY, LogZ placeholders).

A final order N needs the ITER series to (N + q, 3N + p), q and p the largest poles (T^-3.5, B^-10.5 here): 24/72 for N = 20.

    python3 assemble_remainders.py --iter ../families/ITER/dist_24_72 --inputs ../families/ITER/inputs \
        --generators "../../iter int /generators_deq/out/generators.wl" --zgqq ../../Zgqq_finite_remainders \
        --out ../../remainder_series -n 32
"""
import argparse
import math
import multiprocessing as mproc
import os
import re
import sys
import time
from fractions import Fraction
from types import SimpleNamespace

import numpy as np
import sympy as sp
from gmpy2 import mpq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import blowup_array as ba      # noqa: E402
import mparse as mp            # noqa: E402
import phys                    # noqa: E402
import pipe_common as pc       # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('--iter', required=True, help='run directory of the ITER family (rec_c<c>_w<w>.pkl, .order)')
ap.add_argument('--inputs', required=True, help='ITER inputs (labels.m, sol0.m)')
ap.add_argument('--generators', required=True, help='generators.wl: the generators in terms of the functions')
ap.add_argument('--zgqq', required=True, help='Zgqq_finite_remainders directory')
ap.add_argument('--out', required=True)
ap.add_argument('-c', '--charts', default='1,3,5')
ap.add_argument('--order', type=int, default=None, help='final order N (default: the largest the series allow)')
ap.add_argument('--only', default=None, help='comma list of coefficient/layer, e.g. alpha1/C3F0Z0 (default: all 84)')
ap.add_argument('-n', '--jobs', type=int, default=2)
ap.add_argument('--nprimes', type=int, default=12, help='initial number of primes (more are added per layer if needed)')
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
WORK = os.path.join(a.out, '.work')
with open(os.path.join(a.iter, '.order')) as f:
    ORDt, ORDy = (int(x) for x in f.read().split())
NtF, NyF = ORDt + 1, ORDy + 1                      # function grids: n = 0..ORD (rec slots it, iy = n + 1)
POINT = {1: 'c1', 5: 'c2', 3: 'c3'}
VARS = {1: ('x', 'y'), 5: ('y', 'z'), 3: ('x', 'z')}     # (SMALL, BIG) = phys.MAP
t, u, T, B, S, G = sp.symbols('t u T B S G')
U = 1 - B - T * B**2
CHART = {1: {t: B, u: U}, 3: {t: U, u: B}, 5: {t: T * B**2, u: B}}
PHYS = {1: {t: G, u: 1 - S - G}, 3: {t: 1 - S - G, u: G}, 5: {t: S, u: G}}
LAYERNAME = lambda key: '_'.join('%s%s' % (n, v) for n, v in zip(('Nc', 'Nf', 'NfZ'),
                                  re.fullmatch(r'Nc\^(-?\d+)\*Nf\^(\d+)\*NfZ\^(\d+)', key).groups()))

# ================================================================ inputs
labels = re.findall(r'"([^"]+)"', open(os.path.join(a.inputs, 'labels.m')).read())
LIDX = {l: i for i, l in enumerate(labels)}
sol0 = [Fraction(s.strip()) for s in mp.split_top_list(open(os.path.join(a.inputs, 'sol0.m')).read())]
GEN = {}
for line in open(a.generators):
    m = re.match(r'(T\[\d+\]) -> (.*?),?\s*$', line.strip())
    if m:
        GEN[m.group(1)] = [(Fraction(c), LIDX[l]) for c, l in re.findall(r'\((-?[\d/]+)\)\*([A-Z]\[[\d,]+\])', m.group(2))]


def _items(fn):
    src = re.sub(r'\(\*.*?\*\)', '', open(os.path.join(a.zgqq, fn)).read(), flags=re.S).strip()
    return [s.strip() for s in src[1:-1].split(',')]


RFB = [sp.sympify(s.replace('^', '**'), locals={'t': t, 'u': u}) for s in _items('RationalFunctionBasis.wl')]
TM = []                                            # (constant key, power of r in exp2 units, rational factor)
for s in _items('TranscendentalMultipliers.wl'):
    (key, c), = mp.parse_string(s).items()
    d = dict(key)
    kr = d.pop('r', 0)
    TM.append((tuple(sorted(d.items())), kr, c))
LAYERS = []                                         # (coefficient, layer key, rows)
for name in ('alpha1', 'alpha2', 'alpha3', 'gamma1', 'gamma2', 'gamma3'):
    src = open(os.path.join(a.zgqq, 'FiniteRemainders', name + '.wl')).read()
    for key, rows in re.findall(r'"([^"]+)"\s*->\s*\{(.*?)\}\s*(?=,\s*"|\s*\}\s*$)', src[src.index('*)') + 2:], re.S):
        LAYERS.append((name, key, [(int(r) - 1, 'T[%s]' % g, int(m_) - 1, Fraction(q))
                                   for r, g, m_, q in re.findall(r'\{(\d+),(\d+),(\d+),(-?\d+(?:/\d+)?)\}', rows)]))
if a.only:
    want = {tuple(s.split('/')) for s in a.only.split(',')}
    LAYERS = [L for L in LAYERS if (L[0], L[1]) in want or (L[0], LAYERNAME(L[1])) in want or
              (L[0], 'C%sF%sZ%s' % tuple(str(int(v) % 10) for v in re.findall(r'-?\d+', L[1]))) in want]
# per layer: {(constant key, r power): {rational index: {function index: coefficient}}}
COMB = []
for name, key, rows in LAYERS:
    comb = {}
    for r, g, m_, q in rows:
        ck, kr, cm = TM[m_]
        vec = comb.setdefault((ck, kr), {}).setdefault(r, {})
        for c, i in GEN[g]:
            vec[i] = vec.get(i, 0) + q * Fraction(int(cm.numerator), int(cm.denominator)) * c
    COMB.append({g_: {r: {i: v for i, v in vec.items() if v} for r, vec in d.items()} for g_, d in comb.items()})


def weight_of(lab):
    m = re.match(r'[BS]\[(\d+),', lab)
    if m:
        return int(m.group(1))
    return None


# the weight of every function: from the rec files (a function lives in the series of its weight) and B[0,1]
WEIGHT = {}


# ================================================================ rational functions in a chart
def to_dict(expr, x1, x2):
    P = sp.Poly(sp.expand(expr), x1, x2)
    return {mono: mpq(int(c.p), int(c.q)) for mono, c in zip(P.monoms(), P.coeffs())}


def split(d):
    a_ = min(i for i, _ in d); b_ = min(j for _, j in d)
    return a_, b_, {(i - a_, j - b_): q for (i, j), q in d.items()}


def chart_rationals(c):
    """for every rational function: (shifts [(s2T, s2B, coef)] with the monomial, denominator key, Dtil)"""
    out, dt = [], {}
    for e in RFB:
        num, den = sp.fraction(sp.together(e))
        an, bn, Nt = split(to_dict(num.subs(CHART[c], simultaneous=True), T, B))
        ad, bd, Dt = split(to_dict(den.subs(CHART[c], simultaneous=True), T, B))
        if Dt.get((0, 0), 0) == 0:
            sys.exit('chart %d does not resolve the denominator %s' % (c, den))
        dk = tuple(sorted(Dt.items()))
        dt[dk] = Dt
        out.append(([(2 * (an - ad + i), 2 * (bn - bd + j), q) for (i, j), q in Nt.items()], dk))
    return out, dt


def physical_bound(c):
    """lowest powers (exp2 units) of (SMALL, BIG) the physical denominators allow, and the extra -1/2, -1/2 of 1/r"""
    lo = [0, 0]
    for e in RFB:
        den = sp.fraction(sp.together(e))[1]
        ps = [0, 0]
        for f, mult in sp.factor_list(den)[1]:
            d = to_dict(f.subs(PHYS[c], simultaneous=True), S, G)
            if d.get((0, 0), 0) != 0 or len(d) != 1:
                continue                               # a unit, or spurious (vanishes but not a monomial)
            (i, j), = d
            ps[0] -= 2 * i * mult; ps[1] -= 2 * j * mult
        lo = [min(lo[0], ps[0]), min(lo[1], ps[1])]
    return lo


def sqrt_inv_U(ne, nm):
    """U^(-1/2) = (1 - w)^(-1/2), w = B + T B^2, exact on the grid [0, ne) x [0, nm)"""
    out = {(0, 0): mpq(1)}
    wk = {(0, 0): mpq(1)}
    binom = mpq(1)
    for k in range(1, nm + 1):
        nw = {}
        for (i, j), q in wk.items():
            for (di, dj) in ((0, 1), (1, 2)):
                if i + di < ne and j + dj < nm:
                    nw[(i + di, j + dj)] = nw.get((i + di, j + dj), 0) + q
        wk = nw
        if not wk:
            break
        binom = binom * mpq(2 * k - 1, 2 * k)           # (-1)^k binom(-1/2, k)
        for key, q in wk.items():
            out[key] = out.get(key, 0) + binom * q
    return out


# ================================================================ per prime
def inv_series(Dt, p, ne, nm):
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


def toeplitz_stack(S2, ne, nm, p):
    """X -> X * S2 (2-d series product, truncated) row by row: T[di][j, J] = S2[di, J - j] (J >= j), balanced float64"""
    Sb = np.where(S2 > p // 2, S2 - p, S2).astype(np.float64)
    T = np.zeros((ne, nm, nm))
    for dj in range(nm):
        r = np.arange(nm - dj)
        T[:, r, r + dj] = Sb[:, dj][:, None]
    return T


def series_product(Xb, T, p):
    """Xb (ntags, ne, nm) balanced float64 times the series of T -> int64 mod p (each product sums nm terms: exact)"""
    ntg, ne, nm = Xb.shape
    out = np.zeros((ntg, ne, nm), np.int64)
    for di in range(ne):
        if T[di].any():
            part = np.remainder(Xb[:, :ne - di, :] @ T[di], p).astype(np.int64)
            out[:, di:, :] = (out[:, di:, :] + part) % p
    return out


def bal(x, p):
    return np.where(x > p // 2, x - p, x).astype(np.float64)


_G = {}


def fcache(c, w, k):
    return os.path.join(WORK, 'F_c%d_w%d_k%03d.pkl' % (c, w, k))


def fill_cache(job):
    """the function series of (c, w), mod every prime in ks: {tag: int64 (nfun_w, NtF, NyF)}"""
    c, w, ks = job
    idx = _G['idx'][w]
    pos = {i: n for n, i in enumerate(idx)}
    if w == 0:
        bytag = {(0, 0, 0, 0, ()): [(pos[i], 1, 1, mpq(int(v.numerator), int(v.denominator))) for i, v in enumerate(sol0) if v]}
    else:
        bytag = {}
        for (tag, i, it, iy), q in pc.load_pickle(os.path.join(a.iter, 'rec_c%d_w%d.pkl' % (c, w))).items():
            bytag.setdefault(tag, []).append((pos[i], it, iy, q))
    for k in ks:
        p = pc.PRIME_TABLE[k]
        out = {}
        for tag, ents in bytag.items():
            arr = out[tag] = np.zeros((len(idx), NtF, NyF), np.int64)
            for n, it, iy, q in ents:
                arr[n, it - 1, iy - 1] = ba.rmod(q, p)
        pc.save_residues(out, fcache(c, w, k))
    return w


def acache(c, li, k):
    return os.path.join(WORK, 'A_c%d_L%02d_k%03d.pkl' % (c, li, k))


def contract(job):
    """prime k: the residues of the layers li in todo -> cache files"""
    k, todo = job
    p = pc.PRIME_TABLE[k]
    c, NA, MA, QT, QB = _G['c'], _G['NA'], _G['MA'], _G['QT'], _G['QB']
    rat, dtil = _G['rat'], _G['dtil']
    F = {w: pc.load_residues(fcache(c, w, k)) for w in _G['idx']}
    Fb = {w: {tag: bal(v, p) for tag, v in F[w].items()} for w in F}
    try:
        inv = {dk: inv_series(Dt, p, NA, MA) for dk, Dt in dtil.items()}
    except (ZeroDivisionError, ValueError):
        return k, False
    U12 = {kk: int(ba.rmod(q, p)) for kk, q in _G['u12'].items()}
    Cm = {}

    def cmat(dk, half):
        if (dk, half) not in Cm:
            S2 = inv[dk]
            if half:                                   # times U^(-1/2)
                Uh = np.zeros((NA, MA), np.int64)
                for (i, j), v in U12.items():
                    Uh[i, j] = v
                S3 = np.zeros((NA, MA), np.int64)
                for (i, j), v in U12.items():
                    if v:
                        S3[i:, j:] = (S3[i:, j:] + v * S2[:NA - i, :MA - j]) % p
                S2 = S3
            Cm[(dk, half)] = toeplitz_stack(S2, NA, MA, p)
        return Cm[(dk, half)]

    for li in todo:
        X = {}                                         # (dk, half) -> {tag': (NA, MA) int64}
        for (ckm, kr), rvecs in COMB[li].items():
            half = kr != 0
            rows = sorted(rvecs)
            for w, idx in _G['idx'].items():
                if not F[w]:
                    continue
                pos = {i: n for n, i in enumerate(idx)}
                V = np.zeros((len(rows), len(idx)))
                for rr, r in enumerate(rows):
                    for i, v in rvecs[r].items():
                        if i in pos:
                            V[rr, pos[i]] = int(ba.rmod(mpq(v.numerator, v.denominator), p))
                if not V.any():
                    continue
                V = np.where(V > p // 2, V - p, V)
                tags = list(Fb[w])
                Fw = np.stack([Fb[w][tg] for tg in tags], axis=1).reshape(len(idx), -1)
                Gw = np.remainder(V @ Fw, p).astype(np.int64).reshape(len(rows), len(tags), NtF, NyF)
                for rr, r in enumerate(rows):
                    shifts, dk = rat[r]
                    key = (dk, half)
                    Xd = X.setdefault(key, {})
                    for jt, tg in enumerate(tags):
                        sub = Gw[rr, jt]
                        if not sub.any():
                            continue
                        pt, py, la, lb, ck = tg
                        ck2, sg = ba.key_mul(ck, ckm)
                        for s2T, s2B, q in shifts:
                            s2T -= 1 if half else 0
                            s2B -= 3 if half else 0
                            eT, eB = pt + s2T, py + s2B
                            tg2 = (eT % 2, eB % 2, la, lb, ck2)
                            dn, dm = eT // 2 + QT, eB // 2 + QB
                            i0, j0 = max(0, dn), max(0, dm)
                            i1, j1 = min(NA, dn + NtF), min(MA, dm + NyF)
                            if i1 <= i0 or j1 <= j0:
                                continue
                            qv = int(ba.rmod(q, p)) * sg % p
                            acc = Xd.get(tg2)
                            if acc is None:
                                acc = Xd[tg2] = np.zeros((NA, MA), np.int64)
                            acc[i0:i1, j0:j1] = (acc[i0:i1, j0:j1] + qv * sub[i0 - dn:i1 - dn, j0 - dm:j1 - dm]) % p
        A = {}
        for (dk, half), Xd in X.items():
            tg2s = list(Xd)
            R = series_product(np.stack([bal(Xd[tg], p) for tg in tg2s]), cmat(dk, half), p)
            for tg, blk in zip(tg2s, R):
                acc = A.get(tg)
                A[tg] = blk if acc is None else (acc + blk) % p
        pc.save_residues({tg: v[None] for tg, v in A.items() if v.any()}, acache(c, li, k))
    return k, True


# ================================================================ output
def gpow_series(alpha, n):
    g = [mpq((-1) ** j, math.factorial(j + 1)) for j in range(n + 1)]
    f = [mpq(1)]
    for m in range(1, n + 1):
        s = mpq(0)
        for j in range(1, m + 1):
            s += ((alpha + 1) * j - m) * g[j] * f[m - j]
        f.append(s / m)
    return f


def bernoulli(physd, sm, bg, N):
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
        f.write((mp.fmt_poly(poly, order=order) if poly else '0') + "\n")


# ================================================================ main
for c in [int(s) for s in a.charts.split(',')]:
    t0 = time.time()
    pt_, (sm, bg) = POINT[c], VARS[c]
    # functions by weight (from the rec files)
    idx = {0: [LIDX['B[0,1]']]}
    for w in range(1, 7):
        S_ = pc.load_pickle(os.path.join(a.iter, 'rec_c%d_w%d.pkl' % (c, w)))
        idx[w] = sorted({i for (_, i, _, _) in S_} | {LIDX[l] for l in labels if weight_of(l) == w})
    rat, dtil = chart_rationals(c)
    s2T = min(s for sh, _ in rat for s, _, _ in sh) - 1          # -1: T^(-1/2) of 1/r
    s2B = min(s for sh, _ in rat for _, s, _ in sh) - 3          # -3: B^(-3/2) of 1/r
    QT, QB = -(s2T // 2), -(s2B // 2)
    NA, MA = ORDt + 1 + QT, ORDy + 1 + QB                        # n' in [-QT, ORDt], [-QB, ORDy]
    At, Am = ORDt + s2T // 2, ORDy + s2B // 2                   # complete up to these powers
    N = min(At, Am // 3)
    if a.order is not None:
        if a.order > N:
            sys.exit('order %d needs the ITER series to %d/%d (have %d/%d)' % (a.order, a.order - s2T // 2, 3 * a.order - s2B // 2, ORDt, ORDy))
        N = a.order
    assert (max(len(v) for v in idx.values()) + MA) * (pc.PRIME_TABLE[0] // 2) ** 2 < 2 ** 53   # exact float64 products
    lo = physical_bound(c)
    lo = [lo[0] - 1, lo[1] - 1]                                # 1/r: -1/2 in SMALL and in BIG (exp2 units)
    _G.update(c=c, NA=NA, MA=MA, QT=QT, QB=QB, rat=rat, dtil=dtil, idx=idx, u12=sqrt_inv_U(NA, MA))
    os.makedirs(WORK, exist_ok=True)
    print('=== corner %s (chart %d, variables %s, %s): %d layers; poles T^%s B^%s -> final order N = %d; ITER series %d/%d; '
          '%d denominators  (%.1fs)' % (pt_, c, sm, bg, len(LAYERS), mpq(s2T, 2), mpq(s2B, 2), N, ORDt, ORDy, len(dtil), time.time() - t0), flush=True)
    pending = list(range(len(LAYERS)))
    ks = []
    kn = 0
    target = a.nprimes + ba.NV
    while pending:
        new = list(range(kn, kn + target - len(ks)))
        kn += len(new)
        t1 = time.time()
        with mproc.get_context('fork').Pool(min(a.jobs, 7)) as pool:
            list(pool.imap_unordered(fill_cache, [(c, w, new) for w in idx]))
        with mproc.get_context('fork').Pool(min(a.jobs, len(new))) as pool:
            for k, ok in pool.imap_unordered(contract, [(k, pending) for k in new]):
                if ok:
                    ks.append(k)
                else:
                    print('    prime k=%d skipped' % k)
        ks.sort()
        print('  %d primes for %d layers  (%.1fs)' % (len(ks), len(pending), time.time() - t1), flush=True)
        still = []
        for li in pending:
            name, key, _ = LAYERS[li]
            res = [pc.load_residues(acache(c, li, k)) for k in ks]
            out, bad, _ = ba.reconstruct(res, [pc.PRIME_TABLE[k] for k in ks], jobs=1)
            if bad:
                still.append(li)
                continue
            slots = {(tag, 0, ia - QT + 1, ja - QB + 1): v for (tag, _, ia, ja), v in out.items()}
            pd = phys.slots_to_phys(SimpleNamespace(c=c), slots, N, N)
            pd = {kk[1:]: v for kk, v in pd.items()}
            below = sum(1 for kk in pd if kk[0] < lo[0] or kk[1] < lo[1])
            amin = mpq(min((kk[0] for kk in pd), default=0), 2); bmin = mpq(min((kk[1] for kk in pd), default=0), 2)
            poly = {}
            for (ex2, ey2, lx, ly, ck), v in pd.items():
                kk = list(ck)
                if ex2: kk.append((sm, ex2))
                if ey2: kk.append((bg, ey2))
                if lx: kk.append(('Log[%s]' % sm, 2 * lx))
                if ly: kk.append(('Log[%s]' % bg, 2 * ly))
                poly[tuple(sorted(kk))] = v
            d = os.path.join(a.out, name, LAYERNAME(key))
            os.makedirs(d, exist_ok=True)
            write_expr(poly, os.path.join(d, 'series%s.m' % pt_),
                       '%s, colour layer %s, at %s,%s -> 0 (corner %s): rectangle of order %d in %s and %s (lowest powers %s^%s, '
                       '%s^%s); constants Pi, Zeta[3], Zeta[5], I; from the ITER series at blow-up order %d/%d, chart %d '
                       '(assemble_remainders.py)%s'
                       % (name, key, sm, bg, pt_, N, sm, bg, sm, amin, bg, bmin, ORDt, ORDy, c,
                          '' if not below else '; !! %d TERMS BELOW THE PHYSICAL POLES' % below))
            bp = bernoulli(pd, sm, bg, N)
            write_expr(bp, os.path.join(d, 'seriesBer%s.m' % pt_),
                       '%s, colour layer %s, at %s,%s -> 0 (corner %s): Bernoulli variables %s -> 1 - Exp[-s%s], %s -> 1 - Exp[-s%s], '
                       'order %d in s%s and s%s; Log[%s] -> Log%s, Log[%s] -> Log%s are placeholders (substitute at the point)'
                       % (name, key, sm, bg, pt_, sm, sm, bg, bg, N, sm, bg, sm, sm.upper(), bg, bg.upper()))
            print('  %-7s %-18s %s: %5d terms, lowest powers %s^%s %s^%s, %s; Bernoulli %d terms'
                  % (name, key, pt_, len(pd), sm, amin, bg, bmin, 'spurious poles cancelled' if not below else '!! %d BELOW' % below, len(bp)),
                  flush=True)
        if still:
            print('  %d layers need more primes' % len(still), flush=True)
            target += max(4, a.jobs)
        pending = still
    for f in os.listdir(WORK):
        if f.startswith('A_c%d_' % c) or f.startswith('F_c%d_' % c):
            os.remove(os.path.join(WORK, f))
    print('=== corner %s done  (%.1fs)' % (pt_, time.time() - t0), flush=True)
