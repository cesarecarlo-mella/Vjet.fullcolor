"""high-precision evaluation of the iterated-integral system (cpp/iter_quad: double-double, about 32 digits):

  double_series(dist):             the 586 functions of the exact system (out/basis.wl) as one double series around the
        corner (0,0),  sum c x^(a/2) y^(b/2) log^k(x) log^l(y),  with exact coefficients: the blow-up of the ITER series
        is undone (chart 1: x = t y^2, y = B; chart 4: y = t x^2, x = B).  No negative powers appear, and the two charts
        agree exactly where both have a term, so the union of the two truncations is used: with orders t 16 / B 48 it
        holds every term of total degree <= 32 (the first missing one is x^16 y^17), error ~ (|x| + |y|)^33.
  series_values(ser, point):       the functions at a point with 0 < x, y (256-bit fixed point; converted to mp)
  write_data(fn, point, values, r0):   the integer system (out/deq.wl) with the start point and the start values
  run(points, data, mode, tol, delta):  integrate from the start of the data file to the points -> [586 mp values]
  path_margin(S, P, d):            as gen_eval.path_margin, for the path from S to P
  generators(F):                   the 579 generators (out/generators.wl) from the 586 functions
"""
import os, re, struct, subprocess, pickle, sys
from collections import defaultdict
from fractions import Fraction
from math import comb
import numpy as np
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, 'cpp', 'iter_quad')
INPUTS = os.path.join(HERE, '..', '..', 'blowup_series', 'families', 'ITER', 'inputs')
mp.mp.dps = 40
_LAB = os.path.join(HERE, 'out', 'labels.m')                # copy of ITER inputs/labels.m (shipped)
_labels = re.findall(r'"([^"]+)"', open(_LAB if os.path.exists(_LAB) else os.path.join(INPUTS, 'labels.m')).read())
NF = len(_labels)
LIDX = {l: i for i, l in enumerate(_labels)}
PREC = 256                                          # fixed-point bits of the series evaluation


def _const(ck):
    v = mp.mpc(1)
    for sym, e2 in ck:
        v *= {'Pi': mp.pi, 'Zeta[3]': mp.zeta(3), 'Zeta[5]': mp.zeta(5), 'I': mp.mpc(0, 1)}[sym] ** (mp.mpf(e2) / 2)
    return v


def _reversed(dist, chart, w):
    """the terms of rec_c<chart>_w<w> in (x, y): {(function, 2 pow x, 2 pow y, log x, log y, constants): Fraction}"""
    sys.path.insert(0, os.path.join(HERE, '..', '..', 'blowup_series', 'code'))
    import pipe_common as pc
    out = defaultdict(Fraction)
    for (tag, i, it, iy), q in pc.load_pickle(os.path.join(dist, 'rec_c%d_w%d.pkl' % (chart, w))).items():
        pt, py, la, lb, ck = tag
        p2 = 2 * (it - 1) + pt                                      # small variable (x in chart 1): t^p B^b = s^p o^(b - 2p)
        o2 = 2 * (iy - 1) + py - 2 * p2
        for j in range(la + 1):                                     # log t = L(s) - 2 L(o), log B = L(o)
            c = q * comb(la, j) * (-2) ** (la - j)
            out[(i, p2, o2, j, la - j + lb, ck) if chart == 1 else (i, o2, p2, la - j + lb, j, ck)] += c
    return out


def double_series(dist, cache=None):
    """the double series from charts 1 and 4 of dist (cached in out/double_series_<dist>.pkl)"""
    cache = cache or os.path.join(HERE, 'out', 'double_series_%s.pkl' % os.path.basename(os.path.normpath(dist)))
    if os.path.exists(cache):
        return pickle.load(open(cache, 'rb'))
    terms = {}                                                      # constants -> list of (i, a2, b2, k, l, coefficient)
    sol0 = [Fraction(s.strip()) for s in open(os.path.join(INPUTS, 'sol0.m')).read().strip()[1:-1].split(',')]
    for i, v in enumerate(sol0):
        if v:
            terms.setdefault((), []).append((i, 0, 0, 0, 0, v))
    for w in range(1, 7):
        A = _reversed(dist, 1, w)
        amax = max(k[1] for k in A); bmax = max(k[2] + 2 * k[1] for k in A)    # the truncation of chart 1 (doubled)
        for k, v in _reversed(dist, 4, w).items():
            if not (k[1] <= amax and k[2] + 2 * k[1] <= bmax):     # chart 4 adds the terms chart 1 has not
                A[k] += v
        for (i, a2, b2, kx, ky, ck), v in A.items():
            if v:
                terms.setdefault(ck, []).append((i, a2, b2, kx, ky, v))
        print('double series: weight %d done' % w, flush=True)
    ser = {'amax': amax, 'bmax': bmax, 'groups': {}}
    for ck, L in terms.items():                                     # store compactly: index arrays + fixed-point coefficients
        idx = np.array([t[:5] for t in L], dtype=np.int32)
        ser['groups'][ck] = (idx, [(t[5].numerator << PREC) // t[5].denominator for t in L])
    pickle.dump(ser, open(cache, 'wb'), protocol=4)
    return ser


def chart_series(dist, chart, cache=None):
    """the exact series of the ITER family in one blow-up chart (SMALL = t BIG^2), not reversed: the same structure as
    double_series, in the chart variables (t, B) (cached in out/chart_series_<dist>_c<chart>.pkl).  Needed at the
    corners c2, c3, where the functions are not double series in the two small variables."""
    cache = cache or os.path.join(HERE, 'out', 'chart_series_%s_c%d.pkl' % (os.path.basename(os.path.normpath(dist)), chart))
    if os.path.exists(cache):
        return pickle.load(open(cache, 'rb'))
    sys.path.insert(0, os.path.join(HERE, '..', '..', 'blowup_series', 'code'))
    import pipe_common as pc
    terms = {}
    sol0 = [Fraction(s.strip()) for s in open(os.path.join(INPUTS, 'sol0.m')).read().strip()[1:-1].split(',')]
    for i, v in enumerate(sol0):
        if v:
            terms.setdefault((), []).append((i, 0, 0, 0, 0, v))
    for w in range(1, 7):
        for (tag, i, it, iy), q in pc.load_pickle(os.path.join(dist, 'rec_c%d_w%d.pkl' % (chart, w))).items():
            pt, py, la, lb, ck = tag
            if q:
                terms.setdefault(ck, []).append((i, 2 * (it - 1) + pt, 2 * (iy - 1) + py, la, lb, q))
    ser = {'chart': chart, 'groups': {}}
    for ck, L in terms.items():
        idx = np.array([t[:5] for t in L], dtype=np.int32)
        ser['groups'][ck] = (idx, [(t[5].numerator << PREC) // t[5].denominator for t in L])
    pickle.dump(ser, open(cache, 'wb'), protocol=4)
    return ser


def series_values(ser, point, maxdeg=None, phases=(0, 0)):
    """the functions from a series in two variables (u, w) (double_series: (x, y); chart_series: (t, B)) at
    |u| = point[0], |w| = point[1] -> list of NF mp complex.
    phases: (theta_u, theta_w): log u = log|u| + i theta_u, u^(p/2) = |u|^(p/2) exp(i theta_u p/2), the same for w
    (0: the decay; the continuation into a production region: theta = +pi for a variable that crosses zero from
    above, and log t -> log|t| - i pi, log B -> log|B| + i pi in a chart SMALL = t BIG^2).
    maxdeg: only the terms with (a + b)/2 <= maxdeg (to see the convergence)"""
    with mp.workprec(PREC + 64):
        u, w = mp.mpf(point[0]), mp.mpf(point[1])
        assert u > 0 and w > 0
        th = [mp.mpf(phases[0]), mp.mpf(phases[1])]
        fpc = lambda v: (int(mp.nint(mp.re(v) * mp.mpf(2) ** PREC)), int(mp.nint(mp.im(v) * mp.mpf(2) ** PREC)))
        allidx = np.concatenate([g[0] for g in ser['groups'].values()])
        na, nb, nk, nl = (int(allidx[:, j].max()) + 1 for j in (1, 2, 3, 4))
        lu, lw = mp.mpc(mp.log(u), th[0]), mp.mpc(mp.log(w), th[1])
        UL = [[fpc(u ** (mp.mpf(a) / 2) * mp.expj(th[0] * a / 2) * lu ** k) for k in range(nk)] for a in range(na)]
        WL = [[fpc(w ** (mp.mpf(b) / 2) * mp.expj(th[1] * b / 2) * lw ** l) for l in range(nl)] for b in range(nb)]
        real = th[0] == 0 and th[1] == 0
        F = [mp.mpc(0)] * NF
        for ck, (idx, cf) in ser['groups'].items():
            Sr, Si = [0] * NF, [0] * NF
            for (i, a2, b2, k, l), c in zip(idx.tolist(), cf):
                if maxdeg is not None and a2 + b2 > 2 * maxdeg:
                    continue
                ur, ui = UL[a2][k]; wr, wi = WL[b2][l]
                if real:
                    Sr[i] += (((c * ur) >> PREC) * wr) >> PREC
                else:
                    pr, pi = (c * ur) >> PREC, (c * ui) >> PREC
                    Sr[i] += (pr * wr - pi * wi) >> PREC
                    Si[i] += (pr * wi + pi * wr) >> PREC
            cv = _const(ck)
            F = [f + cv * mp.mpc(sr, si) / mp.mpf(2) ** PREC if (sr or si) else f for f, sr, si in zip(F, Sr, Si)]
    return [+f for f in F]


def _dd(v):
    hi = float(v)
    return hi, float(v - mp.mpf(hi))


def write_data(fn, point, values, r0=None):
    """the integer system of out/deq.wl, start point (x, y), start values (NF mp complex), root r0 (default +Sqrt[xyz])"""
    src = open(os.path.join(HERE, 'out', 'deq.wl')).read()
    src = re.sub(r'\(\*.*?\*\)', '', src, flags=re.S)
    mats = re.findall(r'SparseArray\[\{(.*?)\}, \{(\d+), \d+\}\]', src, re.S)
    assert len(mats) == 20 and int(mats[0][1]) == NF
    x, y = mp.mpf(point[0]), mp.mpf(point[1])
    if r0 is None:
        r0 = mp.sqrt(x * y * (1 - x - y))
    r0 = mp.mpc(r0)
    with open(fn, 'wb') as f:
        f.write(b'ITERQ001' + struct.pack('<2i', NF, 20))
        for body, _ in mats:
            ent = [(int(i) - 1, int(j) - 1, float(Fraction(v))) for i, j, v in re.findall(r'\{(\d+),(\d+)\}->(-?\d+(?:/\d+)?)', body)]
            f.write(struct.pack('<i', len(ent)))
            if ent:
                f.write(struct.pack('<%di' % len(ent), *[e[0] for e in ent]))
                f.write(struct.pack('<%di' % len(ent), *[e[1] for e in ent]))
                f.write(struct.pack('<%dd' % len(ent), *[e[2] for e in ent]))
        f.write(struct.pack('<4d', *_dd(x), *_dd(y)))
        f.write(struct.pack('<4d', *_dd(r0.real), *_dd(r0.imag)))
        for v in values:
            v = mp.mpc(v)
            f.write(struct.pack('<4d', *_dd(v.real), *_dd(v.imag)))


def run(points, data, mode='quad', tol=None, delta=(0.1, 0.2), kmax=None):
    """-> [(info, [NF mp complex])] for the points: (x, y) as doubles, or as mp numbers (end point in double-double);
    info['r'] is the root Sqrt[x y z] at the end of the path"""
    tol = tol if tol is not None else (1e-28 if mode == 'quad' else 1e-13)
    def fmt(v):
        if isinstance(v, (mp.mpf, mp.mpc)):
            hi, lo = _dd(mp.re(v))
            return '%.17g %.17g' % (hi, lo)
        return '%.17g 0' % float(v)
    inp = '\n'.join(fmt(p[0]) + ' ' + fmt(p[1]) for p in points)
    args = [BIN, data, mode, repr(float(tol)), repr(float(delta[0])), repr(float(delta[1]))] + ([str(kmax)] if kmax else [])
    r = subprocess.run(args, input=inp, capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError('iter_quad: ' + r.stderr[-300:])
    out, cur = [], None
    for line in r.stdout.split('\n'):
        if line.startswith('#'):
            h = line.split(); cur = []
            out.append((dict(x=float(h[1]), y=float(h[2]), steps=int(h[3]), evals=int(h[4]), time=float(h[5]),
                             r=complex(float(h[6]), float(h[7])) if len(h) > 7 else None), cur))
        elif line:
            a, b, c, d = (mp.mpf(float(s)) for s in line.split())     # each part is a double: round to 53 bits first
            cur.append(mp.mpc(a + b, c + d))
    return out


def path_margin(S, P, d, region='decay'):
    """smallest distance, in the path parameter, of a complex zero of a letter from the open path S -> P"""
    import gen_eval
    ts = np.linspace(0, 1, 5)
    X = S[0] + (ts + 4j * d[0] * ts * (1 - ts)) * (P[0] - S[0])
    Y = S[1] + (ts + 4j * d[1] * ts * (1 - ts)) * (P[1] - S[1])
    best = (np.inf, None)
    for k, f in gen_eval.profile(region)['sing'].items():
        c = np.polyfit(ts, f(X, Y), 4)
        c = c[np.argmax(np.abs(c) > 1e-14 * np.abs(c).max()):]
        for r in (np.roots(c) if len(c) > 1 else []):
            if abs(r - 1) < 1e-9 or abs(r) < 1e-9:
                continue
            dist = abs(r.imag) if 0 <= r.real <= 1 else min(abs(r), abs(r - 1))
            best = min(best, (dist, k), key=lambda b: b[0])
    return best


_GEN = None


def generators(F):
    """the generators T[1..579] (in the order of Generators.wl) from the functions"""
    global _GEN
    if _GEN is None:
        _GEN = []
        for line in open(os.path.join(HERE, 'out', 'generators.wl')):
            m = re.match(r'(T\[\d+\]) -> (.*?),?\s*$', line.strip())
            if m:
                _GEN.append((m.group(1), [(mp.mpf(Fraction(c).numerator) / Fraction(c).denominator, LIDX[l])
                                          for c, l in re.findall(r'\((-?[\d/]+)\)\*([A-Z]\[[\d,]+\])', m.group(2))]))
        _GEN.sort(key=lambda g: int(g[0][2:-1]))
    return [mp.fsum(c * F[i] for c, i in comb_) for _, comb_ in _GEN]
