"""the 84 remainders in double-double (about 32 digits): the quad counterpart of remainder_eval.remainders_from.

    remainders_from_dd(points, T, rsign=1)   T: generators as a ddnp.CD array (npts, 579), or complex double (mixed:
                                             the rational functions and multipliers in quad, the generators as given)
    generators_dd(F)                         the 579 generators from the 586 functions of the exact system (CD arrays)

The coefficients sum_m q * TranscendentalMultipliers[m] are recomputed exactly (mpmath) and rounded to double-double;
the rational functions are evaluated in double-double at the double-precision input point."""
import os, re
from fractions import Fraction
import numpy as np
import mpmath as mp
import remainder_eval as R
from ddnp import DD, CD, dd_sum

HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 40


def _tm_mp(s):
    """multiplier string -> its constant (mp), the power of r being remainder_eval.RPOW"""
    e = s.replace('Zeta[3]', 'z3').replace('Zeta[5]', 'z5').replace('^', '**')
    return mp.mpc(eval(e, dict(Pi=mp.pi, I=mp.mpc(0, 1), z3=mp.zeta(3), z5=mp.zeta(5), r=mp.mpf(1))))


_TM = [_tm_mp(s) for s in R.TM_STR]
_ROWS = []                                                  # per key: (r, g, k index arrays, CD coefficients)
for n in R.NAMES:
    src = open(os.path.join(R.ZG, 'FiniteRemainders', n + '.wl')).read()
    body = src[src.index('*)') + 2:]
    for key, rows in re.findall(r'"([^"]+)"\s*->\s*\{(.*?)\}\s*(?=,\s*"|\s*\}\s*$)', body, re.S):
        acc = {}
        for r, g, m, q in re.findall(r'\{(\d+),(\d+),(\d+),(-?\d+(?:/\d+)?)\}', rows):
            k = (int(r) - 1, int(g) - 1, R.RPOW[int(m) - 1])
            q = Fraction(q)
            acc[k] = acc.get(k, 0) + mp.mpf(q.numerator) / q.denominator * _TM[int(m) - 1]
        assert R.KEYS[len(_ROWS)] == (n, key)
        _ROWS.append(tuple(np.array([k[i] for k in acc]) for i in range(3)) + (CD.from_mp(list(acc.values())),))
_KPOW = sorted({int(k) for row in _ROWS for k in row[2]})


def _rational_dd(points, rsign, region='decay'):
    P = np.asarray(points, float)
    x, y = DD(P[:, 0]), DD(P[:, 1])
    if region != 'decay':                                   # the physical point in double-double from (X, Y)
        X, Y = x, y
        Z = 1.0 - X - Y
        x, y = {'z': (-X / Z, -Y / Z), 'x': (1.0 / Z, -X / Z), 'y': (-X / Z, 1.0 / Z)}[region]
        rsign = -1
    t, u = y, 1.0 - x - y
    one = DD(np.ones(len(P)))
    RF = [one * f(t, u) for f in R._RFB]
    rr = (x * y * u).sqrt()
    if rsign < 0:
        rr = -rr
    return RF, {k: rr ** k for k in _KPOW}


def _stack(L):
    return DD(np.stack([a.hi for a in L], -1), np.stack([a.lo for a in L], -1))


def remainders_from_dd(points, T, rsign=1, chunk=64, region='decay'):
    """-> CD array (npts, 84).  region: points are (X, Y) of a production region (physical point and r = -Sqrt[x y z]
    computed here)"""
    if not isinstance(T, CD):
        T = CD(DD(np.real(T)), DD(np.imag(T)))
    out_re, out_im = [], []
    for c0 in range(0, len(points), chunk):
        sl = slice(c0, c0 + chunk)
        RF, RK = _rational_dd(points[sl], rsign, region)
        RFs, RKs = _stack(RF), _stack([RK[k] for k in _KPOW])   # (n, 129), (n, powers of r)
        Ts = T[sl]
        cols = []
        for r, g, k, c in _ROWS:
            pref = RFs[:, r] * RKs[:, np.searchsorted(_KPOW, k)]   # rational function times r^k: (n, rows)
            cols.append(dd_sum((Ts[:, g] * c) * pref, axis=1))
        out_re.append(_stack([v.re for v in cols])); out_im.append(_stack([v.im for v in cols]))
    cat = lambda L: DD(np.concatenate([a.hi for a in L]), np.concatenate([a.lo for a in L]))
    return CD(cat(out_re), cat(out_im))


_GEN = None


def generators_dd(F):
    """F: CD array (npts, 586) -> CD array (npts, 579)"""
    global _GEN
    if _GEN is None:
        import quad_data as Q
        Q.generators([mp.mpc(0)] * Q.NF)                         # parses out/generators.wl
        _GEN = []
        for _, comb in Q._GEN:
            c = np.array([float(a) for a, _ in comb]); i = np.array([b for _, b in comb])
            assert all(abs(a) < 2 ** 26 and a == int(a) for a, _ in comb)
            _GEN.append((i, c))
    cols = [dd_sum(F[:, i] * c[None, :], axis=1) for i, c in _GEN]
    return CD(_stack([v.re for v in cols]), _stack([v.im for v in cols]))
