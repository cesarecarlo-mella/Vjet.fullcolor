#!/usr/bin/env python3
"""
he_corner.py -- coverage of the production-region triangle by ONE corner expansion, tested on LC
finite remainders (2024_Vjet_3L_LC_ancillary, region y > 1).

Triangle: v = 1/y, u = -z/y, w = -x/y = 1 - u - v  (u, v, w > 0).  Corners:
  th  (threshold, v -> 1, u, w -> 0): continued decay corner; series in (x, y) = (w, u);
  he2 (high energy, w -> 1, u, v -> 0): y -> oo at fixed z = -u/v; series in (x, y) = (u, v);
  he3 (high energy, u -> 1, w, v -> 0): y -> oo at fixed x = -w/v; series in (x, y) = (w, v).
  he2b, he3b: the same high-energy corners in the other ordering (v first), i.e. the square chart
  around s = 1 (resp. s~ = 1): series in (x, y) = (v, u), (v, w).
The series file is the xyseries.py output (iterated: x << y first) of the coefficient written in the
corner's fibration (to_xy.wls for he2, refibrate.wls for th and he3).

Representations of the SAME data, compared with the exact remainder on a lattice of the triangle:
  native : sum c x^m y^n Log[x]^a Log[y]^b, truncated at x^N, y^N.  At he2/he3 this is a series in
           the ratio x/y, radius 1 (branch line x + y = 0 through the corner);
  he2/he3, square chart rho = x + y, s = x/(x+y)  (x = s rho, y = (1-s) rho; the branch line goes
           to s = oo):  square (rho^N, s^N), bern_rho (beta = -log(1-rho)), bern_s
           (sigma = -log(1-s)), bern_both;  log(rho), log(s) kept exactly;
  th:      bern (bx = -log(1-x), by = -log(1-y)), bern_x, bern_y;  Log[x], Log[y] kept exactly.
Each is an exact re-expansion of the series data within its own truncation.  The coverage is
always reported on the (u, v, w) triangle.

exact: the ORIGINAL G-file (fibration u, v) evaluated with gpl_eval.py (long double, ~1e-17).
Constants: Pi, Zeta exact; ii = I; q2 = 1; colour: one layer (--layer, e.g. Nc2Nf1: only the
Nc^2 Nf^1 terms, with Nc = Nf = 1) or the full sum at Nc = 3, Nf = 5.
digits = min(16, -log10|S-E|/|E|), real and imaginary part separately (a part below 1e-6 |E| is
not counted); 'err' holds the absolute error |S - E| (complex modulus) per representation.

    python3 he_corner.py <name> <G-file (u,v)> <series-file> <N> <out.json> [--corner he2] [--K 30] [--layer Nc2Nf1]
"""
import argparse
import json
import os
import re
import sys
import time
from fractions import Fraction
from math import comb, factorial

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gpl_eval
import xyseries

LD, CLD = np.longdouble, np.clongdouble
PI = LD('3.14159265358979323846264338327950288')
ZETA = {'Zeta[3]': LD('1.20205690315959428539973816151144999'),
        'Zeta[5]': LD('1.03692775514336992633136548645703417')}

ap = argparse.ArgumentParser()
ap.add_argument('name'); ap.add_argument('gfile'); ap.add_argument('sfile')
ap.add_argument('N', type=int); ap.add_argument('out')
ap.add_argument('--K', type=int, default=30)
ap.add_argument('--layer', default=None, help='e.g. Nc2Nf1: keep only the Nc^2 Nf^1 terms')
ap.add_argument('--corner', default='he2', choices=['he2', 'he3', 'th', 'he2b', 'he3b'])
ap.add_argument('--grid', default='tri', choices=['tri', 'sp'],
                help='tri: u = i/K, v = j/K inside the triangle;  sp: s = i/K, p = j/K (u = s p, v = (1-s) p)')
a = ap.parse_args()
N = a.N
if a.layer:
    m = re.fullmatch(r'Nc(\d+)Nf(\d+)', a.layer)
    LAYER = (int(m.group(1)), int(m.group(2)))


def const_value(tokens, layer=True):
    """numeric value of a list of constant tokens, or None if the term is not in the layer.
    layer=False: the file holds one colour layer already (refibrate.wls), no Nc, Nf tokens."""
    if a.layer and layer and (tokens.count('Nc'), tokens.count('Nf')) != LAYER:
        return None
    v = CLD(1)
    for t in tokens:
        if t == 'Pi':
            v *= PI
        elif t in ZETA:
            v *= ZETA[t]
        elif t in ('Nc', 'Nf'):
            v *= 1 if a.layer else (3 if t == 'Nc' else 5)
        elif t in ('ii', 'I'):
            v *= CLD(1j)
        elif t == 'Log[2]':
            v *= LD('0.693147180559945309417232121458176568')
        elif t.startswith('q2'):
            pass                                  # q2 = 1
        else:
            raise ValueError('unknown constant %r' % t)
    return v


def tof(q):
    return LD(q.numerator) / LD(q.denominator)


# ------------------------------------------------------------------ series (xyseries output)
TOK = re.compile(r'^(x|y|Log\[x\]|Log\[y\])(?:\^(-?\d+))?$')


def read_series(path):
    """-> {(m, n, a, b): complex long double}: c x^m y^n Log[x]^a Log[y]^b, constants folded in."""
    s = re.sub(r'\(\*.*?\*\)', '', open(path).read(), flags=re.S)
    has_colour = bool(re.search(r'\bN[cf]\b', s))
    body = s.split('=', 1)[1].strip()
    body = body[1:body.rindex(')')]
    out = {}
    for term in body.split('\n + '):
        term = term.strip()
        if not term or term == '0':
            continue
        toks = term.split('*')
        e = {'x': 0, 'y': 0, 'Log[x]': 0, 'Log[y]': 0}
        consts = []
        for t in toks[1:]:
            mt = TOK.match(t)
            if mt:
                e[mt.group(1)] += int(mt.group(2) or 1)
            else:
                consts.append(t)
        cv = const_value(consts, layer=has_colour)
        if cv is None:
            continue
        k = (e['x'], e['y'], e['Log[x]'], e['Log[y]'])
        out[k] = out.get(k, 0) + tof(Fraction(toks[0])) * cv
    return out


def series_mul(p, q, L):
    out = [Fraction(0)] * L
    for i, x in enumerate(p[:L]):
        if x:
            for j, y in enumerate(q[:L - i]):
                out[i + j] += x * y
    return out


def binom_series(n, L):
    """(1 - t)^n for any integer n: coefficients of t^0..t^(L-1)."""
    out, c = [], Fraction(1)
    for k in range(L):
        out.append(c)
        c = -c * (n - k) / (k + 1)
    return out


def gpow(alpha, L):
    """((1 - e^-t)/t)^alpha: coefficients of t^0..t^(L-1) (J.C.P. Miller recurrence)."""
    g = [Fraction((-1) ** k, factorial(k + 1)) for k in range(L)]
    f = [Fraction(1)]
    for m in range(1, L):
        f.append(sum(((alpha + 1) * j - m) * g[j] * f[m - j] for j in range(1, m + 1)) / m)
    return f


def add(out, key, val):
    out[key] = out.get(key, 0) + val


def to_chart(ser, rho_bern, s_bern):
    """high-energy corner: x^m y^n Lx^a Ly^b = rho^(m+n) s^m (1-s)^n (Ls + Lrho)^a (L1s + Lrho)^b,
    L1s = log(1-s), Lrho = log(rho), Ls = log(s) (the last two kept exactly).
      rho part: plain rho^p, or Bernoulli rho = 1 - e^-beta:  rho^p = beta^p ((1-e^-beta)/beta)^p;
      s part  : plain s^m (1-s)^n L1s^j expanded in s, or Bernoulli s = 1 - e^-sigma:
                s^m (1-s)^n L1s^j = sigma^m ((1-e^-sigma)/sigma)^m e^(-n sigma) (-sigma)^j.
    -> {(P, Q, A, B): c}:  X1^P X2^Q Lrho^A Ls^B, X1 = rho or beta, X2 = s or sigma; P, Q <= N."""
    out, rcache, scache = {}, {}, {}
    for (m, n, al, bl), c in ser.items():
        p = m + n
        if p > N or m > N:
            continue
        Lr, Ls = N - p + 1, N - m + 1
        if (p, Lr) not in rcache:
            rcache[(p, Lr)] = gpow(p, Lr) if rho_bern else [Fraction(1)] + [Fraction(0)] * (Lr - 1)
        R = rcache[(p, Lr)]
        for j in range(bl + 1):
            key = (m, n, j, Ls)
            if key not in scache:
                if s_bern:
                    f = series_mul(gpow(m, Ls), [Fraction((-n) ** k, factorial(k)) for k in range(Ls)], Ls)
                    scache[key] = (j, [x * (-1) ** j for x in f])
                else:
                    l1 = [Fraction(0)] + [Fraction(-1, k) for k in range(1, Ls)]
                    f = binom_series(n, Ls)
                    for _ in range(j):
                        f = series_mul(f, l1, Ls)
                    scache[key] = (0, f)
            shift, S = scache[key]
            for i in range(al + 1):
                f = c * comb(al, i) * comb(bl, j)
                A, Bl = (al - i) + (bl - j), i
                for kr, cr in enumerate(R):
                    if not cr:
                        continue
                    for ks, cs in enumerate(S):
                        Q = m + shift + ks
                        if Q > N:
                            break
                        if cs:
                            add(out, (p + kr, Q, A, Bl), f * tof(cr * cs))
    return out


def to_bern2(ser, bx, by):
    """threshold corner: x^m y^n Lx^a Ly^b with x = 1 - e^-bx and/or y = 1 - e^-by (Lx, Ly kept).
    -> {(P, Q, a, b): c}: X^P Y^Q Lx^a Ly^b, X = x or bx, Y = y or by; P, Q <= N."""
    out, gp = {}, {}
    for (m, n, al, bl), c in ser.items():
        if m > N or n > N:
            continue
        parts = []
        for e, flag in ((m, bx), (n, by)):
            L = N - e + 1
            if flag:
                if (e, L) not in gp:
                    gp[(e, L)] = gpow(e, L)
                parts.append(gp[(e, L)])
            else:
                parts.append([Fraction(1)])
        for kx, cx in enumerate(parts[0]):
            if not cx:
                continue
            for ky, cy in enumerate(parts[1]):
                if cy:
                    add(out, (m + kx, n + ky, al, bl), c * tof(cx * cy))
    return out


def evaluate(series, X1, X2, L1, L2):
    tot = np.zeros(len(X1), dtype=CLD)
    for (p, q, A, B), c in series.items():
        tot += c * X1 ** p * X2 ** q * L1 ** A * L2 ** B
    return tot


# ------------------------------------------------------------------ exact value
def exact_values(gfile, U, V):
    terms = []
    for coef, cst, gx, gy in xyseries.parse_gfile(gfile):
        cv = const_value(cst.split('*') if cst != '1' else [])
        if cv is not None:
            terms.append((coef, cv, gx, gy))
    coefs = {coef: eval(coef.replace('^', '**'), {'x': U, 'y': V}) for coef, _, _, _ in terms}
    wx = sorted({gx for _, _, gx, _ in terms if gx}); wy = sorted({gy for _, _, _, gy in terms if gy})
    out = np.zeros(len(U), dtype=CLD)
    for k in range(len(U)):
        vf = float(V[k])
        Gx = gpl_eval.gpl_values(wx, U[k], {'0': 0.0, '1': 1.0, '1 - y': 1 - vf, '-y': -vf})
        Gy = gpl_eval.gpl_values(wy, V[k], {'0': 0.0, '1': 1.0})
        Gx[()] = LD(1); Gy[()] = LD(1)
        out[k] = sum(coefs[coef][k] * cv * Gx[gx] * Gy[gy] for coef, cv, gx, gy in terms)
    return out, len(terms)


def digits(S, E):
    res = []
    for part in (np.real, np.imag):
        s, e = part(S), part(E)
        d = np.full(len(e), np.nan)
        ok = np.abs(e) > LD('1e-6') * np.abs(E)
        d[ok] = np.minimum(16, -np.log10(np.abs((s[ok] - e[ok]) / e[ok]).astype(float) + 1e-300))
        res.append(d)
    return res


t0 = time.time()
K = a.K
if a.grid == 'tri':
    pts = [(Fraction(i, K), Fraction(j, K)) for i in range(1, K) for j in range(1, K) if i + j <= K - 1]
else:
    pts = [(Fraction(i, K) * Fraction(j, K), (1 - Fraction(i, K)) * Fraction(j, K)) for i in range(1, K) for j in range(1, K)]
U = np.array([tof(u) for u, _ in pts]); V = np.array([tof(v) for _, v in pts])
W = 1 - U - V
X, Y = {'he2': (U, V), 'he3': (W, V), 'th': (W, U), 'he2b': (V, U), 'he3b': (V, W)}[a.corner]   # series variables
ser = read_series(a.sfile)
vals = {'native': evaluate(ser, X, Y, np.log(X), np.log(Y))}
if a.corner == 'th':
    BX, BY = -np.log(1 - X), -np.log(1 - Y)
    for name, bx, by, X1, X2 in (('bern', True, True, BX, BY), ('bern_x', True, False, BX, Y),
                                 ('bern_y', False, True, X, BY)):
        vals[name] = evaluate(to_bern2(ser, bx, by), X1, X2, np.log(X), np.log(Y))
else:
    rho, s = X + Y, X / (X + Y)
    beta, sigma = -np.log(1 - rho), -np.log(1 - s)
    Lr, Ls = np.log(rho), np.log(s)
    for name, rb, sb, X1, X2 in (('square', False, False, rho, s), ('bern_rho', True, False, beta, s),
                                 ('bern_s', False, True, rho, sigma), ('bern_both', True, True, beta, sigma)):
        vals[name] = evaluate(to_chart(ser, rb, sb), X1, X2, Lr, Ls)
print("%s%s corner %s: N = %d, series %d terms  (%.0fs)"
      % (a.name, ' [%s]' % a.layer if a.layer else '', a.corner, N, len(ser), time.time() - t0), flush=True)
E, nt = exact_values(a.gfile, U, V)
print("   exact: %d terms at %d points  (%.0fs)" % (nt, len(pts), time.time() - t0), flush=True)
D = {k: digits(v, E) for k, v in vals.items()}
rows = []
for k, (u, v) in enumerate(pts):
    rows.append({'u': str(u), 'v': str(v), 'exact': [float(E[k].real), float(E[k].imag)],
                 **{rep: [None if np.isnan(D[rep][0][k]) else round(float(D[rep][0][k]), 2),
                          None if np.isnan(D[rep][1][k]) else round(float(D[rep][1][k]), 2)] for rep in D},
                 'err': {rep: float(abs(vals[rep][k] - E[k])) for rep in vals}})
json.dump({'name': a.name, 'layer': a.layer, 'corner': a.corner, 'N': N, 'K': K, 'grid': a.grid, 'rows': rows}, open(a.out, 'w'))
for rep in D:
    mn = np.fmin(D[rep][0], D[rep][1])
    print("   %-9s min(Re,Im) digits over the triangle: median %4.1f | >=4: %3d%% | >=8: %3d%% | >=12: %3d%%"
          % (rep, np.nanmedian(mn), 100 * np.mean(mn >= 4), 100 * np.mean(mn >= 8), 100 * np.mean(mn >= 12)), flush=True)
