#!/usr/bin/env python3
# =============================================================================
#  gpl_eval.py -- fast numerical GPLs G(a1,...,an; X) for REAL indices that are not in (0, X]
#  (here: G(..;x) with indices {0, 1, 1-y, -y} and G(..;y) with indices {0, 1}).
#
#  G(a, w; t) = Int_0^t ds/(s - a) G(w; s), computed for every needed suffix on a common set of
#  Gauss-Legendre nodes in [0, X] by spectral cumulative integration, in numpy longdouble
#  (64-bit mantissa).  Panels:
#    * geometric towards 0 (ratio 1/3) down to t_min = 1e-40 X: resolves the log t behaviour
#      of trailing zeros and every scale of a small index such as -y;
#    * geometric towards X when an index lies just beyond X (distance d): widths ~ d/2.
#  G(0,...,0; t) = log^n t / n! exactly (tangential regularisation, as in PolyLogTools/GiNaC);
#  the part of the integral below t_min is dropped (it is O(t_min log^k t_min)).
#  Accuracy (checked against GiNaC at 30 digits): ~1e-17 relative.
# =============================================================================
import numpy as np
import mpmath as mp

LD = np.longdouble
NG = 24
_mpdps = 40


def _gl_tables(n):
    mp.mp.dps = _mpdps
    xs, ws = mp.gauss_quadrature(n, 'legendre')
    order = sorted(range(n), key=lambda k: xs[k])
    xs = [xs[k] for k in order]; ws = [ws[k] for k in order]
    # cumulative integration matrix: Int_{-1}^{x_i} f = sum_j C[i,j] f(x_j)
    # via Lagrange basis: C[i,j] = Int_{-1}^{x_i} l_j(t) dt (exact with high-order GL quadrature)
    C = mp.matrix(n, n)
    for j in range(n):
        def lj(t, j=j):
            p = mp.mpf(1)
            for m in range(n):
                if m != j:
                    p *= (t - xs[m]) / (xs[j] - xs[m])
            return p
        for i in range(n):
            # integrate polynomial of degree n-1 exactly: map GL nodes onto [-1, x_i]
            a, b = mp.mpf(-1), xs[i]
            s = 0
            for xk, wk in zip(xs, ws):
                s += wk * lj((b - a) / 2 * xk + (a + b) / 2)
            C[i, j] = s * (b - a) / 2
    X = np.array([LD(mp.nstr(v, 36)) for v in xs])
    W = np.array([LD(mp.nstr(v, 36)) for v in ws])
    Cm = np.array([[LD(mp.nstr(C[i, j], 36)) for j in range(n)] for i in range(n)])
    return X, W, Cm


GLX, GLW, GLC = _gl_tables(NG)


def panels(X, sing):
    """panel edges on [t_min, X]; sing: real indices (singular points) of the functions."""
    X = float(X)
    d = min([s - X for s in sing if s > X] + [np.inf])
    up = [X / 2, X]
    w = d / 2
    while X - w > X / 2:              # geometric towards X, smallest width d/2
        up.append(X - w); w *= 2
    low, t = [], X / 2
    while t > X * 1e-40:              # geometric towards 0, ratio 1/3
        t /= 3; low.append(t)
    return np.array([LD(e) for e in sorted(set(low + up))])


class GPLGrid:
    """all G(w; t) for t on the nodes in [0, X]; value at X is the last cumulative entry."""

    def __init__(self, X, sing):
        self.X = LD(X)
        e = panels(X, sing)
        e[-1] = LD(X)
        self.a, self.b = e[:-1], e[1:]
        h = (self.b - self.a) / 2
        self.h = h
        self.t = ((self.a + self.b) / 2)[:, None] + h[:, None] * GLX[None, :]   # (P, NG)
        self.memo = {}

    def cumint(self, f):
        """F(t) = Int_{t_min}^{t} f on the nodes, and F(X)"""
        loc = (f @ GLC.T) * self.h[:, None]               # within-panel cumulative
        tot = (f @ GLW) * self.h                          # panel totals
        start = np.concatenate(([LD(0)], np.cumsum(tot)[:-1]))
        return loc + start[:, None], np.sum(tot)

    def G(self, word, idx):
        """word: tuple of index labels; idx: label -> real value. returns (values on nodes, value at X)"""
        if word in self.memo:
            return self.memo[word]
        if not word:
            r = (np.ones_like(self.t), LD(1))
        elif all(idx[w] == 0 for w in word):
            n = len(word)
            from math import factorial
            L = np.log(self.t)
            r = (L ** n / factorial(n), np.log(self.X) ** n / factorial(n))
        else:
            a = LD(idx[word[0]])
            f, _ = self.G(word[1:], idx)
            r = self.cumint(f / (self.t - a))
        self.memo[word] = r
        return r


def gpl_values(words, X, idx):
    """{word: G(word; X)} for a list of words (tuples of labels); idx: label -> float value"""
    sing = sorted(set(float(v) for v in idx.values()))
    g = GPLGrid(X, sing)
    return {w: g.G(w, idx)[1] for w in words}
