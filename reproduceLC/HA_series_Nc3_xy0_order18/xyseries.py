#!/usr/bin/env python3
# =============================================================================
#  xyseries.py -- expansion of the GPL form of an HA coefficient around (x,y) = (0,0)
#
#  Same procedure as convertLC.nb:  ExpandPolyLogs[ExpandPolyLogs[expr, {x,0,N}], {y,0,N}],
#  i.e. the ITERATED expansion: first in x at fixed y (x << y), then every coefficient in y.
#  Everything is expanded: the G functions AND the rational prefactors.
#
#  Result:  sum  c * x^m Log[x]^kx * y^n Log[y]^ky * (constant)   with m <= N, n <= N
#  (negative powers of x and y are kept: they come from 1/x, 1/y and from the x-expansion
#  of 1/(x+y) and of G(..,-y,..;x), i.e. powers of x/y).
#
#  Input : <layer>_HA<i>_G.m from ii2g.py (G[..., x] with indices {0, 1, 1 - y, -y},
#          G[..., y] with indices {0, 1}, rational prefactors in x, y).
#  Exact rational arithmetic (python-flint fmpq_poly for the y-dependence).
#
#  usage: python3 xyseries.py <in_G.m> <out.m> <name> [N=15]
# =============================================================================
import re
import sys
import time
from fractions import Fraction
from math import comb, factorial

import flint
import sympy as sp

N = 15            # final order in x and in y
YW = N + 90       # working truncation in y (absolute power), must exceed N + |most negative y power|
XNEG = 0          # set from the data: largest power of 1/x in a prefactor
P = flint.fmpq_poly
Q = flint.fmpq


# ------------------------------------------------------------------ y-series: {ky: (off, poly)}
def ys_const(c, off=0, ky=0):
    if isinstance(c, Fraction):
        c = Q(c.numerator, c.denominator)
    return {ky: (off, P([c]))} if c != 0 else {}


def ys_norm(off, p):
    """drop terms beyond YW and leading zeros (keep exact)"""
    if p.is_zero():
        return None
    p = p.truncate(YW - off + 1) if YW - off + 1 > 0 else P([])
    if p.is_zero():
        return None
    # strip low zeros to raise the offset
    cs = p.coeffs()
    k = 0
    while k < len(cs) and cs[k] == 0:
        k += 1
    if k:
        p = p.right_shift(k); off += k
    return off, p


def ys_add_into(A, B, scale=None):
    for ky, (ob, pb) in B.items():
        if scale is not None:
            pb = pb * scale
        if ky in A:
            oa, pa = A[ky]
            o = min(oa, ob)
            s = pa.left_shift(oa - o) + pb.left_shift(ob - o)
            r = ys_norm(o, s)
            if r is None:
                del A[ky]
            else:
                A[ky] = r
        else:
            r = ys_norm(ob, pb)
            if r is not None:
                A[ky] = r
    return A


def ys_mul(A, B):
    out = {}
    for ka, (oa, pa) in A.items():
        for kb, (ob, pb) in B.items():
            o = oa + ob
            n = YW - o + 1
            if n <= 0:
                continue
            pr = pa.mul_low(pb, n) if hasattr(pa, 'mul_low') else (pa * pb).truncate(n)
            ys_add_into(out, {ka + kb: (o, pr)})
    return out


def ys_minpow(A):
    return min((o for o, p in A.values()), default=10 ** 9)


# ------------------------------------------------------------------ x-series: {(m, kx): ys}
def xs_add_into(A, B, scale=None):
    for k, v in B.items():
        if k in A:
            ys_add_into(A[k], v, scale)
            if not A[k]:
                del A[k]
        else:
            t = ys_add_into({}, v, scale)
            if t:
                A[k] = t
    return A


def xs_mul(A, B, mmax):
    out = {}
    for (m1, k1), v1 in A.items():
        for (m2, k2), v2 in B.items():
            if m1 + m2 > mmax:
                continue
            pr = ys_mul(v1, v2)
            if pr:
                xs_add_into(out, {(m1 + m2, k1 + k2): pr})
    return out


def xs_mul_ys(A, B):
    out = {}
    for k, v in A.items():
        pr = ys_mul(v, B)
        if pr:
            out[k] = pr
    return out


# ------------------------------------------------------------------ G functions
ONE_MINUS_Y_INV = {}       # n -> ys of (1-y)^(-n)


def inv_1my(n):
    if n not in ONE_MINUS_Y_INV:
        ONE_MINUS_Y_INV[n] = {0: (0, P([comb(n + j - 1, j) for j in range(YW + 1)]))} if n > 0 else ys_const(1)
    return ONE_MINUS_Y_INV[n]


def inv_pow(a, n):
    """a^(-n) as a y-series, a in {'1', '1 - y', '-y'}"""
    if a == '1':
        return ys_const(1)
    if a == '1 - y':
        return inv_1my(n)
    if a == '-y':
        return ys_const((-1) ** n, off=-n)
    raise ValueError(a)


def int_log(m, k):
    """Int_0^x t^(m-1) log^k t dt, m >= 1  ->  list of (power of x, log power, rational)"""
    return [(m, k - j, Fraction((-1) ** j * factorial(k), factorial(k - j) * m ** (j + 1))) for j in range(k + 1)]


GMEMO = {}


def g_series(word, xmax):
    """G(word; x) as x-series to x^xmax; indices '0', '1', '1 - y', '-y'.
    (Also used for G(...;y) with indices '0','1': then the 'x' of the series is y.)"""
    key = (word, xmax)
    if key in GMEMO:
        return GMEMO[key]
    if not word:
        res = {(0, 0): ys_const(1)}
    else:
        a, f = word[0], g_series(word[1:], xmax)
        res = {}
        if a == '0':
            for (m, k), c in f.items():
                if m == 0:
                    xs_add_into(res, {(0, k + 1): c}, Q(1, k + 1))
                else:
                    for (mm, kk, r) in int_log(m, k):
                        xs_add_into(res, {(mm, kk): c}, Q(r.numerator, r.denominator))
        else:
            # 1/(t-a) = - sum_p t^p a^(-(p+1)); then integrate t^q log^k t
            for (m, k), c in f.items():
                for p in range(0, xmax - m):
                    q = m + p
                    g = ys_mul(c, inv_pow(a, p + 1))
                    if not g:
                        continue
                    for (mm, kk, r) in int_log(q + 1, k):
                        if mm <= xmax:
                            xs_add_into(res, {(mm, kk): g}, -Q(r.numerator, r.denominator))
    GMEMO[key] = res
    return res


def gy_series(word):
    """G(word; y), indices '0','1', as a y-series {ky: (off, poly)}"""
    s = g_series(word, YW)
    out = {}
    for (m, k), c in s.items():
        # c is a constant ys {0:(0,P([q]))}
        q = c[0][1].coeffs()[0] if c else 0
        if q != 0:
            ys_add_into(out, {k: (m, P([q]))})
    return out


# ------------------------------------------------------------------ rational prefactors
x_, y_ = sp.symbols('x y')


def prefactor_series(expr, xmax):
    """rational function of x, y -> x-series with x^m, m <= xmax (negative m allowed)."""
    num, den = sp.fraction(sp.together(expr))
    cden, facs = sp.factor_list(den, x_, y_)
    pn = sp.Poly(sp.expand(num), x_, y_)
    shift_x = 0
    series = {}
    # numerator
    for (i, j), c in pn.terms():
        xs_add_into(series, {(i, 0): ys_const(Q(int(sp.numer(c)), int(sp.denom(c))), off=j)})
    lead = Fraction(1) / Fraction(int(sp.numer(cden)), int(sp.denom(cden)))
    post = []
    for f, e in facs:
        f = sp.expand(f); e = int(e)
        if f == x_:
            shift_x += e; continue
        if f == -x_:
            shift_x += e; lead *= (-1) ** e; continue
        if f == y_ or f == -y_:
            if f == -y_:
                lead *= (-1) ** e
            post.append(('y', e)); continue
        for sign in (1, -1):
            g = sp.expand(sign * f)
            if g == x_ - 1:
                kind = 'x-1'
            elif g == y_ - 1:
                kind = 'y-1'
            elif g == x_ + y_:
                kind = 'x+y'
            elif g == x_ + y_ - 1:
                kind = 'x+y-1'
            else:
                continue
            lead *= sign ** e
            post.append((kind, e)); break
        else:
            raise ValueError('unexpected denominator factor %s' % f)
    xm = xmax + shift_x
    for kind, e in post:
        if kind == 'y':
            fs = {(0, 0): ys_const(1, off=-e)}
        elif kind == 'x-1':        # (x-1)^-e = (-1)^e sum C(e+p-1,p) x^p
            fs = {(p, 0): ys_const((-1) ** e * comb(e + p - 1, p)) for p in range(xm + 1)}
        elif kind == 'y-1':        # (y-1)^-e = (-1)^e (1-y)^-e
            fs = {(0, 0): ys_add_into({}, inv_1my(e), Q((-1) ** e))}
        elif kind == 'x+y':        # (x+y)^-e = sum_p C(-e,p) x^p y^(-e-p)
            fs = {(p, 0): ys_const((-1) ** p * comb(e + p - 1, p), off=-e - p) for p in range(xm + 1)}
        elif kind == 'x+y-1':      # (x+y-1)^-e = sum_p C(-e,p) x^p (y-1)^(-e-p)
            fs = {}
            for p in range(xm + 1):
                n = e + p
                fs[(p, 0)] = ys_add_into({}, inv_1my(n), Q((-1) ** p * comb(e + p - 1, p) * (-1) ** n))
        series = xs_mul(series, fs, xm)
    out = {}
    for (m, k), v in series.items():
        if m - shift_x <= xmax:
            out[(m - shift_x, k)] = ys_add_into({}, v, Q(lead.numerator, lead.denominator))
    return out, shift_x


# ------------------------------------------------------------------ parsing the G files
def split_star(s):
    out, d, cur = [], 0, ''
    for ch in s:
        if ch == '[':
            d += 1
        elif ch == ']':
            d -= 1
        if ch == '*' and d == 0:
            out.append(cur); cur = ''
        else:
            cur += ch
    out.append(cur)
    return [o for o in out if o]


def parse_gfile(path):
    s = open(path).read()
    s = re.sub(r'\(\*.*?\*\)', '', s, flags=re.S)
    body = s.split('=', 1)[1].strip()
    if body.startswith('('):          # name = ( ... );
        body = body[1:].rstrip()
        assert body.endswith(');'), 'unexpected file end'
        body = body[:-2] + ';'
    terms = []
    for line in body.split('\n'):
        line = line.strip()
        if not line.startswith('+ ('):
            continue
        line = line[2:].rstrip(';')
        # coefficient = first balanced (...)
        d = 0
        for i, ch in enumerate(line):
            if ch == '(':
                d += 1
            elif ch == ')':
                d -= 1
                if d == 0:
                    break
        coef, rest = line[1:i], line[i + 1:]
        toks = split_star(rest)
        gx = gy = ()
        cst = []
        for t in toks:
            m = re.fullmatch(r'G\[(.*), ([xy])\]', t)
            if m:
                idx = tuple(a.strip() for a in m.group(1).split(','))
                if m.group(2) == 'x':
                    gx = idx
                else:
                    gy = idx
            else:
                cst.append(t)
        terms.append((coef, '*'.join(cst) or '1', gx, gy))
    return terms


def mma_coef(q):
    q = Fraction(int(q.p), int(q.q)) if not isinstance(q, Fraction) else q
    return str(q.numerator) if q.denominator == 1 else '%d/%d' % (q.numerator, q.denominator)


def main():
    src, out, name = sys.argv[1], sys.argv[2], sys.argv[3]
    global N, YW
    if len(sys.argv) > 4:
        N = int(sys.argv[4]); YW = N + 90
    t0 = time.time()
    terms = parse_gfile(src)
    # prefactors
    pref = {}
    xneg = 0
    for coef, cst, gx, gy in terms:
        if coef not in pref:
            e = sp.sympify(coef.replace('^', '**'), locals={'x': x_, 'y': y_})
            pref[coef] = prefactor_series(e, N)
            xneg = max(xneg, pref[coef][1])
    XG = N + xneg
    print('%s: %d terms, %d distinct prefactors, max 1/x^%d  (%.0fs)' % (src.split('/')[-1], len(terms), len(pref), xneg, time.time() - t0), flush=True)
    # group: (const, gx) -> sum prefactor * G(y)
    groups = {}
    gycache = {}
    for coef, cst, gx, gy in terms:
        ps = pref[coef][0]
        if gy not in gycache:
            gycache[gy] = gy_series(gy)
        acc = groups.setdefault((cst, gx), {})
        xs_add_into(acc, xs_mul_ys(ps, gycache[gy]))
    print('   %d (constant, G(x)) groups, %d G(y) atoms  (%.0fs)' % (len(groups), len(gycache), time.time() - t0), flush=True)
    total = {}
    minpow = 0
    for (cst, gx), acc in groups.items():
        gxs = g_series(gx, XG)
        minpow = min([minpow] + [ys_minpow(v) for v in gxs.values()] + [ys_minpow(v) for v in acc.values()])
        pr = xs_mul(acc, gxs, N)
        xs_add_into(total.setdefault(cst, {}), pr)
    assert YW - N >= -2 * minpow, 'working y-depth too small: min power %d' % minpow
    # write
    nterm = 0
    nmin = 10 ** 9
    lines = []
    for cst in sorted(total):
        for (m, kx) in sorted(total[cst]):
            for ky, (off, p) in sorted(total[cst][(m, kx)].items()):
                for j, c in enumerate(p.coeffs()):
                    n = off + j
                    if c == 0 or n > N:
                        continue
                    nmin = min(nmin, n)
                    fac = [mma_coef(c)]
                    if cst != '1':
                        fac.append(cst)
                    if m:
                        fac.append('x^%d' % m if m != 1 else 'x')
                    if n:
                        fac.append('y^%d' % n if n != 1 else 'y')
                    if kx:
                        fac.append('Log[x]^%d' % kx if kx != 1 else 'Log[x]')
                    if ky:
                        fac.append('Log[y]^%d' % ky if ky != 1 else 'Log[y]')
                    lines.append('*'.join(fac))
                    nterm += 1
    with open(out, 'w') as f:
        f.write('(* %s: iterated expansion around (x,y) = (0,0) of %s:\n' % (name, src.split('/')[-1]))
        f.write('   first in x to x^%d at fixed y, then every coefficient in y to y^%d\n' % (N, N))
        f.write('   (= ExpandPolyLogs[ExpandPolyLogs[expr, {x,0,%d}], {y,0,%d}]); exact coefficients.\n' % (N, N))
        f.write('   %d terms; lowest power of y: %d. *)\n' % (nterm, nmin))
        f.write('%s = (' % name + ('\n + '.join(lines) if lines else '0') + ');\n')
    print('   wrote %s: %d terms, x^<=%d, y^%d..%d, min y power in factors %d  (%.0fs)' % (
        out, nterm, N, nmin, N, minpow, time.time() - t0), flush=True)


if __name__ == '__main__':
    main()
