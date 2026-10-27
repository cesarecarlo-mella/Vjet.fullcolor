#!/usr/bin/env python3
"""
build_letters.py -- dlog series of the 20 letters in the blow-up charts of charts.py
(exact rationals), at any order.

Replaces all_families/pipeline_master/00_build_letters.wl (Mathematica) and
writes the same files the pipeline reads (<tag> = c1..c6, or a continued chart such as RZ_yv):

    dlog_dt_<tag>.m : { l[k] -> d/dt   Log[l_k(x(t,B), y(t,B))] }
    dlog_dy_<tag>.m : { l[k] -> d/dB   Log[l_k(x(t,B), y(t,B))] }   (B = the chart's BIG variable)
    info_<tag>.m

Charts (z = 1 - x - y), see charts.py:

    c1 : x = t*y^2 ,  BIG = y ,  z = 1 - y - t y^2
    c2 : z = t*y^2 ,  BIG = y ,  x = 1 - y - t y^2
    c3 : x = t*z^2 ,  BIG = z ,  y = 1 - z - t z^2
  and their mirrors (SMALL <-> BIG), used for the crossed families:
    c4 : y = t*x^2 ,  BIG = x        (mirror of c1)
    c5 : y = t*z^2 ,  BIG = z        (mirror of c2)
    c6 : z = t*x^2 ,  BIG = x        (mirror of c3)
  continued regions A > 1 (v = 1/A, bh = -b/A): R<A>_<b>v (bh = t v^2) and R<A>_v<b> (v = t bh^2).

  plain letter : arg = t^a B^b Q with Q(0,0) != 0 (asserted: the blow-up
                 resolves the letter), dlog = monomial poles + Q'/Q
  l7, l8       : f = Log[(a - I r)/(a + I r)], r = Sqrt[x y z],
                 a = x z (l7), a = x y (l8):
                 d f = I (2 a' b - a b') / (r (a^2 + b)),  b = x y z
                 r is the principal root at the decay corners, as in the
                 existing letter files, and -Sqrt|xyz| in the continued
                 regions (charts.odd_sign).

Truncation: t^e with e <= ORDT, B^e with e <= ORDB (e half-integer for l7,
l8), as in the pipeline.  A numeric self-check (mpmath, 40 digits) compares
every series with the derivative of the true composed letter.

    python3 build_letters.py ORDT ORDB [--corners 1,2,3,RZ_yv] [--out ../letters]
    python3 build_letters.py 15 45 --out /tmp/l --compare /path/to/old/letters   # exact check
"""
import argparse
import os
import re
import sys
import time

import mpmath as mpm
import sympy as sp
from gmpy2 import mpq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from series2 import Ser
import charts
import mparse as mp

ap = argparse.ArgumentParser()
ap.add_argument('ordt', type=int)
ap.add_argument('ordb', type=int)
ap.add_argument('--corners', default='1,2,3,4,5,6', help='chart ids of charts.py, e.g. 1,2,3 or RZ_yv,RZ_vy')
ap.add_argument('--out', default=os.path.join(HERE, '..', 'letters'))
ap.add_argument('--compare', default=None,
                help='dir with existing dlog_d{t,y}_<tag>.m: check that both agree exactly up to '
                     'the common order (terms above it are ignored on both sides)')
a = ap.parse_args()
OT, OB = a.ordt, a.ordb
XT, XB = 4, 8                      # extra orders kept internally, cut at the end

x, y, t, B = sp.symbols('x y t B')
z = 1 - x - y

# --------------------------------------------------------------- alphabet (letters.dict.m)
PLAIN = {1: x, 2: y, 3: 1 - x - y, 4: 1 - x, 5: 1 - y, 6: x + y,
         9: 1 - (1 + x) * y, 10: 1 - x * (1 + y), 11: x - y + y**2, 12: -x + x**2 + y,
         13: (x - 1)**2 - y, 14: -y + (x + y)**2, 15: -x + (x + y)**2, 16: -x + (1 - y)**2,
         17: x + x * y + y**2, 18: x**2 + y + x * y, 19: (y - 1)**2 + x * y, 20: (x - 1)**2 + x * y}
ODD = {7: x * z, 8: x * y}          # a; b = x y z
LETTERS = list(range(1, 21))

# sign of r relative to the principal root, per chart (charts.odd_sign): +1 at the decay
# corners, which reproduces the existing letter files (00_build_letters.wl) exactly (checked
# at 15/45 against all_families_array/letters, all 20 letters, c1-c3); -1 in the continued
# regions, r = -Sqrt|xyz| (Log of the two negative variables -> log|.| + I Pi).


# --------------------------------------------------------------- helpers
def monomial_split(p):
    """polynomial p(t,B) -> (a, b, Q) with p = t^a B^b Q."""
    P = sp.Poly(sp.expand(p), t, B)
    a_ = min(m[0] for m in P.monoms())
    b_ = min(m[1] for m in P.monoms())
    Q = sp.expand(P.as_expr() / (t**a_ * B**b_))
    return a_, b_, Q


def to_ser(p, NT, NB):
    P = sp.Poly(sp.expand(p), t, B)
    s = Ser(NT, NB)
    for (i, j), c in zip(P.monoms(), P.coeffs()):
        if i <= NT and j <= NB:
            c = sp.Rational(c)
            s.c[i][j] += mpq(int(c.p), int(c.q))
    return s


def to_poly(s, e2t, e2b, imag, BIGNAME, coef=1):
    """Ser * t^(e2t/2) B^(e2b/2) [* I] -> mparse Poly, truncated to t^OT, B^OB."""
    out = {}
    for i in range(s.NT + 1):
        for j in range(s.NY + 1):
            c = s.c[i][j]
            if not c:
                continue
            et, eb = 2 * i + e2t, 2 * j + e2b
            if et > 2 * OT or eb > 2 * OB:
                continue
            assert et >= -2 and eb >= -2, "pole beyond 1/t or 1/B"
            k = []
            if imag:
                k.append(('I', 2))
            if et:
                k.append(('t', et))
            if eb:
                k.append((BIGNAME, eb))
            key = tuple(sorted(k))
            out[key] = out.get(key, 0) + c * coef
    return {k: v for k, v in out.items() if v}


def padd(p, q):
    r = dict(p)
    for k, v in q.items():
        r[k] = r.get(k, 0) + v
    return {k: v for k, v in r.items() if v}


def plain_dlog(arg, BIGNAME):
    """dlog_t, dlog_B of a rational arg(t, B)."""
    num, den = sp.fraction(sp.cancel(sp.together(arg)))
    DT, DB = {}, {}
    for part, sgn in ((num, 1), (den, -1)):
        if not part.free_symbols:
            continue
        at, bb, Q = monomial_split(part)
        assert Q.subs({t: 0, B: 0}) != 0, "blow-up does not resolve this letter: Q(0,0) = 0  (%s)" % Q
        NT, NB = OT + XT, OB + XB
        iQ = to_ser(Q, NT, NB).inv()
        dt = to_poly(to_ser(sp.diff(Q, t), NT, NB) * iQ, 0, 0, False, BIGNAME, sgn)
        db = to_poly(to_ser(sp.diff(Q, B), NT, NB) * iQ, 0, 0, False, BIGNAME, sgn)
        if at:
            dt = padd(dt, {(('t', -2),): mpq(sgn * at)})
        if bb:
            db = padd(db, {((BIGNAME, -2),): mpq(sgn * bb)})
        DT, DB = padd(DT, dt), padd(DB, db)
    return DT, DB


def odd_dlog(A, Bq, BIGNAME, sign):
    """d f, f = Log[(A - I r)/(A + I r)], r = sign*Sqrt[Bq] (principal Sqrt > 0)."""
    nb, db_ = sp.fraction(sp.cancel(sp.together(Bq)))
    pn, qn, hn = monomial_split(nb)
    pd_, qd, hd = monomial_split(db_)
    c_n, c_d = hn.subs({t: 0, B: 0}), hd.subs({t: 0, B: 0})
    if c_n < 0:                                 # b > 0 with both parts negative: flip both
        hn, hd, c_n, c_d = -hn, -hd, -c_n, -c_d
    assert c_n > 0 and c_d > 0, "b does not have a positive limit shape"
    NT, NB = OT + XT, OB + XB
    sq_n, sq_d = to_ser(hn, NT, NB).sqrt(), to_ser(hd, NT, NB).sqrt()
    out = []
    for var in (t, B):
        g = sp.cancel(sp.together((2 * sp.diff(A, var) * Bq - A * sp.diff(Bq, var)) / (A**2 + Bq)))
        gn, gd = sp.fraction(g)
        if gn == 0:
            out.append({}); continue
        an, bn, Gn = monomial_split(gn)
        ad, bd, Gd = monomial_split(gd)
        assert Gd.subs({t: 0, B: 0}) != 0, "a^2 + b: non-monomial zero at the corner"
        # d f = I * g / r = I * t^e2t/2 B^e2b/2 * Gn * Sqrt[hd] / (Gd * Sqrt[hn])
        e2t = 2 * (an - ad) - (pn - pd_)
        e2b = 2 * (bn - bd) - (qn - qd)
        S = to_ser(Gn, NT, NB) * sq_d * (to_ser(Gd, NT, NB) * sq_n).inv()
        out.append(to_poly(S, e2t, e2b, True, BIGNAME, sign))
    return out


def evalpoly(p, tv, bv, BIGNAME):
    s = mpm.mpc(0)
    for k, c in p.items():
        term = mpm.mpf(c.numerator) / c.denominator
        for at, e in k:
            if at == 'I':
                term *= 1j
            elif at == 't':
                term *= tv ** (mpm.mpf(e) / 2)
            elif at == BIGNAME:
                term *= bv ** (mpm.mpf(e) / 2)
            else:
                raise ValueError(at)
        s += term
    return s


def letter_value(k, X, Y):
    """true letter l_k(X, Y) (mpmath), principal branches as in letters.dict.m."""
    if k in PLAIN:
        return mpm.log(sp.lambdify((x, y), PLAIN[k], 'mpmath')(X, Y))
    Z = 1 - X - Y
    r = mpm.sqrt(X * Z * Y)
    A = X * Z if k == 7 else X * Y
    return mpm.log((A - 1j * r) / (A + 1j * r))


def truncate(p, BIGNAME, ot, ob):
    """keep t^e B^f with e <= ot, f <= ob (the pipeline's load-time cut)."""
    out = {}
    for k, v in p.items():
        d = dict(k)
        if d.get('t', 0) <= 2 * ot and d.get(BIGNAME, 0) <= 2 * ob:
            out[k] = v
    return out


# --------------------------------------------------------------- build
os.makedirs(a.out, exist_ok=True)
mpm.mp.dps = 40
failed = False
for c in charts.parse_list(a.corners):
    TAG, (SMALL, BIGNAME), SIGN = charts.tag(c), charts.small_big(c), charts.odd_sign(c)
    X_, Y_ = charts.param(c, t, B)
    t0 = time.time()
    sub = {x: X_, y: Y_}
    print("=== %s : %s, BIG = %s, order t^%d %s^%d" % (TAG, charts.blowup(c), BIGNAME, OT, BIGNAME, OB), flush=True)
    DT, DB = {}, {}
    for k in LETTERS:
        if k in PLAIN:
            DT[k], DB[k] = plain_dlog(PLAIN[k].subs(sub, simultaneous=True), BIGNAME)
        else:
            A = ODD[k].subs(sub, simultaneous=True)
            Bq = (x * y * z).subs(sub, simultaneous=True)
            DT[k], DB[k] = odd_dlog(A, Bq, BIGNAME, SIGN)
    print("  built in %.1fs" % (time.time() - t0), flush=True)

    # ---- numeric self-check vs the true composed letter
    tv, bv = mpm.mpf(3) / 100, mpm.mpf(7) / 100
    chart = sp.lambdify((t, B), (X_, Y_), 'mpmath')
    worst = 0
    for k in LETTERS:
        L = lambda T, Bv: letter_value(k, *chart(T, Bv))
        rt = mpm.diff(lambda T: L(T, bv), tv)
        rb = mpm.diff(lambda Bv: L(tv, Bv), bv)
        if k in ODD:
            rt, rb = rt * SIGN, rb * SIGN
        for name, ser_, ref in (('dt', DT[k], rt), ('d' + BIGNAME, DB[k], rb)):
            v = evalpoly(ser_, tv, bv, BIGNAME)
            rel = abs(v - ref) / max(abs(ref), mpm.mpf(10) ** -30)
            worst = max(worst, rel)
            if rel > 1e-12:
                print("  !! l%d %s: series %s  vs  letter %s  (rel %.1e)" % (k, name, mpm.nstr(v, 15), mpm.nstr(ref, 15), rel))
    print("  numeric check at (t, %s) = (3/100, 7/100): max rel. difference %s" % (BIGNAME, mpm.nstr(worst, 3)), flush=True)
    if worst > 1e-12:
        sys.exit("numeric check FAILED")

    # ---- optional exact comparison with existing letter files
    if a.compare:
        with open(os.path.join(a.compare, 'info_%s.m' % TAG)) as fh:
            info = fh.read()
        old_ot = int(re.search(r'"ORDt" -> (\d+)', info).group(1))
        old_ob = int(re.search(r'"ORDy" -> (\d+)', info).group(1))
        o_t, o_b = min(OT, old_ot), min(OB, old_ob)          # common order
        for name, D in (('dlog_dt_%s.m' % TAG, DT), ('dlog_dy_%s.m' % TAG, DB)):
            old = {}
            for _, lhs, rhs in mp.parse_file(os.path.join(a.compare, name)):
                (lk, _), = next(iter(lhs))
                old[int(lk[2:-1])] = rhs
            bad = [k for k in LETTERS if k in old
                   and truncate(D[k], BIGNAME, o_t, o_b) != truncate(old[k], BIGNAME, o_t, o_b)]
            missing = [k for k in LETTERS if k not in old]
            print("  compare %s at %d/%d: %d letters in old file, identical: %s%s"
                  % (name, o_t, o_b, len(old), "ALL" if not bad else "NO, differ: %s" % bad,
                     ("  (not in old file: %s)" % missing) if missing else ""), flush=True)
            failed = failed or bool(bad)

    for name, D in (('dlog_dt_%s.m' % TAG, DT), ('dlog_dy_%s.m' % TAG, DB)):
        with open(os.path.join(a.out, name), 'w') as fh:
            fh.write("{" + ",\n ".join("l[%d] -> %s" % (k, mp.fmt_poly(D[k]) if D[k] else "0")
                                       for k in LETTERS) + "}\n")
    with open(os.path.join(a.out, 'info_%s.m' % TAG), 'w') as fh:
        fh.write('<|"blowup" -> "%s", "corner" -> %s, "W" -> 2, "BIG" -> %s, "SMALL" -> %s, '
                 '"ORDt" -> %d, "ORDy" -> %d, "active" -> {%s}, "oddSign" -> <|7 -> %d, 8 -> %d|>, '
                 '"maxReldiff" -> %s, "builder" -> "build_letters.py"|>\n'
                 % (charts.blowup(c), c if isinstance(c, int) else '"%s"' % c, BIGNAME, SMALL, OT, OB,
                    ", ".join(map(str, LETTERS)), SIGN, SIGN, mpm.nstr(worst, 3)))
    print("  wrote dlog_dt_%s.m, dlog_dy_%s.m, info_%s.m -> %s" % (TAG, TAG, TAG, os.path.abspath(a.out)), flush=True)

sys.exit(1 if failed else 0)
