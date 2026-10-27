#!/usr/bin/env python3
"""Exact (coefficient-by-coefficient) check at corner 1 to order 18:

    HA3 = sum_k R_k(t,u) T_k  +  I Pi sum_k R'_k(t,u) T_k        (C3F0Z0_HA3.wl, t = y, u = 1-x-y)

The R_k have the spurious denominator (x+y)^3 (= (1-u)^3), so  sum R_k T_k  cannot be re-expanded
as a double series directly.  Multiply by the common denominator
    D = t^2 u (u-1)^3 (t-1) (t+u-1) (t+u)^3 :
    D * HA3 = sum_k P_k(x,y) T_k(x,y),   P_k = D R_k polynomials.
The coefficient of x^a y^b (any logs/constants) of the right-hand side only needs the T_k
coefficients with powers <= (a, b), so it is exact for a, b <= 20 with the order-20 trident series.
The left-hand side uses HAser["C3F0Z0",3] (order 18, lowest power y^-1, from the GPL form);
D * HAser is exact for a, b <= 18.  Both sides are compared for all a, b <= 18, separately for
every monomial Log[x]^i Log[y]^j Pi^p Zeta3^q Zeta5^r I^s, modulo two 31-bit primes.
"""
import re
import sys
import time
from fractions import Fraction

import numpy as np
import sympy as sp

from wparse import parse
from tparse import KEY, parse_term

PRIMES = (2147483647, 2147483629)
N = 18
SER = sys.argv[1] if len(sys.argv) > 1 else '/mnt/user-data/outputs/HA_series_Nc3_xy0_order18/C3F0Z0_HA3_series18.m'
t0 = time.time()
x, y, t, u = sp.symbols('x y t u')

# ---------------------------------------------------------------- P_k = D R_k in (x, y)
D_tu = t**2 * u * (u - 1)**3 * (t - 1) * (t + u - 1) * (t + u)**3
D_xy = sp.Poly(sp.expand(D_tu.subs({t: y, u: 1 - x - y})), x, y)
wl = parse('C3F0Z0_HA3.wl')
Pk = []          # (key, part, {(a,b): Fraction})
DTU = sp.Poly(D_tu, t, u, domain='QQ')
for part, lst in wl.items():
    for key, num, den in lst:
        pn = sp.Poly(sum(sp.Rational(c.numerator, c.denominator) * t**a * u**b for (a, b), c in num.items()), t, u, domain='QQ')
        pd = sp.Poly(sum(sp.Rational(c.numerator, c.denominator) * t**a * u**b for (a, b), c in den.items()), t, u, domain='QQ')
        q, r = (DTU * pn).div(pd)
        assert r.is_zero, ('not divisible', key)
        pxy = sp.Poly(q.as_expr().subs({t: y, u: 1 - x - y}), x, y, domain='QQ')
        Pk.append((key, part, {m: Fraction(int(c.p), int(c.q)) for m, c in pxy.terms()}))
print('P_k built: %d (%.0fs); deg D = %s' % (len(Pk), time.time() - t0, D_xy.total_degree()), flush=True)


def modq(q, p):
    return (q.numerator % p) * pow(q.denominator % p, -1, p) % p


# ---------------------------------------------------------------- T_k (corner 1), exact mod p
s = open('sel_c1.m').read().strip()
Tk = {}            # key -> {struct: [ (a, b, n, d) ]}
pos, end = 2, len(s) - 2
need = set(k for k, _, _ in Pk)
while pos < end:
    nxt = s.find(',\n ', pos)
    if nxt < 0 or nxt > end:
        nxt = end
    ent = s[pos:nxt]; pos = nxt + 3
    m = KEY.match(ent)
    key = (int(m.group(1)), int(m.group(2)))
    if key not in need:
        continue
    body = re.sub(r'\s+', ' ', ent[m.end():]).strip()
    d = {}
    if body != '0':
        for tt in body.replace(' + ', ' \x00+').replace(' - ', ' \x00-').split(' \x00'):
            e, n, dd = parse_term(tt.lstrip('+'), 'x', 'y')
            if e[0] > N or e[1] > N:
                continue
            d.setdefault(tuple(e[2:]), []).append((e[0], e[1], n, dd))
    Tk[key] = d
del s
print('T_k read (%.0fs)' % (time.time() - t0), flush=True)

results = {}
for p in PRIMES:
    # RHS: struct (lx, ly, pi, z3, z5, I) -> (N+1)x(N+1) matrix mod p
    rhs = {}
    for key, part, P in Pk:
        extra_pi, I = (1, 1) if part == 'ipi' else (0, 0)
        for st, terms in Tk[key].items():
            T = np.zeros((N + 1, N + 1), dtype=np.int64)
            for a, b, n, dd in terms:
                T[a, b] = (T[a, b] + modq(Fraction(n, dd), p)) % p
            full = (st[0], st[1], st[2] + extra_pi, st[3], st[4], I)
            acc = rhs.setdefault(full, np.zeros((N + 1, N + 1), dtype=np.int64))
            for (da, db), c in P.items():
                if da > N or db > N:
                    continue
                cm = modq(c, p)
                acc[da:, db:] = (acc[da:, db:] + cm * T[:N + 1 - da, :N + 1 - db]) % p
    # LHS: D * HAser
    b = re.sub(r'\(\*.*?\*\)', '', open(SER).read(), flags=re.S).split('=', 1)[1].strip()
    b = b[1:].rstrip()[:-2]
    ser = {}
    for term in b.split('\n + '):
        st = [0, 0, 0, 0, 0, 0]; a = bb = 0; c = Fraction(1)
        for f in term.strip().split('*'):
            if f[0] == '-' and not f[1].isdigit():
                c = -c; f = f[1:]
            base, _, e = f.partition('^'); e = int(e) if e else 1
            if base == 'x': a = e
            elif base == 'y': bb = e
            elif base == 'Log[x]': st[0] = e
            elif base == 'Log[y]': st[1] = e
            elif base == 'Pi': st[2] = e
            elif base == 'Zeta[3]': st[3] = e
            elif base == 'Zeta[5]': st[4] = e
            elif base == 'I': st[5] = 1
            else: c *= Fraction(f)
        ser.setdefault(tuple(st), {})[(a, bb)] = ser.get(tuple(st), {}).get((a, bb), 0) + c
    lhs = {}
    for st, cs in ser.items():
        acc = lhs.setdefault(st, np.zeros((N + 1, N + 1), dtype=np.int64))
        for (a, bb), c in cs.items():
            for (da, db), dc in D_xy.terms():
                A, B = a + da, bb + db
                if 0 <= A <= N and 0 <= B <= N:
                    acc[A, B] = (acc[A, B] + modq(c * Fraction(int(dc.p), int(dc.q)), p)) % p
                elif B < 0 and A <= N:
                    raise ValueError('negative power survives in D*HAser')
    keys = set(rhs) | set(lhs)
    z = np.zeros((N + 1, N + 1), dtype=np.int64)
    bad = {k: int(np.count_nonzero((rhs.get(k, z) - lhs.get(k, z)) % p)) for k in keys}
    nz = sum(int(np.count_nonzero(lhs.get(k, z))) for k in keys)
    results[p] = (len(keys), nz, {k: v for k, v in bad.items() if v})
    print('prime %d: %d log/constant structures, %d nonzero coefficients of D*HA3 (a,b<=%d); mismatches: %s  (%.0fs)'
          % (p, len(keys), nz, N, results[p][2] or 'none', time.time() - t0), flush=True)
    lx = [k for k in rhs if k[0] > 0 and np.count_nonzero(rhs[k] % p)]
    print('   Log[x] structures with nonzero RHS coefficients:', lx or 'none (Log[x] cancels)')
