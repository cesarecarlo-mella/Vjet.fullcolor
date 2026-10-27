#!/usr/bin/env python3
"""Full expansion (Trident functions AND rational prefactors) of C3F0Z0 HA3 around corner cK.

    HA3 = sum_k R_k(t,u) T_k + I Pi sum_k R'_k(t,u) T_k          (C3F0Z0_HA3.wl; t = y, u = 1-x-y)

Corner variables (p, q): c1 (x, y), c2 (y, z), c3 (x, z),  z = 1-x-y.
The R_k contain spurious poles that do not factorise into monomials ((x+y)^3 at c1, (y+z)^3 at c2,
(x+z) at c3), so the product cannot be expanded term by term.  Instead (exact arithmetic, FLINT):
  1. D = t^2 u (u-1)^3 (t-1) (t+u-1) (t+u)^3 (common denominator) is factorised in (p, q) as
     D = M * B * U, with M a monomial, B the homogeneous spurious factor, U a unit (U(0) != 0);
  2. Q = sum_k (D R_k) T_k: polynomials times the trident series -> exact for total degree <= 20
     (the trident series are exact in the rectangle p^a q^b, a, b <= 20);
  3. W = Q / U (series inverse), exact to total degree 20;
  4. degree by degree: G_n = W_{n + deg B} / B (exact division of homogeneous polynomials -- the
     remainder must vanish: this is the check that the spurious pole cancels), H = G / M.
  H is then exact for total degree n <= 20 - deg M - deg B: this is the order of the output.
Every coefficient is kept separately for each structure Log[p]^i Log[q]^j Pi^k Zeta3^l Zeta5^m I^s.

usage: python3 expand_corner.py <corner> [outdir]
"""
import re
import sys
import time
from fractions import Fraction

import flint

from wparse import parse
from tparse import KEY, parse_term

NT = 20                       # trident series: rectangle a, b <= 20  ->  complete in total degree <= 20
corner = sys.argv[1]
OUTDIR = sys.argv[2] if len(sys.argv) > 2 else '.'
VN = {'c1': ('x', 'y'), 'c2': ('y', 'z'), 'c3': ('x', 'z')}[corner]
t0 = time.time()

ctx = flint.fmpq_mpoly_ctx.get(('p', 'q'), 'lex')
P_, Q_ = ctx.gens()
ctu = flint.fmpq_mpoly_ctx.get(('t', 'u'), 'lex')
T_, U_ = ctu.gens()
one = ctx.from_dict({(0, 0): 1})
# (t, u) in terms of the corner variables (p, q)
TU = {'c1': (Q_, one - P_ - Q_),          # x = p, y = q
      'c2': (P_, Q_),                     # y = p, z = q
      'c3': (one - P_ - Q_, Q_)}[corner]  # x = p, z = q,  y = 1 - p - q


def trunc(f, n):
    return ctx.from_dict({m: c for m, c in f.to_dict().items() if m[0] + m[1] <= n}) if not f.is_zero() else f


def tu_poly(d):
    return ctu.from_dict({m: flint.fmpq(c.numerator, c.denominator) for m, c in d.items()})


# ---------------------------------------------------------------- D = M B U
D_tu = T_ ** 2 * U_ * (U_ - 1) ** 3 * (T_ - 1) * (T_ + U_ - 1) * (T_ + U_) ** 3
D = D_tu.compose(*TU, ctx=ctx)
cst, facs = D.factor()
M, B, U = one, one, one * cst
for f, e in facs:
    dct = f.to_dict()
    if (0, 0) in dct:
        U *= f ** e
    elif len(dct) == 1:
        M *= f ** e
    else:
        assert len(set(a + b for a, b in dct)) == 1, ('non-homogeneous spurious factor', f)
        B *= f ** e
degM, degB = M.total_degree(), B.total_degree()
NH = NT - degM - degB
mono = list(M.to_dict())[0]
print('%s: M = %s, B = %s, U(0) = %s  ->  exact to total degree %d  (%.0fs)'
      % (corner, M, B, U.to_dict().get((0, 0)), NH, time.time() - t0), flush=True)

# ---------------------------------------------------------------- P_k = D R_k  (polynomials in p, q)
wl = parse('C3F0Z0_HA3.wl')
Pk = []
for part, lst in wl.items():
    for key, num, den in lst:
        qq, rr = divmod(D_tu * tu_poly(num), tu_poly(den))
        assert rr.is_zero()
        Pk.append((key, part, qq.compose(*TU, ctx=ctx)))
need = set(k for k, _, _ in Pk)
print('P_k: %d  (%.0fs)' % (len(Pk), time.time() - t0), flush=True)

# ---------------------------------------------------------------- trident series (exact, total degree <= NT)
s = open('sel_%s.m' % corner).read().strip()
Tk = {}
pos, end = 2, len(s) - 2
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
            e, n, dd = parse_term(tt.lstrip('+'), VN[0], VN[1])
            if e[0] + e[1] > NT:
                continue
            d.setdefault(tuple(e[2:]), {})[(e[0], e[1])] = flint.fmpq(n, dd)
    Tk[key] = {st: ctx.from_dict(c) for st, c in d.items()}
del s
print('trident series read (%.0fs)' % (time.time() - t0), flush=True)

# ---------------------------------------------------------------- Q = sum P_k T_k
Q = {}
for key, part, P in Pk:
    extra = (1, 1) if part == 'ipi' else (0, 0)
    for st, T in Tk[key].items():
        full = (st[0], st[1], st[2] + extra[0], st[3], st[4], extra[1])
        Q[full] = Q.get(full, ctx.from_dict({})) + trunc(P * T, NT)
Q = {st: f for st, f in Q.items() if not f.is_zero()}
print('Q: %d structures  (%.0fs)' % (len(Q), time.time() - t0), flush=True)

# ---------------------------------------------------------------- 1/U to total degree NT
u0 = U.to_dict()[(0, 0)]
Rest = one - U * flint.fmpq(1) / u0                 # U = u0 (1 - Rest)
Uinv, pw = one, one
for k in range(1, NT + 1):
    pw = trunc(pw * Rest, NT)
    Uinv += pw
Uinv = Uinv * (flint.fmpq(1) / u0)

# ---------------------------------------------------------------- H = (Q / U) / B / M  degree by degree
H = {}
bad = 0
for st, f in Q.items():
    W = trunc(f * Uinv, NT)
    byd = {}
    for mexp, c in W.to_dict().items():
        byd.setdefault(mexp[0] + mexp[1], {})[mexp] = c
    for n in range(0, NT + 1):
        if n not in byd:
            continue
        g, r = divmod(ctx.from_dict(byd[n]), B)
        if not r.is_zero():
            bad += 1
            print('   NOT DIVISIBLE: structure %s, degree %d' % (st, n))
            continue
        for (a, b), c in g.to_dict().items():
            a2, b2 = a - mono[0], b - mono[1]
            if a2 + b2 <= NH:
                H.setdefault(st, {})[(a2, b2)] = Fraction(int(c.p), int(c.q))
print('H: %d structures, %d coefficients, divisibility failures: %d  (%.0fs)'
      % (len(H), sum(len(v) for v in H.values()), bad, time.time() - t0), flush=True)
mins = [min(a for st in H for a, b in H[st]), min(b for st in H for a, b in H[st])]
print('lowest powers: %s^%d, %s^%d' % (VN[0], mins[0], VN[1], mins[1]))

# ---------------------------------------------------------------- write Mathematica + a pickle for the numerics
import pickle
pickle.dump({'corner': corner, 'vars': VN, 'order': NH, 'H': H}, open('%s/HA3exp_%s.pkl' % (OUTDIR, corner), 'wb'))
name = 'HA3exp["C3F0Z0", 3, "%s"]' % corner


def term(st, a, b, c):
    f = []
    for sym, e in ((VN[0], a), (VN[1], b), ('Log[%s]' % VN[0], st[0]), ('Log[%s]' % VN[1], st[1]),
                   ('Pi', st[2]), ('Zeta[3]', st[3]), ('Zeta[5]', st[4]), ('I', st[5])):
        if e:
            f.append(sym if e == 1 else '%s^%d' % (sym, e))
    return '*'.join([str(c)] + f)


with open('%s/C3F0Z0_HA3_expanded_%s.m' % (OUTDIR, corner), 'w') as fo:
    fo.write('(* %s: full expansion (Trident functions and rational prefactors) of C3F0Z0_HA3.wl\n'
             '   around corner %s, (%s, %s) -> 0, z = 1-x-y; exact coefficients.\n'
             '   Truncation: total degree a + b <= %d in %s^a %s^b (all terms of lower total degree are exact;\n'
             '   this is the highest order fixed by the order-20 trident series, see README).\n'
             '   Lowest powers: %s^%d, %s^%d.  Built with expand_corner.py. *)\n'
             % (name, corner, VN[0], VN[1], NH, VN[0], VN[1], VN[0], mins[0], VN[1], mins[1]))
    fo.write('%s = (0' % name)
    for st in sorted(H):
        for (a, b) in sorted(H[st], key=lambda m: (m[0] + m[1], m)):
            c = H[st][(a, b)]
            if c:
                fo.write('\n + %s' % term(st, a, b, c))
    fo.write(');\n')
print('written (%.0fs)' % (time.time() - t0))
