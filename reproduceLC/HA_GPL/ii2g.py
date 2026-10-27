#!/usr/bin/env python3
# =============================================================================
#  ii2g.py -- iterated integrals II[...] -> GPLs in (x, y), as in convertLC.nb
#
#  Input : <layer>_HA<i>.wl from remainder_denominator_bracket/:
#          sum of blocks  (0 + (q)*(poly(u,t))*(const*II[a1,...,an]) + ...)/(den(u,t))
#  Output: the same function written with G functions of x and y.
#
#  Conventions (identical to convertLC.nb):
#    * letters l[k] from all.ints/a_tildes/letters.dict.m (variables x, y);
#    * dII[a1,...,an] = dlog l[an] II[a1,...,a(n-1)]   (a1 innermost);
#    * II = path-ordered integral from the base point (0,0) along
#      (0,0) -> (0,y) -> (x,y), with the tangential (shuffle) regularisation
#      log x, log y -> 0 at the base point: this is what the notebook's
#      recursive GIntegrate does (y-leg at x = 0 with G(...;0) -> 0, then
#      x-leg at fixed y);
#    * rational prefactors: u -> 1-x-y, t -> y   (as in the notebook).
#
#  Chen's formula for the composed path:
#     II[a1..an] = sum_k  G(c_x(an),...,c_x(a_k+1); x) * G(c_y(ak),...,c_y(a1); y)
#  with c_y(l) / c_x(l) the GPL index of dlog l restricted to the y-leg (x=0) /
#  x-leg (y fixed); a letter whose restriction vanishes kills the term.
#
#  Only letters 1..6 are supported (C1F2Z0 and the other layers built from
#  them); the script stops on any other letter.
#
#  usage: python3 ii2g.py <input.wl> <output.m> [name]
# =============================================================================
import re
import sys
import time
from collections import defaultdict

import sympy as sp

x, y, u, t = sp.symbols('x y u t')

# restriction of dlog l[k] to the two legs: list of GPL indices (as strings)
#   l1 = Log[x], l2 = Log[y], l3 = Log[1-x-y], l4 = Log[1-x], l5 = Log[1-y], l6 = Log[x+y]
YLEG = {1: None, 2: '0', 3: '1', 4: None, 5: '1', 6: '0'}           # at x = 0, variable y
XLEG = {1: '0', 2: None, 3: '1 - y', 4: '1', 5: None, 6: '-y'}      # at fixed y, variable x

LINE = re.compile(r'^\+\(([-0-9/]+)\)\*\((.*)\)\*\((.*)\)$')


def gstr(idx, var):
    return 'G[%s, %s]' % (', '.join(idx), var) if idx else ''


def chen(word):
    """II[word] -> list of (Gx string, Gy string); '' stands for 1."""
    out = []
    n = len(word)
    for k in range(n + 1):
        yw, xw = word[:k], word[k:]
        if any(YLEG[a] is None for a in yw) or any(XLEG[a] is None for a in xw):
            continue
        out.append((gstr([XLEG[a] for a in reversed(xw)], 'x'),
                    gstr([YLEG[a] for a in reversed(yw)], 'y')))
    return out


def to_sympy(s):
    s = s.replace('^', '**')
    return sp.sympify(s, locals={'u': u, 't': t})


def convert(path):
    sub = {u: 1 - x - y, t: y}
    acc = defaultdict(lambda: defaultdict(lambda: sp.Integer(0)))   # (const, Gx, Gy) -> den -> numerator
    block = []
    nterms = 0
    for line in open(path):
        line = line.strip()
        if not line or line in ('(0', '+(0', ')'):
            continue
        if line.startswith(')/('):
            den = sp.factor(to_sympy(line[2:]).subs(sub))
            for q, poly, cst, word in block:
                num = sp.Rational(q) * to_sympy(poly).subs(sub)
                for gx, gy in chen(word):
                    acc[(cst, gx, gy)][den] += num
            block = []
            continue
        m = LINE.match(line)
        if not m:
            raise ValueError('cannot parse: ' + line)
        q, poly, X = m.groups()
        mw = re.fullmatch(r'(?:(.*)\*)?II\[([0-9, ]*)\]', X)
        if not mw:
            raise ValueError('cannot parse factor: ' + X)
        cst = mw.group(1) or '1'
        word = tuple(int(a) for a in mw.group(2).split(',') if a.strip())
        bad = [a for a in word if a not in XLEG]
        if bad:
            raise ValueError('letter(s) %s not supported (only 1..6): %s' % (bad, X))
        block.append((q, poly, cst, word))
        nterms += 1
    return acc, nterms


def mma(e):
    return sp.printing.mathematica.mathematica_code(e)


def write(acc, out, name, src):
    keys = sorted(acc, key=lambda k: (k[0], len(k[1]), k[1], len(k[2]), k[2]))
    nonzero = 0
    with open(out, 'w') as f:
        f.write('(* %s: GPL form of %s, built with ii2g.py\n' % (name, src))
        f.write('   variables x, y (u -> 1-x-y, t -> y); II base point (0,0), path (0,0)->(0,y)->(x,y),\n')
        f.write('   tangential regularisation; G[a1,...,an, z] = Int_0^z dt/(t-a1) G[a2,...,an, t]. *)\n')
        f.write('%s = (0' % name)
        for k in keys:
            cst, gx, gy = k
            coef = 0
            for den, num in acc[k].items():
                num = sp.expand(num)
                if num != 0:
                    coef += num / den
            if coef == 0:
                continue
            coef = sp.factor(sp.together(coef))
            if coef == 0:
                continue
            nonzero += 1
            fac = [c for c in (cst if cst != '1' else '', gx, gy) if c]
            f.write('\n + (%s)%s' % (mma(coef), ''.join('*' + c for c in fac)))
        f.write(');\n')
    return nonzero


if __name__ == '__main__':
    src, out = sys.argv[1], sys.argv[2]
    name = sys.argv[3] if len(sys.argv) > 3 else 'HAG'
    t0 = time.time()
    acc, nterms = convert(src)
    n = write(acc, out, name, src.split('/')[-1])
    print('%s: %d II terms -> %d G-monomials  (%.1fs)' % (src.split('/')[-1], nterms, n, time.time() - t0))
