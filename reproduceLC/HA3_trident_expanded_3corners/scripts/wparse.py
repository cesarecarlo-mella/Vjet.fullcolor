#!/usr/bin/env python3
"""Parse C3F0Z0_HA3.wl:  HA3 = sum_k R_k(t,u) Trident[k]  +  I*Pi*( sum_k R'_k(t,u) Trident[k] ).

Returns {'re': [(key, num, den)], 'ipi': [...]} with num/den as {(i,j): Fraction} polynomials in t,u.
"""
import re
from fractions import Fraction

LINE = re.compile(r'^\(\((.*)\)/\((.*)\)\)\*Trident\[(\d+), (\d+)\]\s*\+?\s*$')
MONO = re.compile(r'([+-]?)\s*([0-9]+(?:/[0-9]+)?)?\*?((?:[tu](?:\^\d+)?\*?)*)')


def poly(s):
    """'a/b*t^2*u - c*t + d' -> {(pt, pu): Fraction}"""
    s = s.replace(' ', '')
    out = {}
    for term in re.findall(r'[+-]?[^+-]+', s):
        sg = -1 if term[0] == '-' else 1
        term = term.lstrip('+-')
        c = Fraction(1); pt = pu = 0
        for f in term.split('*'):
            if f[0].isdigit():
                c *= Fraction(f)
            else:
                b, _, p = f.partition('^'); p = int(p) if p else 1
                if b == 't': pt += p
                elif b == 'u': pu += p
                else: raise ValueError(f)
        k = (pt, pu)
        out[k] = out.get(k, 0) + sg * c
    return out


def parse(path):
    res = {'re': [], 'ipi': []}
    part = 're'
    for line in open(path):
        line = line.strip()
        if not line or line == ')':
            continue
        if line == 'I*Pi*(':
            part = 'ipi'; continue
        m = LINE.match(line)
        if not m:
            raise ValueError(line[:200])
        n, d, i, w = m.groups()
        res[part].append(((int(i), int(w)), poly(n), poly(d)))
    return res


def peval(p, t, u):
    return sum(c * t ** a * u ** b for (a, b), c in p.items())


if __name__ == '__main__':
    import sys
    r = parse(sys.argv[1])
    print({k: len(v) for k, v in r.items()})
