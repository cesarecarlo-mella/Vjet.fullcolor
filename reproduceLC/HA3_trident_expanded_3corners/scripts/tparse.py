#!/usr/bin/env python3
"""Parse trident_selected_cK_ord20.m (<|Trident[i, w] -> series, ...|>) into numpy arrays.

Each term is  coeff * u^a v^b Log[u]^la Log[v]^lb Pi^p Zeta[3]^z3 Zeta[5]^z5   (u,v = corner variables).
Output (cached .npz): for every term: fid, a, b, la, lb, p, z3, z5, num (str), den (str).
"""
import re
import sys
import time
import numpy as np

VARS = {'c1': ('x', 'y'), 'c2': ('y', 'z'), 'c3': ('x', 'z')}
KEY = re.compile(r'Trident\[(\d+), (\d+)\] -> ')


def parse_term(t, u, v):
    sgn = 1
    if t[0] == '-':
        sgn = -1; t = t[1:]
    den = 1
    if t[0] == '(':
        j = t.rindex(')/')
        den = int(t[j + 2:]); t = t[1:j]
    else:
        m = re.fullmatch(r'(.*)/(\d+)', t)
        if m:
            t, den = m.group(1), int(m.group(2))
    num = sgn
    e = [0] * 7
    for f in t.split('*'):
        if f.isdigit():
            num *= int(f); continue
        if '^' in f:
            base, p = f.split('^'); p = int(p)
        else:
            base, p = f, 1
        if base == u: e[0] += p
        elif base == v: e[1] += p
        elif base == 'Log[%s]' % u: e[2] += p
        elif base == 'Log[%s]' % v: e[3] += p
        elif base == 'Pi': e[4] += p
        elif base == 'Zeta[3]': e[5] += p
        elif base == 'Zeta[5]': e[6] += p
        else:
            raise ValueError('unknown factor %r in %r' % (f, t))
    return e, num, den


LIM = 2 ** 63


def ldq(n, d):
    if abs(n) < LIM and d < LIM:
        return np.longdouble(n) / np.longdouble(d)
    import mpmath as mp
    mp.mp.dps = 30
    return np.longdouble(mp.nstr(mp.mpf(n) / d, 25))


def parse_file(path, corner):
    """streams the file; returns keys, exps (int16, [fid,a,b,la,lb,p,z3,z5]), coef (longdouble, rational part)"""
    u, v = VARS[corner]
    s = open(path).read().strip()
    assert s.startswith('<|') and s.endswith('|>')
    keys, EX, CO = [], [], []
    pos = 2
    end = len(s) - 2
    fid = 0
    while pos < end:
        nxt = s.find(',\n ', pos)
        if nxt < 0 or nxt > end:
            nxt = end
        ent = s[pos:nxt]
        pos = nxt + 3
        m = KEY.match(ent)
        keys.append((int(m.group(1)), int(m.group(2))))
        body = re.sub(r'\s+', ' ', ent[m.end():]).strip()
        if body != '0':
            body = body.replace(' + ', ' \x00+').replace(' - ', ' \x00-')
            terms = body.split(' \x00')
            ex = np.zeros((len(terms), 8), dtype=np.int16)
            co = np.empty(len(terms), dtype=np.longdouble)
            for j, t in enumerate(terms):
                e, n, d = parse_term(t.lstrip('+'), u, v)
                ex[j, 0] = fid; ex[j, 1:] = e
                co[j] = ldq(n, d)
            EX.append(ex); CO.append(co)
        fid += 1
    return keys, np.concatenate(EX), np.concatenate(CO)


if __name__ == '__main__':
    corner = sys.argv[1]
    t0 = time.time()
    keys, ex, co = parse_file('sel_%s.m' % corner, corner)
    np.save('parsed_%s_coef.npy' % corner, co)
    np.savez('parsed_%s.npz' % corner, keys=np.array(keys), ex=ex)
    print(corner, len(keys), 'functions', len(ex), 'terms',
          'max exps', ex.max(axis=0), '%.0fs' % (time.time() - t0))
