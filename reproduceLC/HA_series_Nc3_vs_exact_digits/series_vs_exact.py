#!/usr/bin/env python3
# =============================================================================
#  series_vs_exact.py -- Nc^3 layer (C3F0Z0), HA1..6: expansion around corner 1
#  (x, y) -> 0 versus the exact GPL expression, on a lattice of points in the triangle.
#
#  exact  : alpha[i,3] / gamma[i,3] of 2024_Vjet_3L_LC_ancillary/FiniteRemainders/decay.m,
#           coefficient of Nc^3 (= value at Nc = 1, Nf = 0), q2 -> 1, z = 1-x-y.
#           Identical to our HAG["C3F0Z0", i] (checked to 1e-36); its G functions are written
#           in z (fibration {z, y}), which GiNaC evaluates fast everywhere in the triangle,
#           while the {x, y} form is very slow where x -> 1-y.  All G's by GiNaC (ginsh).
#  native : HAser["C3F0Z0", i] (order N in x and y), evaluated at (x, y)
#  Bernoulli: the same series re-expanded in bx = -log(1-x), by = -log(1-y) to order N
#           (x = 1-e^-bx, y = 1-e^-by; Log[y] kept), evaluated at (bx, by)
#  digits : min(16, -log10 |S - E| / |E|), SEPARATELY for the real and the imaginary part;
#           a part with |E| < chop (1e-4) is not counted (NaN).
#
#  usage: python3 series_vs_exact.py <series_dir> <N> <decay.m> <outdir> [K=30] [ncores=2]
#  grid : x = i/K, y = j/K, i, j >= 1, i + j <= K-1   (K = 30 -> 406 points)
# =============================================================================
import os
import re
import sys
import json
import time
import subprocess
from fractions import Fraction
from multiprocessing import Pool

import mpmath as mp
import numpy as np
LD = np.longdouble

mp.mp.dps = 40
SERDIR, N, DECAY, OUT = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
K = int(sys.argv[5]) if len(sys.argv) > 5 else 30
NCORES = int(sys.argv[6]) if len(sys.argv) > 6 else 2
CHOP = 1e-4
NAMES = {1: 'alpha[1, 3]', 2: 'alpha[2, 3]', 3: 'alpha[3, 3]', 4: 'gamma[1, 3]', 5: 'gamma[2, 3]', 6: 'gamma[3, 3]'}
GRE = re.compile(r'(?<![A-Za-z])G\[[^\[\]]*\]')


def split_top(s):
    out, d, cur = [], 0, ''
    for ch in s:
        if ch in '[({': d += 1
        elif ch in '])}': d -= 1
        if ch == ',' and d == 0:
            out.append(cur); cur = ''
        else:
            cur += ch
    out.append(cur)
    return out


def pyexpr(e):
    e = e.replace('\n', ' ').replace('\\', '').replace('^', '**')
    e = e.replace('Pi', 'PI').replace('Zeta[3]', 'Z3').replace('Zeta[5]', 'Z5')
    e = re.sub(r'\bI\b', 'IU', e)
    return re.sub(r'(?<![\w.])(\d+)(?![\w.])', r'Q(\1)', e)


def parse_c(v):
    v = v.replace(' ', '')
    if not v.endswith('I'):
        return mp.mpc(mp.mpf(v), 0)
    body = v[:-1].rstrip('*'); cut = -1
    for k in range(len(body) - 1, 0, -1):
        if body[k] in '+-' and body[k - 1] not in 'eE':
            cut = k; break
    re_, im_ = ('0', body) if cut < 0 else (body[:cut], body[cut:])
    if im_ in ('', '+'): im_ = '1'
    if im_ == '-': im_ = '-1'
    return mp.mpc(mp.mpf(re_), mp.mpf(im_))


BASE = {'Q': mp.mpf, 'PI': mp.pi, 'Z3': mp.zeta(3), 'Z5': mp.zeta(5), 'IU': mp.mpc(0, 1)}

# ------------------------------------------------------------------ exact (published form)
s = re.sub(r'\(\*.*?\*\)', '', open(DECAY).read(), flags=re.S).strip()[1:-1]
pub = {}
for r in split_top(s):
    lhs, rhs = r.split('->', 1)
    pub[lhs.strip()] = rhs.strip()
EX_ATOMS = sorted(set(a for i in NAMES for a in GRE.findall(pub[NAMES[i]])), key=len, reverse=True)
EX_CODE = {}
for i, nm in NAMES.items():
    e = pub[nm]
    for k, g in enumerate(EX_ATOMS):
        e = e.replace(g, ' GG%d ' % k)
    EX_CODE[i] = compile(pyexpr(e), nm, 'eval')


# our own GPL form (fibration {x, y}): fast where x is small compared with 1-y
GDIR = os.environ.get('GDIR', '/mnt/user-data/outputs/HA_GPL_leading')
MY = {}
for i in NAMES:
    b = re.sub(r'\(\*.*?\*\)', '', open(os.path.join(GDIR, 'C3F0Z0_HA%d_G.m' % i)).read(), flags=re.S).split('=', 1)[1].strip()
    MY[i] = b[1:].rstrip()[:-2] if b.startswith('(') else b.rstrip(';')
MY_ATOMS = sorted(set(a for i in NAMES for a in GRE.findall(MY[i])), key=len, reverse=True)
MY_CODE = {}
for i in NAMES:
    e = MY[i]
    for k, g in enumerate(MY_ATOMS):
        e = e.replace(g, ' GG%d ' % k)
    MY_CODE[i] = compile(pyexpr(e), 'HAG%d' % i, 'eval')


def ginsh_values(atoms, xq, yq, digits=30, timeout=None):
    lines = ['Digits=%d;' % digits, 'x=%s:' % xq, 'y=%s:' % yq, 'z=1-x-y:']
    for g in atoms:
        a = split_top(g[2:-1])
        lines.append('evalf(G({%s},%s));' % (','.join(t.strip() for t in a[:-1]), a[-1].strip()))
    r = subprocess.run(['ginsh'], input='\n'.join(lines) + '\n', capture_output=True, text=True, timeout=timeout)
    out = [o.strip() for o in r.stdout.split('\n') if o.strip()]
    out = out[-len(atoms):]
    assert len(out) == len(atoms) and not r.stderr.strip(), (r.stderr[:300], len(out))
    return [parse_c(v) for v in out]


sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gpl_eval as ge


def _word(g):
    a = [t.strip() for t in split_top(g[2:-1])]
    return tuple(a[:-1]), a[-1]


MY_WORDS = [_word(g) for g in MY_ATOMS]


def ld2mp(v):
    return mp.mpf(np.format_float_scientific(v, precision=21))


def exact(xq, yq):
    """our GPL form HAG["C3F0Z0", i] (identical to the published alpha/gamma[i,3] Nc^3 part,
    checked to 1e-36), with every G evaluated by gpl_eval.py (longdouble, ~1e-19)."""
    X, Y = Fraction(xq), Fraction(yq)
    xv, yv = mp.mpf(X.numerator) / X.denominator, mp.mpf(Y.numerator) / Y.denominator
    xl = LD(X.numerator) / X.denominator; yl = LD(Y.numerator) / Y.denominator
    gx = ge.gpl_values([w for w, v in MY_WORDS if v == 'x'], xl, {'0': LD(0), '1': LD(1), '1 - y': 1 - yl, '-y': -yl})
    gy = ge.gpl_values([w for w, v in MY_WORDS if v == 'y'], yl, {'0': LD(0), '1': LD(1)})
    env = dict(BASE, x=xv, y=yv)
    for k, (w, v) in enumerate(MY_WORDS):
        env['GG%d' % k] = ld2mp((gx if v == 'x' else gy)[w])
    return [eval(MY_CODE[i], env) for i in range(1, 7)], 'xy-gpl_eval'


def exact_ginsh(xq, yq, digits=30):
    """independent cross-check: published form, all G's by GiNaC"""
    X, Y = Fraction(xq), Fraction(yq)
    xv, yv = mp.mpf(X.numerator) / X.denominator, mp.mpf(Y.numerator) / Y.denominator
    env = dict(BASE, x=xv, y=yv, z=1 - xv - yv, q2=mp.mpf(1), Nc=mp.mpf(1), Nf=mp.mpf(0))
    gv = ginsh_values(EX_ATOMS, xq, yq, digits)
    env.update({'GG%d' % k: v for k, v in enumerate(gv)})
    return [eval(EX_CODE[i], env) for i in range(1, 7)]


LD = np.longdouble


# ------------------------------------------------------------------ series tables
def read_series(i):
    f = os.path.join(SERDIR, 'C3F0Z0_HA%d_series%d.m' % (i, N))
    b = re.sub(r'\(\*.*?\*\)', '', open(f).read(), flags=re.S).split('=', 1)[1].strip()
    b = b[1:].rstrip()[:-2] if b.startswith('(') else b.rstrip(';')
    tab = {}
    for t in b.split('\n + '):
        t = t.strip()
        m = n = kx = ky = 0
        fac = []
        for p in t.split('*'):
            if p == 'x': m = 1
            elif p.startswith('x^'): m = int(p[2:])
            elif p == 'y': n = 1
            elif p.startswith('y^'): n = int(p[2:])
            elif p == 'Log[x]': kx = 1
            elif p.startswith('Log[x]^'): kx = int(p[7:])
            elif p == 'Log[y]': ky = 1
            elif p.startswith('Log[y]^'): ky = int(p[7:])
            else: fac.append(p)
        assert kx == 0, 'Log[x] present'
        c = eval(pyexpr('*'.join(fac)), dict(BASE))
        tab[(m, n, ky)] = tab.get((m, n, ky), 0) + c
    return tab


def u_powers(nmax, order):
    """(1 - e^-b)^p as power series in b to b^order, p = -1..nmax (Laurent for p = -1)."""
    e = [Fraction((-1) ** k, 1) / Fraction(__import__('math').factorial(k)) for k in range(order + 3)]
    u = [Fraction(0)] + [-e[k] for k in range(1, order + 3)]      # 1 - e^-b = b - b^2/2 + ...
    mul = lambda a, b: [sum(a[i] * b[k - i] for i in range(k + 1)) for k in range(order + 2)]
    P = {0: {0: Fraction(1)}}
    cur = [Fraction(1)] + [Fraction(0)] * (order + 1)
    for p in range(1, nmax + 1):
        cur = mul(cur, u)
        P[p] = {k: c for k, c in enumerate(cur) if c != 0 and k <= order}
    # 1/u = (1/b) * 1/(u/b);  u/b = 1 - b/2 + b^2/6 - ...
    w = u[1:] + [Fraction(0)]
    inv = [Fraction(1)]
    for k in range(1, order + 2):
        inv.append(-sum(w[j] * inv[k - j] for j in range(1, k + 1)))
    P[-1] = {k - 1: c for k, c in enumerate(inv) if c != 0 and k - 1 <= order}
    return P


def bernoulli(tab):
    P = u_powers(N, N)
    out = {}
    for (m, n, ky), c in tab.items():
        for j, a in P[m].items():
            for k, b in P[n].items():
                if j <= N and k <= N:
                    key = (j, k, ky)
                    out[key] = out.get(key, 0) + c * mp.mpf(a.numerator) / a.denominator * mp.mpf(b.numerator) / b.denominator
    return out


def evaluate(tab, u, v, L):
    pu, pv, pl = {}, {}, {}
    tot = mp.mpc(0)
    for (m, n, ky), c in tab.items():
        if m not in pu: pu[m] = u ** m
        if n not in pv: pv[n] = v ** n
        if ky not in pl: pl[ky] = L ** ky
        tot += c * pu[m] * pv[n] * pl[ky]
    return tot


NAT = BER = None


def init():
    global NAT, BER
    NAT = {i: read_series(i) for i in range(1, 7)}
    BER = {i: bernoulli(NAT[i]) for i in range(1, 7)}


def digits(s, e):
    if abs(e) < CHOP:
        return float('nan')
    d = abs(s - e)
    return 16.0 if d < 1e-16 * abs(e) else float(min(16, -mp.log10(d / abs(e))))


def work(pt):
    xq, yq = pt
    t0 = time.time()
    E, form = exact(xq, yq)
    xv, yv = mp.mpf(Fraction(xq).numerator) / Fraction(xq).denominator, mp.mpf(Fraction(yq).numerator) / Fraction(yq).denominator
    L = mp.log(yv)
    row = {'x': float(xv), 'y': float(yv), 'exact': [], 'native': [], 'bernoulli': []}
    for i in range(1, 7):
        Sn = evaluate(NAT[i], xv, yv, L)
        Sb = evaluate(BER[i], -mp.log(1 - xv), -mp.log(1 - yv), L)
        e = E[i - 1]
        row['exact'].append([str(mp.re(e)), str(mp.im(e))])
        row['native'].append([digits(mp.re(Sn), mp.re(e)), digits(mp.im(Sn), mp.im(e))])
        row['bernoulli'].append([digits(mp.re(Sb), mp.re(e)), digits(mp.im(Sb), mp.im(e))])
    row['t'] = time.time() - t0
    row['form'] = form
    return row


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    pts = [('%d/%d' % (i, K), '%d/%d' % (j, K)) for i in range(1, K) for j in range(1, K) if i + j <= K - 1]
    print('%d points, %d exact G atoms, order %d' % (len(pts), len(EX_ATOMS), N), flush=True)
    init()
    print('series tables: native %s, Bernoulli %s terms' % ([len(NAT[i]) for i in NAT], [len(BER[i]) for i in BER]), flush=True)
    rows = []
    T0 = time.time()
    with Pool(NCORES, initializer=init) as pool:
        for k, r in enumerate(pool.imap_unordered(work, pts)):
            rows.append(r)
            if (k + 1) % 25 == 0:
                print('  %d/%d points  %.0fs' % (k + 1, len(pts), time.time() - T0), flush=True)
    rows.sort(key=lambda r: (r['x'], r['y']))
    json.dump({'order': N, 'K': K, 'chop': CHOP, 'rows': rows}, open(os.path.join(OUT, 'digits_order%d.json' % N), 'w'))
    print('done %.0fs' % (time.time() - T0))
