#!/usr/bin/env python3
"""Symbol-level test: can the eps^0 remainders be written with CLASSICAL polylogs?

Criterion (Goncharov): F of weight n is a polynomial in classical Li_k (any arguments)
only if the cobracket components delta_{n-k,k}(F), 2 <= k <= n/2, vanish.
  delta_{a,b}(w) = rho(w[:a]) (x) rho(w[a:])  -  rho(w[b:]) (x) rho(w[:b])
  rho = projector killing shuffle products (Duhr-Gangl-Rhodes recursion).
Symbol of II[a1..an] = a1 (x) ... (x) an (a1 = first integration = first entry);
letters are treated as multiplicatively independent.
Pure functions: F_P = sum_terms r_term(P) * II[word] with the rational prefactor
evaluated at a random point P (mod p = 2^61-1); vanishing must hold for every P.

Usage:
  python3 classical_obstruction.py test            # per amplitude, 2 random points
  python3 classical_obstruction.py rank C W A B N  # rank of delta_{A,B} obstruction space,
        e.g. rank 1 4 2 2 20   (C = constant prefactor: 1, Pi^2, ...; N = points/amplitude;
        add 'lin' as 7th arg to keep only words with letters 1-6)
"""
import re, glob, os, sys, random, collections, time
from functools import lru_cache
P = (1 << 61) - 1
inv = lambda x: pow(x, P - 2, P)
D = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
QUAD, ALG = set(range(9, 21)), {7, 8}
by = collections.defaultdict(list)          # amp -> [(const, word, den, coef, poly)]
for f in sorted(glob.glob(os.path.join(D, '*.wl'))):
    amp = os.path.basename(f)[:-3]; L = open(f).read().splitlines()
    i = next(k for k, l in enumerate(L) if l.startswith('+eps^(0)*')); buf = []
    for l in L[i + 1:]:
        m = re.match(r'^\+\(([-0-9/]+)\)\*\((.*)\)\*\((.*)\)$', l)
        if m:
            c = m.group(3); w = re.search(r'II\[([0-9,]*)\]', c).group(1)
            buf.append((c[:c.index('II[')].rstrip('*'), tuple(int(x) for x in w.split(',')) if w else (), m.group(1), m.group(2)))
        elif l.startswith(')/('):
            by[amp] += [(c, w, l[3:-1], co, p) for c, w, co, p in buf]; buf = []

@lru_cache(None)
def rho(w):
    n = len(w)
    if n == 1: return {w: 1}
    out = collections.Counter(); f = (n - 1) * inv(n) % P
    for x, c in rho(w[:-1]).items(): out[x + (w[-1],)] += c * f
    for x, c in rho(w[1:]).items(): out[x + (w[0],)] -= c * f
    return {k: v % P for k, v in out.items() if v % P}

def delta(F, a, b):
    Z = collections.defaultdict(collections.Counter)
    for w, c in F.items(): Z[w[:a]][w[a:]] += c; Z[w[b:]][w[:b]] -= c
    out = collections.Counter()
    for Lw, R in Z.items():
        RR = collections.Counter()
        for r, c in R.items():
            for y, d in rho(r).items(): RR[y] += c * d
        RR = {y: v % P for y, v in RR.items() if v % P}
        for x, c in (rho(Lw).items() if RR else ()):
            for y, d in RR.items(): out[(x, y)] += c * d
    return {k: v % P for k, v in out.items() if v % P}

def groups(amp, const, w, lin=False):
    g = collections.defaultdict(collections.Counter)
    for c, wd, d, co, p in by[amp]:
        if c == const and len(wd) == w and (not lin or max(wd) <= 6):
            a_, _, b_ = co.partition('/'); g[(p, d)][wd] += int(a_) * inv(int(b_) if b_ else 1)
    return [(compile(p.replace('^', '**'), 'p', 'eval'), compile(d.replace('^', '**'), 'd', 'eval'), v) for (p, d), v in g.items()]

def pure(G, rnd):
    env = {'t': rnd.randrange(P), 'u': rnd.randrange(P)}; F = collections.Counter()
    for pc, dc, vec in G:
        r = eval(pc, env) % P * inv(eval(dc, env) % P) % P
        for x, c in vec.items(): F[x] += r * c
    return {x: c % P for x, c in F.items() if c % P}

def sector(k):
    s = set(k[0]) | set(k[1]); return ('quad' if s & QUAD else '') + ('+' if s & QUAD and s & ALG else '') + ('alg' if s & ALG else '') or 'lin'

TESTS = [('', 4, (2, 2)), ('', 5, (3, 2)), ('', 6, (4, 2)), ('', 6, (3, 3)), ('Pi^2', 4, (2, 2))]
if sys.argv[1:2] == ['test']:
    rnd = random.Random(7)
    for amp in sorted(by):
        res = []
        for const, w, (a, b) in TESTS:
            G = groups(amp, const, w)
            if not G: continue
            Ds = [delta(pure(G, rnd), a, b) for _ in range(2)]
            res.append(f"{const or '1'}*w{w} d{a},{b}: {len(Ds[0])}/{len(Ds[1])} " + str(dict(collections.Counter(sector(k) for k in Ds[0]))))
        print(amp, ' | '.join(res), flush=True)
elif sys.argv[1:2] == ['rank']:
    const = '' if sys.argv[2] == '1' else sys.argv[2]; w, a, b, N = map(int, sys.argv[3:7]); lin = sys.argv[7:8] == ['lin']
    basis = {}; rnd = random.Random(11)
    def add(v):
        v = dict(v)
        while v:
            piv = min(v)
            if piv not in basis: f = inv(v[piv]); basis[piv] = {k: c * f % P for k, c in v.items()}; return 1
            r = basis[piv]; f = v[piv]
            for k, c in r.items():
                x = (v.get(k, 0) - f * c) % P
                if x: v[k] = x
                else: v.pop(k, None)
        return 0
    for rnd_pass in range(2):
        new = 0
        for amp in sorted(by):
            G = groups(amp, const, w, lin)
            for _ in range(N if G else 0): new += add(delta(pure(G, rnd), a, b))
        print(f"{sys.argv[2]}*w{w} delta_{a},{b} {'lin' if lin else 'all'}: pass {rnd_pass+1} rank {len(basis)} (+{new})", flush=True)
else:
    print(__doc__)
