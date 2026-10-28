"""exact canonical DE of the generators in an integral basis (replaces the generator basis of exact_deq.py, whose
connection had rationals of ~100 digits):

    at each weight w (from 6 down), L_w = the Z-span of the generators of weight w and of the strips S_k b of the basis
    functions b of weight w + 1.  Its basis B_w is found by HNF + LLL on the coordinates at r pivot words (r = rank;
    the projection is injective on the span, so it preserves the lattice), then reconstructed exactly in word space.
    Every generator and every strip is then an integer combination of B_w: all entries of M[k] are integers.

    dF = Sum_k dlog W_k M[k] . F,   F = (B[6,1..], ..., B[1,1..], One)

    At the weights listed as redundant (option, e.g. 6,5) no basis is chosen: the functions are the distinct generating
    vectors themselves (the generators and the strips from above, up to a scale factor), linearly dependent.  The
    connection into such a weight then only holds the scale factors, and its generators are functions of the system.
    (Large integers at weights 5-6 come from relations among the generators with large coefficients; a redundant set
    never solves for them.)

writes basis.wl (label -> II[...] combination), deq.wl ({labels, {M[1], ..., M[20]}}) and generators.wl (each generator
T[g] as a combination of the functions); every relation is checked exactly, word by word.
usage: integral_basis.py <Generators.wl> <outdir> [redundant weights, e.g. 6,5]"""
import sys, time, collections
from fractions import Fraction
from math import lcm
import flint
import iisys

src, out = sys.argv[1:3]
REDUNDANT = {int(a) for a in sys.argv[3].split(',')} if len(sys.argv) > 3 else set()
P = iisys.P
t0 = time.time()
gens = iisys.parse_generators(src)
byw = collections.defaultdict(list)
for name, f in gens:
    byw[len(next(iter(f)))].append((name, f))
top = max(byw)


def pivot_words(rows):
    """r pivot words (mod p row echelon) and r independent rows: the r x r minor on them is invertible"""
    piv, indep = [], []
    for n_, f in enumerate(rows):
        r = {w: iisys.modp(c) for w, c in f.items()}
        for col, pr in piv:
            if r.get(col):
                a = r[col]
                for c2, v2 in pr.items():
                    r[c2] = (r.get(c2, 0) - a * v2) % P
        r = {c: v for c, v in r.items() if v}
        if r:
            col = min(r)
            inv = pow(r[col], P - 2, P)
            piv.append((col, {c: v * inv % P for c, v in r.items()}))
            indep.append(n_)
    return [c for c, _ in piv], indep


basis = {0: [{(): Fraction(1)}]}
names = {}                                   # w -> labels of the functions of weight w
coords = {}                                  # generator name -> {(w, i): int}
M = {}                                       # (k, (w, i), (w - 1, j)) -> int
for w in range(top, 0, -1):
    S = [('gen', name, f) for name, f in byw[w]]
    for i, b in enumerate(basis.get(w + 1, [])):
        for k in range(1, iisys.NL + 1):
            g = iisys.strip(b, k)
            if g:
                S.append(('strip', (i, k), g))
    if w in REDUNDANT:
        members, keymap, nm = [], {}, []
        for kind, key, f in S:
            k0 = min(f)
            kk = tuple(sorted((wd, c / f[k0]) for wd, c in f.items()))
            if kk not in keymap:
                keymap[kk] = len(members); members.append(f)
                nm.append(key if kind == 'gen' else 'S[%d,%d]' % (w, len(members)))
            j = keymap[kk]
            sc = f[k0] / members[j][k0]
            assert all(f[wd] == sc * c for wd, c in members[j].items()) and len(f) == len(members[j])
            if kind == 'gen':
                coords[key] = {(w, j): sc}
            else:
                i, k = key
                M[k, (w + 1, i), (w, j)] = sc
        basis[w], names[w] = members, nm
        big = max(max(abs(x.numerator), x.denominator) for (_, a, b), x in M.items() if a[0] == w + 1) if w < top else 0
        print('weight %d (redundant): %d generating vectors -> %d functions; connection from weight %d: scale factors up to %s  (%.0f s)'
              % (w, len(S), len(members), w + 1, big, time.time() - t0), flush=True)
        continue
    D = 1
    for _, _, f in S:
        for c in f.values():
            D = lcm(D, c.denominator)
    rows = [{wd: int(c * D) for wd, c in f.items()} for _, _, f in S]
    cols, indep = pivot_words(rows)
    r = len(cols)
    Pm = flint.fmpz_mat([[row.get(c, 0) for c in cols] for row in rows])
    H = Pm.hnf()
    H = flint.fmpz_mat([[H[i, j] for j in range(r)] for i in range(r)])          # the non-zero rows: a lattice basis
    R = H.lll()
    # coordinates of every generating vector in R:  Pm = C . R  (must be integral)
    Cq = flint.fmpq_mat(Pm) * flint.fmpq_mat(R).inv()
    C = [[Cq[i, j] for j in range(r)] for i in range(len(S))]
    assert all(c.q == 1 for row in C for c in row), 'coordinates not integral at weight %d' % w
    C = [[int(c.p) for c in row] for row in C]
    # the basis in word space: R = X . Pm_indep for r independent generating vectors (exact rational solve)
    Pi = flint.fmpq_mat([[rows[j].get(c, 0) for c in cols] for j in indep])
    X = flint.fmpq_mat(R) * Pi.inv()
    B = []
    for a in range(r):
        v = collections.defaultdict(Fraction)
        for b_, j in enumerate(indep):
            x = X[a, b_]
            if x != 0:
                xf = Fraction(int(x.p), int(x.q))
                for wd, c in rows[j].items():
                    v[wd] += xf * c
        B.append({wd: c / D for wd, c in v.items() if c})
    basis[w] = B
    names[w] = ['B[%d,%d]' % (w, i + 1) for i in range(r)]
    # exact checks (word by word) and bookkeeping
    for (kind, key, f), crow in zip(S, C):
        chk = dict(f)
        for j, c in enumerate(crow):
            if c:
                for wd, v in B[j].items():
                    chk[wd] = chk.get(wd, 0) - c * v
        assert not any(chk.values()), 'relation fails at weight %d for %s %s' % (w, kind, key)
        if kind == 'gen':
            coords[key] = {(w, j): c for j, c in enumerate(crow) if c}
        else:
            i, k = key
            for j, c in enumerate(crow):
                if c:
                    M[k, (w + 1, i), (w, j)] = c
    bg = max((abs(c) for (kind, _, _), row in zip(S, C) if kind == 'gen' for c in row), default=0)
    bs = max((abs(c) for (kind, _, _), row in zip(S, C) if kind == 'strip' for c in row), default=0)
    print('weight %d: %d generating vectors, rank %d, D = %d; largest entry: M (strips from weight %d) %d, generators %d  (%.0f s)'
          % (w, len(S), r, D, w + 1, bs, bg, time.time() - t0), flush=True)
# weight 0: the lattice of constants spanned by 1 and the strips of the weight-1 basis is g Z, g = gcd; B[0,1] = g
consts = [Fraction(1)] + [iisys.strip(b, k)[()] for b in basis[1] for k in range(1, iisys.NL + 1) if iisys.strip(b, k)]
den = 1
for c in consts:
    den = lcm(den, c.denominator)
num = 0
for c in consts:
    num = __import__('math').gcd(num, int(c * den))
g0 = Fraction(num, den)
basis[0] = [{(): g0}]
for i, b in enumerate(basis[1]):
    for k in range(1, iisys.NL + 1):
        st = iisys.strip(b, k)
        if st:
            q = st[()] / g0
            assert q.denominator == 1
            M[k, (1, i), (0, 0)] = int(q)
for name, f in byw.get(0, []):
    q = f[()] / g0
    assert q.denominator == 1
    coords[name] = {(0, 0): int(q)}
print('weight 0: B[0,1] = %s (the constant g; One = %s B[0,1])' % (g0, 1 / g0))

order = [(w, i) for w in range(top, -1, -1) for i in range(len(basis[w]))]
idx = {key: n for n, key in enumerate(order)}
names[0] = ['B[0,1]']
label = [names[w][i] for w, i in order]
n = len(order)
print('%d functions; %d non-zero entries, largest |numerator|, denominator of M: %d, %d'
      % (n, len(M), max(abs(Fraction(v).numerator) for v in M.values()), max(Fraction(v).denominator for v in M.values())))


def wl(x):
    x = Fraction(x)
    return str(x.numerator) if x.denominator == 1 else '%d/%d' % (x.numerator, x.denominator)


with open('%s/basis.wl' % out, 'w') as fh:
    fh.write('(* integral basis of the closed system, weight 6 -> 0 (integral_basis.py): B[w,i] -> II[...] combination;\n'
             '   B[0,1] is a constant (II[] = 1).  Every generator is an integer combination of these (generators.wl). *)\n{\n')
    fh.write(',\n'.join('%s -> %s' % (lab, '+'.join('(%s)*II[%s]' % (wl(c), ','.join(map(str, wd)))
                                                    for wd, c in sorted(basis[w][i].items())))
                        for lab, (w, i) in zip(label, order)) + '\n}\n')
with open('%s/deq.wl' % out, 'w') as fh:
    fh.write('(* dF = Sum[dlog W[k] M[k] . F, {k, 20}],  F = labels of basis.wl in that order (%d functions), W[k] the\n'
             '   letters of alphabet.wl.  {labels, {M[1], ..., M[20]}} *)\n' % n)
    fh.write('{{%s},\n{' % ', '.join(label))
    mats = []
    for k in range(1, iisys.NL + 1):
        ent = sorted((idx[a] + 1, idx[b] + 1, v) for (kk, a, b), v in M.items() if kk == k)
        mats.append('SparseArray[{%s}, {%d, %d}]' % (', '.join('{%d,%d}->%s' % (a, b, wl(v)) for a, b, v in ent), n, n))
    fh.write(',\n'.join(mats) + '}}\n')
with open('%s/generators.wl' % out, 'w') as fh:
    fh.write('(* the generators of Generators.wl as integer combinations of the basis B[w,i] of basis.wl *)\n{\n')
    fh.write(',\n'.join('%s -> %s' % (name, ' + '.join('(%s)*%s' % (wl(c), label[idx[key]]) for key, c in sorted(coords[name].items())))
                        for name, _ in gens) + '\n}\n')
print('wrote basis.wl, deq.wl, generators.wl in %s  (%.0f s)' % (out, time.time() - t0))
