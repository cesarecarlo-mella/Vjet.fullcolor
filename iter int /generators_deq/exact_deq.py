"""exact canonical DE of the generators: F = (weight 6, ..., weight 1, 1), at each weight the generators first, then the
extra functions (strips) needed to close the system;  dF = sum_k dlog W_k M[k] . F  with rational M[k].
Writes basis.wl (the functions F[i] in II[...]) and deq.wl (M[k] as SparseArrays).   usage: exact_deq.py system.pkl outdir"""
import pickle, sys
from fractions import Fraction
import flint
import iisys

S = pickle.load(open(sys.argv[1], 'rb'))
out = sys.argv[2]
names, wt = S['names'], list(S['wt'])
basis = S['basis']
P = iisys.P

# labels: generators keep their names; extras F[w, j]
order, label = [], []
for w in range(6, -1, -1):
    gens = [n for n, x in zip(names, wt) if x == w]
    for j, f in enumerate(basis.get(w, [])):
        order.append((w, f))
        label.append(gens[j] if j < len(gens) else ('One' if w == 0 else 'F[%d,%d]' % (w, j + 1 - len(gens))))
index = {}
for i, (w, _) in enumerate(order):
    index.setdefault(w, []).append(i)
print('basis sizes:', {w: len(v) for w, v in sorted(index.items(), reverse=True)}, 'total', len(order))


def pivots(vecs):
    """echelon pivot words (mod p) of the vectors: the square submatrix on them is invertible"""
    piv = []
    for f in vecs:
        r = {word: iisys.modp(c) for word, c in f.items()}
        for col, pr in piv:
            if r.get(col):
                a = r[col]
                for c2, v2 in pr.items():
                    r[c2] = (r.get(c2, 0) - a * v2) % P
        r = {c: v for c, v in r.items() if v}
        col = min(r)
        inv = pow(r[col], P - 2, P)
        piv.append((col, {c: v * inv % P for c, v in r.items()}))
    return [c for c, _ in piv]


M = {}                                          # (k, i, j) -> Fraction
for w in range(6, 0, -1):
    low = [order[i][1] for i in index[w - 1]]
    pw = pivots(low)
    d = len(low)
    B = flint.fmpq_mat(d, d, [flint.fmpq(low[j].get(pw[r], 0).numerator, low[j].get(pw[r], 0).denominator)
                              for r in range(d) for j in range(d)])
    Binv = B.inv()
    for i in index[w]:
        f = order[i][1]
        for k in range(1, iisys.NL + 1):
            g = iisys.strip(f, k)
            if not g:
                continue
            rhs = flint.fmpq_mat(d, 1, [flint.fmpq(g.get(p, Fraction(0)).numerator, g.get(p, Fraction(0)).denominator) for p in pw])
            c = Binv * rhs
            cf = [Fraction(int(c[j, 0].p), int(c[j, 0].q)) for j in range(d)]
            # exact check on every word
            chk = dict(g)
            for j, x in enumerate(cf):
                if x:
                    for word, v in low[j].items():
                        chk[word] = chk.get(word, 0) - x * v
            assert not any(chk.values()), 'not closed at weight %d' % w
            for j, x in enumerate(cf):
                if x:
                    M[k, i, index[w - 1][j]] = x
    print('weight %d -> %d: done' % (w, w - 1), flush=True)

n = len(order)
hgt = max(max(abs(x.numerator), x.denominator) for x in M.values())
print('%d non-zero entries over the 20 letters; largest numerator/denominator %d' % (len(M), hgt))


def wl(x):
    return str(x.numerator) if x.denominator == 1 else '%d/%d' % (x.numerator, x.denominator)


with open('%s/basis.wl' % out, 'w') as fh:
    fh.write('(* F = basis of the closed system, weight 6 -> 0; generators first at each weight, F[w,j] = extra\n'
             '   functions (combinations of strips of the generators), One = II[] = 1.  List of label -> function. *)\n{\n')
    rows = []
    for (w, f), lab in zip(order, label):
        terms = '+'.join('(%s)*II[%s]' % (wl(c), ','.join(map(str, word))) for word, c in sorted(f.items()))
        rows.append('%s -> %s' % (lab, terms))
    fh.write(',\n'.join(rows) + '\n}\n')
with open('%s/deq.wl' % out, 'w') as fh:
    fh.write('(* dF = Sum[dlog W[k] M[k] . F, {k, 20}],  F = labels of basis.wl in that order (%d functions),\n'
             '   W[k] the letters of 2610.anc/alphabet.wl.  {labels, {M[1], ..., M[20]}} *)\n' % n)
    fh.write('{{%s},\n{' % ', '.join(label))
    mats = []
    for k in range(1, iisys.NL + 1):
        ent = ['{%d,%d}->%s' % (i + 1, j + 1, wl(x)) for (kk, i, j), x in sorted(M.items()) if kk == k]
        mats.append('SparseArray[{%s}, {%d, %d}]' % (', '.join(ent), n, n))
    fh.write(',\n'.join(mats) + '}}\n')
print('wrote %s/basis.wl and %s/deq.wl' % (out, out))
