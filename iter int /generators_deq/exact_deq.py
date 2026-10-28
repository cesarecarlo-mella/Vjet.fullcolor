"""exact canonical DE of the generators: F = (weight 6, ..., weight 1, 1), at each weight the linearly independent
generators first, then the extra functions (strips) needed to close the system;  dF = sum_k dlog W_k M[k] . F  with
rational M[k].  Writes basis.wl (the functions F[i] in II[...]), deq.wl (M[k] as SparseArrays) and dependent.wl (the
generators that are not in F, as exact combinations of F).   usage: exact_deq.py system.pkl Generators.wl outdir"""
import pickle, sys
from fractions import Fraction
import flint
import iisys

S = pickle.load(open(sys.argv[1], 'rb'))
gens = dict(iisys.parse_generators(sys.argv[2]))
out = sys.argv[3]
names, wt = S['names'], list(S['wt'])
assert list(gens) == names, 'Generators.wl does not match system.pkl'
basis, origin = S['basis'], S['origin']
P = iisys.P

# labels: generators in the basis keep their names; extras F[w, j]
order, label = [], []
for w in range(6, -1, -1):
    nx = 0
    for f, o in zip(basis.get(w, []), origin.get(w, [])):
        order.append((w, f))
        if o is not None:
            label.append(names[o])
        elif w == 0:
            label.append('One')
        else:
            nx += 1
            label.append('F[%d,%d]' % (w, nx))
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


class Solver:
    """exact c with g = sum_j c_j vecs[j], checked on every word"""

    def __init__(self, vecs):
        self.vecs = vecs
        self.pw = pivots(vecs)
        d = len(vecs)
        B = flint.fmpq_mat(d, d, [flint.fmpq(vecs[j].get(self.pw[r], 0).numerator, vecs[j].get(self.pw[r], 0).denominator)
                                  for r in range(d) for j in range(d)])
        self.Binv = B.inv()

    def __call__(self, g, what):
        d = len(self.vecs)
        rhs = flint.fmpq_mat(d, 1, [flint.fmpq(g.get(p, Fraction(0)).numerator, g.get(p, Fraction(0)).denominator)
                                    for p in self.pw])
        c = self.Binv * rhs
        cf = [Fraction(int(c[j, 0].p), int(c[j, 0].q)) for j in range(d)]
        chk = dict(g)
        for j, x in enumerate(cf):
            if x:
                for word, v in self.vecs[j].items():
                    chk[word] = chk.get(word, 0) - x * v
        assert not any(chk.values()), what
        return cf


M = {}                                          # (k, i, j) -> Fraction
for w in range(6, 0, -1):
    solve = Solver([order[i][1] for i in index[w - 1]])
    for i in index[w]:
        f = order[i][1]
        for k in range(1, iisys.NL + 1):
            g = iisys.strip(f, k)
            if not g:
                continue
            for j, x in enumerate(solve(g, 'not closed at weight %d' % w)):
                if x:
                    M[k, i, index[w - 1][j]] = x
    print('weight %d -> %d: done' % (w, w - 1), flush=True)

n = len(order)
hgt = max(max(abs(x.numerator), x.denominator) for x in M.values())
print('%d non-zero entries over the 20 letters; largest numerator/denominator %d' % (len(M), hgt))

# generators that are not basis functions: exact combinations of the basis functions of their weight
inbasis = {o for w in origin for o in origin[w] if o is not None}
dep = []
for w in sorted(set(wt)):
    todo = [i for i, x in enumerate(wt) if x == w and i not in inbasis]
    if not todo:
        continue
    solve = Solver([order[i][1] for i in index[w]])
    for i in todo:
        cf = solve(gens[names[i]], '%s is not in the span of the basis' % names[i])
        dep.append((names[i], [(index[w][j], x) for j, x in enumerate(cf) if x]))
print('%d generators are not in the basis (dependent)' % len(dep))


def wl(x):
    return str(x.numerator) if x.denominator == 1 else '%d/%d' % (x.numerator, x.denominator)


with open('%s/basis.wl' % out, 'w') as fh:
    fh.write('(* F = basis of the closed system, weight 6 -> 0; the independent generators first at each weight, F[w,j] =\n'
             '   extra functions (combinations of strips of the generators), One = II[] = 1.  List of label -> function. *)\n{\n')
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
with open('%s/dependent.wl' % out, 'w') as fh:
    fh.write('(* generators that are linearly dependent on the others, as exact combinations of the basis F of basis.wl.\n'
             '   List of name -> combination of labels. *)\n{\n')
    fh.write(',\n'.join('%s -> %s' % (nm, ' + '.join('(%s)*%s' % (wl(x), label[j]) for j, x in cf)) for nm, cf in dep))
    fh.write('\n}\n')
print('wrote %s/basis.wl, %s/deq.wl and %s/dependent.wl' % (out, out, out))
