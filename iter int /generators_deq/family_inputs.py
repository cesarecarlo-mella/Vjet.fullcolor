"""the iterated-integral system as a 23rd family for the blow-up pipeline, from out/deq.wl:
    atilde.m      Atilde_ij = sum_k M[k]_ij l[k]   (dF = dAtilde . F, l[k] = Log[W_k] of alphabet.wl; the weight plays the
                  role of the eps order: M[k] lowers the weight by one)
    labels.m      the functions in the order of atilde (basis.wl: their definition in II[...]; generators.wl: the
                  generators in terms of them)
    sol0.m        weight 0: the constant B[0,1] (= 1/96 in the integral basis), all others 0
    bdy_c1_w<w>.m the constants at the corner c1 (x, y -> 0) in chart 1 (x = t y^2): 0 for w = 1..6 (base point,
                  shuffle regularised: only the log-free constant term vanishes, the logs come from the DE)
    edge_x0.m     the connection restricted to the edge x = 0 (s = y) and edge_y0.m to the edge y = 0 (s = x):
                  {A0, A1} with dF/ds = (A0/s + A1/(s - 1)) . F, for the transport from c1 to c3 (x = 0, s -> 1)
                  and to c2 (y = 0, s -> 1).  Every letter in use is c s^a (1 - s)^b on the edge; a letter that vanishes
                  identically there (x on x = 0, y on y = 0) and the odd letters 7, 8 (-> constants) drop out.
usage: family_inputs.py <deq.wl> <outdir>"""
import os, re, sys
from fractions import Fraction
import sympy as sp

deq, out = sys.argv[1:3]
os.makedirs(out, exist_ok=True)
src = open(deq).read()
src = re.sub(r'\(\*.*?\*\)', '', src, flags=re.S)
labels = re.findall(r'T\[\d+\]|[BFS]\[\d+,\d+\]|One', src[src.index('{{') + 2: src.index('SparseArray')])
n = len(labels)
M = {}
for k, body in enumerate(re.findall(r'SparseArray\[\{(.*?)\}, \{\d+, \d+\}\]', src, re.S), 1):
    for i, j, v in re.findall(r'\{(\d+),(\d+)\}->(-?\d+(?:/\d+)?)', body):
        M[k, int(i) - 1, int(j) - 1] = Fraction(v)
print('%d functions, %d non-zero entries, letters %s' % (n, len(M), sorted({k for k, _, _ in M})))


def wl(x):
    return str(x.numerator) if x.denominator == 1 else '%d/%d' % (x.numerator, x.denominator)


# atilde: dense rows (as the other families), entries sum_k c l[k]
rows = [[[] for _ in range(n)] for _ in range(n)]
for (k, i, j), c in sorted(M.items()):
    rows[i][j].append('(%s)*l[%d]' % (wl(c), k))
with open(os.path.join(out, 'atilde.m'), 'w') as f:
    f.write('{' + ',\n '.join('{' + ', '.join(' + '.join(e) if e else '0' for e in r) + '}' for r in rows) + '}\n')
with open(os.path.join(out, 'labels.m'), 'w') as f:
    f.write('(* the functions of atilde.m, in order; their definitions in II[...]: generators_deq/out/basis.wl *)\n{'
            + ', '.join('"%s"' % s for s in labels) + '}\n')
# weight 0: One = 1, or the constant B[0,1] of the integral basis (its value from basis.wl)
w0 = re.search(r'B\[0,1\] -> \((-?\d+(?:/\d+)?)\)\*II\[\]', open(os.path.join(os.path.dirname(deq), 'basis.wl')).read())
c0 = w0.group(1) if w0 else '1'
with open(os.path.join(out, 'sol0.m'), 'w') as f:
    f.write('{' + ', '.join(c0 if s in ('One', 'B[0,1]') else '0' for s in labels) + '}\n')
for w in range(1, 7):
    with open(os.path.join(out, 'bdy_c1_w%d.m' % w), 'w') as f:
        f.write('{' + ', '.join(['0'] * n) + '}\n')

# edge connections
x, y, s = sp.symbols('x y s')
LET = {1: x, 2: y, 3: 1 - x - y, 4: 1 - x, 5: 1 - y, 6: x + y, 9: 1 - (1 + x) * y, 10: 1 - x * (1 + y), 11: x - y + y**2,
       12: -x + x**2 + y, 13: (x - 1)**2 - y, 14: -y + (x + y)**2, 15: -x + (x + y)**2, 16: -x + (1 - y)**2,
       17: x + x * y + y**2, 18: x**2 + y + x * y, 19: (y - 1)**2 + x * y, 20: (x - 1)**2 + x * y}
used = sorted({k for k, _, _ in M})
for name, sub, desc in (('edge_x0', {x: 0, y: s}, 'x = 0, s = y (from c1 to c3)'), ('edge_y0', {x: s, y: 0}, 'y = 0, s = x (from c1 to c2)')):
    ab, table = {}, []
    for k in used:
        if k in (7, 8):
            ab[k] = (0, 0); table.append('W%d -> constant (odd letter, Sqrt[x y z] -> 0)' % k)
            continue
        e = sp.factor(LET[k].subs(sub))
        if e == 0:
            ab[k] = (0, 0); table.append('W%d -> 0 identically (dropped: the regularisation)' % k)
            continue
        c, fac = sp.factor_list(e)
        a = b = 0
        for p, m in fac:
            if sp.expand(p - s) == 0:
                a += m
            elif sp.expand(p - (1 - s)) == 0 or sp.expand(p - (s - 1)) == 0:
                b += m
            else:
                sys.exit('!! letter %d on %s: factor %s not in {s, 1 - s}' % (k, name, p))
        ab[k] = (a, b)
        table.append('W%d -> %s' % (k, sp.sstr(e)))
    A = {0: {}, 1: {}}
    for (k, i, j), c in M.items():
        for side in (0, 1):
            if ab[k][side]:
                A[side][(i, j)] = A[side].get((i, j), 0) + ab[k][side] * c
    with open(os.path.join(out, name + '.m'), 'w') as f:
        f.write('(* the connection on the edge %s:  dF/ds = (A0/s + A1/(s - 1)) . F,  F in the order of labels.m.\n'
                '   A0 = Sum_k a_k M[k], A1 = Sum_k b_k M[k], with W_k = c s^a_k (1 - s)^b_k on the edge:\n   %s\n'
                '   {A0, A1} *)\n' % (desc, '\n   '.join(table)))
        f.write('{' + ',\n'.join('SparseArray[{%s}, {%d, %d}]' % (', '.join('{%d,%d}->%s' % (i + 1, j + 1, wl(v))
                for (i, j), v in sorted(A[side].items()) if v), n, n) for side in (0, 1)) + '}\n')
    print('%s: %s; A0 %d, A1 %d non-zero entries' % (name, ', '.join('W%d:%s' % (k, ab[k]) for k in used),
                                                      sum(1 for v in A[0].values() if v), sum(1 for v in A[1].values() if v)))
print('wrote', out, '; sol0: weight-0 function =', c0)
