"""validation of iisys (II convention, c1 regularisation, series start, export, BS) on a subset of the 1318 tridents
against the exact corner series (reference_series_all_weights.npz).   usage: trident_check.py [ntri]"""
import collections, re, subprocess, sys, time
from fractions import Fraction
import numpy as np
import mpmath as mp
import iisys

A = '/Users/cesarecarlomella/Desktop/Physics/PhD/XjetProject/INTEGRALS/4P1M/3L1M/trident_deq/2610.anc/trident_basis/'
REF = '/Users/cesarecarlomella/Desktop/Physics/PhD/XjetProject/INTEGRALS/4P1M/3L1M/trident_deq/trident_de/bs/reference_series_all_weights.npz'
BS = 'cpp/tridents_bs'
ntri = int(sys.argv[1]) if len(sys.argv) > 1 else 24
SQRT = {1, 2, 3, 15, 16, 17, 1283}
P0 = (0.3, 0.3)

t0 = time.time()
lab = open(A + 'iterated_integral_token_labels.wl').read().replace('\n', '').replace(' ', '')
lab = re.findall(r'((?:eps\^\d+\*)?(?:[A-Za-z]+(?:\[\d+\])?(?:\^\d+)?\*)*)II\[([\d,]*)\]', lab)
toks = []
for pre, w in lab:
    m = re.match(r'(?:eps\^(\d+)\*)?(.*)', pre)
    toks.append((int(m.group(1) or 0), m.group(2).rstrip('*'), tuple(int(a) for a in w.split(',')) if w else ()))
pos = np.array(re.findall(r'\{(\d+),(\d+)\}', open(A + 'trident_to_token_positions.wl').read()), dtype=np.int64)
val = re.findall(r'-?\d+(?:/\d+)?', open(A + 'trident_to_token_values.wl').read())
assert len(pos) == len(val), (len(pos), len(val))
print('%d tokens, %d entries (%.0fs)' % (len(toks), len(pos), time.time() - t0))
rows = collections.defaultdict(list)
for (i, j), v in zip(pos, val):
    rows[int(i)].append((int(j) - 1, v))

# subset: cover every letter at the highest weights, then fill at random
rng = np.random.default_rng(3)
letters_of = {i: {a for j, _ in r for a in toks[j][2]} for i, r in rows.items()}
chosen = []
for k in list(range(7, 21)) + list(range(1, 7)):
    c = [i for i in rows if i not in SQRT and k in letters_of[i] and i not in chosen]
    if c:
        chosen.append(int(rng.choice(c)))
rest = [i for i in rows if i not in SQRT and i not in chosen]
chosen += [int(i) for i in rng.choice(rest, max(0, ntri - len(chosen)), replace=False)]
chosen = sorted(chosen)[:max(ntri, 20)]
print('tridents:', chosen)

CONST = {'': mp.mpf(1)}
def cval(c):
    if c not in CONST:
        v = mp.mpf(1)
        for f in c.split('*'):
            m = re.match(r'(Pi|Zeta\[(\d+)\])(?:\^(\d+))?$', f)
            base = mp.pi if m.group(1) == 'Pi' else mp.zeta(int(m.group(2)))
            v *= base ** int(m.group(3) or 1)
        CONST[c] = v
    return CONST[c]

funcs, parts = [], []                                  # parts: (trident, eps order, const) per function
for i in chosen:
    g = collections.defaultdict(lambda: collections.defaultdict(Fraction))
    for j, v in rows[i]:
        e, c, w = toks[j]
        g[e, c][w] += Fraction(v)
    for (e, c), f in g.items():
        f = {w: x for w, x in f.items() if x}
        if f:
            funcs.append(('T%d_e%d_%s' % (i, e, c), f)); parts.append((i, e, c))
print('%d pure parts' % len(funcs))
S = iisys.build(funcs)
NT = len(funcs)


def run(fn, pts, tol=1e-12):
    inp = '\n'.join('%.17g %.17g' % tuple(p) for p in pts)
    r = subprocess.run([BS, fn, str(tol), str(tol)], input=inp, capture_output=True, text=True)
    out, cur = [], None
    for l in r.stdout.split('\n'):
        if l.startswith('#'):
            cur = []; out.append((l, cur))
        elif l:
            a, b = l.split(); cur.append(complex(float(a), float(b)))
    return [(h, np.array(v)) for h, v in out]


pt, yv, T = iisys.start_values(S, 'x')
nt = iisys.export(S, 'out/tri_start.bin', pt, yv, T, readout_y=True)
(h, V), = run('out/tri_start.bin', [P0], 1e-13)
V = V.reshape(7, nt)
dims = S['dims']
yP = {0: np.array([1.0 + 0j])}
for w in range(1, 6):
    o = NT + sum(dims[1:w]); yP[w] = V[w, o:o + dims[w]]
iisys.export(S, 'out/tri_bs.bin', P0, yP, V[:, :NT])

R = np.load(REF)
res = run('out/tri_bs.bin', R['pts'])
print('point          corner  steps   worst rel. error per eps order 0..6 over the %d tridents' % len(chosen))
worst = np.zeros(7)
for (x, y), K, ref, (h, Vp) in zip(R['pts'], R['corner'], R['T'], res):
    Vp = Vp.reshape(7, NT)
    tri = collections.defaultdict(complex)
    for n, (i, e, c) in enumerate(parts):
        tri[i, e] += complex(cval(c)) * Vp[S["wt"][n], n]
    err = np.zeros(7)
    for i in chosen:
        for e in range(7):
            r = ref[e, i - 1]
            err[e] = max(err[e], abs(tri[i, e] - r) / max(abs(r), 1))
    worst = np.maximum(worst, err)
    print('(%.2f, %.2f)  %s  %5s   %s' % (x, y, K, h.split()[3], ' '.join('%.0e' % v for v in err)))
print('worst: ' + ' '.join('%.1e' % v for v in worst))
