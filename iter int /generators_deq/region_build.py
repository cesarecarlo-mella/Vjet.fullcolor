"""the numerical evaluator in a production region A > 1 (A = x, y or z; the other two variables negative), as a decay-like
triangle (notes/continuation_notes.tex):

    A = 1/Z and the other two variables -X/Z, -Y/Z  (X + Y + Z = 1, all positive; Z = v the inverse energy):
        z > 1:  (x, y) = (-X/Z, -Y/Z)        corner X, Y -> 0 = c1 continued (chart 1)
        x > 1:  (x, y) = (1/Z, -X/Z)         corner X, Y -> 0 = c2 continued (chart 5), z = -Y/Z
        y > 1:  (x, y) = (-X/Z, 1/Z)         corner X, Y -> 0 = c3 continued (chart 3), z = -Y/Z

 1. the letter map: with the continued root Sqrt[x y z] -> -Sqrt[X Y Z]/Z^2, dlog W_k(x, y) = sum_j c_kj dlog W_j(X, Y)
    with integer c (fitted at 40 digits, checked);  the connection and read-off of the region: A'_j = sum_k c_kj A_k;
 2. the values at the corner P*: the t^0 B^0 part of the ITER series at the attached decay corner with
    Log[t] -> Log[t] - I Pi, Log[B] -> Log[B] + I Pi (the half-turn of the Feynman prescription), for the 586 functions
    of the exact system -> the generators (generators.wl) and the orthonormal basis of the numerical system;
 3. the series start from P* (iisys.start_values, X -> 0 first, as chart 1/5/3), carried to the hub P0' = (0.3, 0.3)
    of the region's triangle and stored: out/region_<A>_bs.bin;
 4. check: the generators at points of the region against a direct integration from the decay hub P0 with the decay
    data, crossing the zeros of the two small variables from above (an independent continuation).

usage: region_build.py <x|y|z> <ITER dist dir> [<ITER inputs dir>]"""
import json, os, pickle, re, subprocess, sys, time
from fractions import Fraction
import numpy as np
import mpmath as mp
import iisys

HERE = os.path.dirname(os.path.abspath(__file__))
region = sys.argv[1]
dist = sys.argv[2]
inputs = sys.argv[3] if len(sys.argv) > 3 else os.path.join(HERE, '..', '..', 'blowup_series', 'families', 'ITER', 'inputs')
CHART = {'z': 1, 'x': 5, 'y': 3}[region]
BS = os.path.join(HERE, 'cpp', 'tridents_bs')
t0 = time.time()


def phys_xy(X, Y):
    Z = 1 - X - Y
    return {'z': (-X / Z, -Y / Z), 'x': (1 / Z, -X / Z), 'y': (-X / Z, 1 / Z)}[region]


# ---------------------------------------------------------------- 1. the letter map
def letters(X, Y, R):
    Z = 1 - X - Y
    return {1: X, 2: Y, 3: Z, 4: 1 - X, 5: 1 - Y, 6: X + Y, 7: (X * Z - 1j * R) / (X * Z + 1j * R),
            8: (X * Y - 1j * R) / (X * Y + 1j * R), 9: 1 - (1 + X) * Y, 10: 1 - X * (1 + Y), 11: X - Y + Y**2,
            12: -X + X**2 + Y, 13: (X - 1)**2 - Y, 14: -Y + (X + Y)**2, 15: -X + (X + Y)**2, 16: -X + (1 - Y)**2,
            17: X + X * Y + Y**2, 18: X**2 + Y + X * Y, 19: (Y - 1)**2 + X * Y, 20: (X - 1)**2 + X * Y}


def old_at(X, Y, k):
    x, y = phys_xy(X, Y)
    Z = 1 - X - Y
    return letters(x, y, -mp.sqrt(X * Y * Z) / Z**2)[k]


def new_at(X, Y, k):
    return letters(X, Y, mp.sqrt(X * Y * (1 - X - Y)))[k]


def grad(f, X, Y, h=mp.mpf('1e-18')):
    return [(mp.log(f(X + h, Y)) - mp.log(f(X - h, Y))) / (2 * h), (mp.log(f(X, Y + h)) - mp.log(f(X, Y - h))) / (2 * h)]


mp.mp.dps = 45
pts = [(mp.mpf(a) / 97, mp.mpf(b) / 89) for a, b in ((13, 21), (27, 11), (31, 30), (19, 47), (41, 17), (9, 9), (23, 38), (35, 25),
                                                      (14, 52), (50, 20), (7, 33), (44, 9))]
Gn = [{j: grad(lambda a, b, j=j: new_at(a, b, j), X, Y) for j in range(1, 21)} for X, Y in pts]
C = np.zeros((21, 21), dtype=object)                     # C[k][j]: dlog W_k(old) = sum_j C[k][j] dlog W_j(new)
for k in range(1, 21):
    rows, rhs = [], []
    for (X, Y), gn in zip(pts, Gn):
        go = grad(lambda a, b: old_at(a, b, k), X, Y)
        for comp in (0, 1):
            rows.append([gn[j][comp] for j in range(1, 21)]); rhs.append(go[comp])
    A = mp.matrix(rows); b = mp.matrix(rhs)
    sol = mp.lu_solve(A.T * A, A.T * b)
    res = mp.norm(A * sol - b) / mp.norm(b)
    ints = [int(mp.nint(mp.re(s))) for s in sol]
    err = max(abs(s - n) for s, n in zip(sol, ints))
    if res > 1e-25 or err > 1e-20:
        sys.exit('!! letter %d does not map to an integer combination (residual %.1e, distance to integers %.1e)' % (k, res, err))
    for j, n in enumerate(ints, 1):
        C[k][j] = n
cmap = {k: {j: int(C[k][j]) for j in range(1, 21) if C[k][j]} for k in range(1, 21)}
print('region %s > 1: letter map (dlog W_k -> sum_j c dlog W_j\'), the letters in use:' % region)
for k in (1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16):
    print('   W%-2d -> %s' % (k, ' '.join('%+d W%d\'' % (c, j) for j, c in sorted(cmap[k].items()))))

# ---------------------------------------------------------------- the region's system
S = pickle.load(open(os.path.join(HERE, 'out', 'system.pkl'), 'rb'))
NL = iisys.NL
S2 = dict(S)
S2['A'] = {(j, w): sum(cmap[k].get(j, 0) * S['A'][k, w] for k in range(1, NL + 1)) for j in range(1, NL + 1) for w in range(1, 6)}
S2['R'] = {(j, w): sum(cmap[k].get(j, 0) * S['R'][k, w] for k in range(1, NL + 1)) for j in range(1, NL + 1) for w in range(1, 7)}
for key in S2['A']:
    S2['A'][key] = np.asarray(S2['A'][key], float)
for key in S2['R']:
    S2['R'][key] = np.asarray(S2['R'][key], float)

# ---------------------------------------------------------------- 2. the values at the corner
labels = re.findall(r'"([^"]+)"', open(os.path.join(inputs, 'labels.m')).read())
LIDX = {l: i for i, l in enumerate(labels)}
sol0 = [Fraction(s.strip()) for s in open(os.path.join(inputs, 'sol0.m')).read().strip()[1:-1].split(',')]
CST = {'Pi': np.pi, 'Zeta[3]': 1.2020569031595942854, 'Zeta[5]': 1.0369277551433699263, 'I': 1j}
FP = np.array([complex(float(v)) for v in sol0])
sys.path.insert(0, os.path.join(HERE, '..', '..', 'blowup_series', 'code'))
import pipe_common as pc           # noqa: E402
for w in range(1, 7):
    for (tag, i, it, iy), q in pc.load_pickle(os.path.join(dist, 'rec_c%d_w%d.pkl' % (CHART, w))).items():
        pt, py, la, lb, ck = tag
        if it == 1 and iy == 1 and pt == 0 and py == 0:
            c = complex(float(q))
            for sym, e2 in ck:
                c *= CST[sym] ** (e2 / 2)
            FP[i] += c * (-1j * np.pi) ** la * (1j * np.pi) ** lb
GEN = {}
for line in open(os.path.join(HERE, 'out', 'generators.wl')):
    m = re.match(r'(T\[\d+\]) -> (.*?),?\s*$', line.strip())
    if m:
        GEN[m.group(1)] = [(float(Fraction(c)), LIDX[l]) for c, l in re.findall(r'\((-?[\d/]+)\)\*([A-Z]\[[\d,]+\])', m.group(2))]
names, wt = S['names'], S['wt']
TP = np.array([sum(c * FP[i] for c, i in GEN[n]) for n in names])
NT = len(names)
T_init = np.zeros((7, NT), complex)
for g, w in enumerate(wt):
    T_init[w, g] = TP[g]
# the orthonormal basis of the numerical system at the corner: from the exact functions (basis.wl)
fvec = {}
for line in open(os.path.join(HERE, 'out', 'basis.wl')):
    m = re.match(r'([A-Z]\[[\d,]+\]) -> (.*?),?\s*$', line.strip())
    if m:
        fvec[m.group(1)] = {tuple(int(a) for a in wd.split(',')) if wd else (): float(Fraction(c))
                            for c, wd in re.findall(r'\((-?[\d/]+)\)\*II\[([\d,]*)\]', m.group(2))}
y_init = {0: np.array([1.0 + 0j])}
for w in range(1, 6):
    words = sorted({wd for f in S['basis'][w] for wd in f})
    sup = {wd: i for i, wd in enumerate(words)}
    Mb = np.zeros((len(words), len(S['basis'][w])))
    for j, f in enumerate(S['basis'][w]):
        for wd, c in f.items():
            Mb[sup[wd], j] = float(c)
    Q, _ = np.linalg.qr(Mb)
    rows, rhs = [], []
    for lab, f in fvec.items():
        if len(next(iter(f))) != w:
            continue
        v = np.zeros(len(words))
        for wd, c in f.items():
            v[sup[wd]] = c
        rows.append(Q.T @ v); rhs.append(FP[LIDX[lab]])
    yv, *_ = np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)
    resid = np.abs(np.array(rows) @ yv - np.array(rhs)).max()
    y_init[w] = yv
    print('   corner values, weight %d: %d functions -> %d basis values, residual %.1e, largest |value| %.3g'
          % (w, len(rows), len(yv), resid, np.abs(rhs).max()))

# ---------------------------------------------------------------- 3. start series, transport to the hub, export
P0 = (0.3, 0.3)


def run(fn, pts, tol=1e-13, delta=(0.1, 0.2)):
    inp = '\n'.join('%.17g %.17g' % p for p in pts)
    r = subprocess.run([BS, fn, str(tol), str(tol), '%.17g' % delta[0], '%.17g' % delta[1]], input=inp, capture_output=True, text=True)
    if r.returncode:
        sys.exit('!! tridents_bs: ' + r.stderr[-300:])
    res, cur = [], None
    for l in r.stdout.split('\n'):
        if l.startswith('#'):
            cur = []; res.append(cur)
        elif l:
            a_, b_ = l.split(); cur.append(complex(float(a_), float(b_)))
    return [np.array(v) for v in res]


out = os.path.join(HERE, 'out')
pt, yv, T = iisys.start_values(S2, 'x', y_init=y_init, T_init=T_init)
fn = os.path.join(out, 'region_%s_start.tmp' % region)
ntot = iisys.export(S2, fn, pt, yv, T, readout_y=True)
dims = S['dims']
V = {}
for delta in ((0.1, 0.2), (-0.1, -0.2)):
    V[delta] = run(fn, [P0], delta=delta)[0].reshape(7, ntot)
os.remove(fn)
gens = lambda Vv: np.array([Vv[w, i] for i, w in enumerate(wt)])
d = np.abs(gens(V[(0.1, 0.2)]) - gens(V[(-0.1, -0.2)])) / np.maximum(np.abs(gens(V[(0.1, 0.2)])), 1)
print('   start at %s (series), to the hub %s: the two deformation sides agree to %.1e' % (pt, P0, d.max()))
Vv = V[(0.1, 0.2)]
yh = {0: np.array([1.0 + 0j])}
for w in range(1, 6):
    o = NT + sum(dims[1:w])
    yh[w] = Vv[w, o:o + dims[w]]
datafile = os.path.join(out, 'region_%s_bs.bin' % region)
iisys.export(S2, datafile, P0, yh, Vv[:, :NT].copy())
json.dump(dict(region=region, chart=CHART, letter_map={str(k): v for k, v in cmap.items()}, hub=P0,
               map='(x, y) = ' + {'z': '(-X/Z, -Y/Z)', 'x': '(1/Z, -X/Z)', 'y': '(-X/Z, 1/Z)'}[region]),
          open(os.path.join(out, 'region_%s.json' % region), 'w'), indent=1)
print('   wrote %s  (%.0f s)' % (datafile, time.time() - t0), flush=True)

# ---------------------------------------------------------------- 4. check against a direct i0 continuation
decay = os.path.join(out, 'generators_bs.bin')
NTd = NT
test = [(0.3, 0.3), (0.2, 0.5), (0.5, 0.2), (0.1, 0.15), (0.15, 0.7), (0.6, 0.3)]
reg = run(datafile, test, tol=1e-12)
worst = 0
for (X, Y), vr in zip(test, reg):
    x, y = phys_xy(X, Y)
    Tr = np.array([vr.reshape(7, NT)[w, i] for i, w in enumerate(wt)])
    # deformation: the small variables (those that change sign) cross zero with a positive imaginary part
    best = None
    for dx in (-0.3, -0.2, -0.1, 0.1, 0.2, 0.3):
        for dy in (-0.3, -0.2, -0.1, 0.1, 0.2, 0.3):
            ts = np.linspace(0, 1, 20001)
            gx = P0[0] + (ts + 4j * dx * ts * (1 - ts)) * (x - P0[0])
            gy = P0[1] + (ts + 4j * dy * ts * (1 - ts)) * (y - P0[1])
            gz = 1 - gx - gy
            ok = True
            for g in (gx, gy, gz):
                s = np.sign(g.real)
                for i in np.flatnonzero(s[1:] != s[:-1]):
                    if s[i] > 0 and g[i].imag <= 0:          # a variable leaving the decay region: from above
                        ok = False
            if ok and (best is None or abs(dx) + abs(dy) < abs(best[0]) + abs(best[1])):
                best = (dx, dy)
    vd = run(decay, [(x, y)], tol=1e-12, delta=best)[0].reshape(7, NTd)
    Td = np.array([vd[w, i] for i, w in enumerate(wt)])
    dd = np.abs(Tr - Td) / np.maximum(np.abs(Td), 1)
    worst = max(worst, dd.max())
    print('   (X, Y) = %s -> (x, y) = (%.3g, %.3g): region evaluator vs direct continuation from P0 (deformation %s): '
          'max rel diff %.1e (%s), median %.1e' % ((X, Y), x, y, best, dd.max(), names[int(dd.argmax())], np.median(dd)), flush=True)
print('region %s done: worst %.1e  (%.0f s)' % (region, worst, time.time() - t0))
