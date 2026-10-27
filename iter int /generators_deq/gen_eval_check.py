"""generators_bs.bin over the triangle: timing, and accuracy from (i) tolerance 1e-10 vs 1e-13, (ii) the other side of
the path deformation, (iii) the direct c1 series at points near the edges x = 0 and y = 0 (independent of P0)."""
import pickle, subprocess
import numpy as np
import iisys

S = pickle.load(open('out/system.pkl', 'rb'))
wt = S['wt']; NT = len(wt)
BS = 'cpp/tridents_bs'


def run(pts, tol=1e-10, delta=(0.1, 0.2)):
    inp = '\n'.join('%.17g %.17g' % tuple(p) for p in pts)
    r = subprocess.run([BS, 'out/generators_bs.bin', str(tol), str(tol), str(delta[0]), str(delta[1])],
                       input=inp, capture_output=True, text=True)
    out, cur = [], None
    for l in r.stdout.split('\n'):
        if l.startswith('#'):
            cur = []; out.append((l.split(), cur))
        elif l:
            a, b = l.split(); cur.append(complex(float(a), float(b)))
    return [(float(h[5]), int(h[3]), np.array([np.array(v).reshape(7, NT)[w, i] for i, w in enumerate(wt)])) for h, v in out]


def rel(a, b):
    return np.abs(a - b) / np.maximum(np.abs(b), 1)


rng = np.random.default_rng(7)
pts = []
while len(pts) < 60:
    x, y = rng.random(2)
    if x + y < 1:
        pts.append((x, y))
edge = [(1e-2, 0.5), (0.5, 1e-2), (0.5, 0.49), (1e-3, 0.3), (0.3, 1e-3), (0.6, 0.399), (0.02, 0.02), (0.96, 0.02), (0.02, 0.96)]
for label, P in (('60 random points', pts), ('near edges/corners', edge)):
    a = run(P); b = run(P, 1e-13); c = run(P, 1e-10, (-0.1, -0.2))
    t = [r[0] for r in a]
    d1 = np.array([rel(x[2], y[2]) for x, y in zip(a, b)])
    d2 = np.array([rel(x[2], y[2]) for x, y in zip(c, b)])
    print('%s: time median %.3f s max %.3f s' % (label, np.median(t), np.max(t)))
    for w in sorted(set(wt)):
        sel = wt == w
        print('  weight %d: tol 1e-10 vs 1e-13: median %.0e worst %.0e | other deformation side: worst %.0e'
              % (w, np.median(d1[:, sel]), d1[:, sel].max(), d2[:, sel].max()))
    mag = np.array([np.abs(r[2]) for r in b])
    print('  |generator| over these points: median %.2g, min %.2g, max %.2g' % (np.median(mag), mag.min(), mag.max()))

# (iii) direct series at the near-edge points of the start construction
for base, eps, edge_v in (('x', 0.01, 0.4), ('y', 0.01, 0.4), ('x', 0.005, 0.7)):
    pt, yv, T = iisys.start_values(S, base, edge=edge_v, eps=eps)
    ser = np.array([T[w, i] for i, w in enumerate(wt)])
    (_, _, v), = run([pt], 1e-13)
    print('direct c1 series (%s-first) at %s vs BS from P0: worst rel diff %.1e' % (base, pt, rel(v, ser).max()))
