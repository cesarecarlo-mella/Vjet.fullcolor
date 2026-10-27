"""generators: build the closed system, start values at c1 (both nested regularisations), carry them to P0 with
tridents_bs, compare, write the production data file.   usage: gen_run.py Generators.wl outdir"""
import pickle, subprocess, sys, time
import numpy as np
import iisys

src, out = sys.argv[1:3]
BS = 'cpp/tridents_bs'
P0 = (0.3, 0.3)
t0 = time.time()
funcs = iisys.parse_generators(src)
S = iisys.build(funcs)
NT = len(funcs)
print('built in %.1fs' % (time.time() - t0))
pickle.dump(S, open('%s/system.pkl' % out, 'wb'))


def run(fn, pts, tol=1e-13, delta=(0.1, 0.2)):
    inp = '\n'.join('%.17g %.17g' % p for p in pts)
    r = subprocess.run([BS, fn, str(tol), str(tol), str(delta[0]), str(delta[1])], input=inp, capture_output=True, text=True)
    lines = r.stdout.split('\n')
    res, i = [], 0
    while i < len(lines) and lines[i].startswith('#'):
        h = lines[i].split(); n = 0; i += 1
        V = []
        while i < len(lines) and lines[i] and not lines[i].startswith('#'):
            a, b = lines[i].split(); V.append(complex(float(a), float(b))); i += 1
        res.append((h, np.array(V)))
    return res


at = {}
for base in ('x', 'y'):
    t0 = time.time()
    pt, yv, T = iisys.start_values(S, base)
    fn = '%s/start_%s.bin' % (out, base)
    ntot = iisys.export(S, fn, pt, yv, T, readout_y=True)
    for delta in ((0.1, 0.2), (-0.1, -0.2)):
        (h, V), = run(fn, [P0], delta=delta)
        V = V.reshape(7, ntot)
        at[base, delta] = V
        print('base %s: start %s (series %.1fs), to P0 with delta %s: %s steps %s evaluations'
              % (base, pt, time.time() - t0, delta, h[3], h[4]))


def gens(V):
    return np.array([V[w, i] for i, w in enumerate(S['wt'])])


ref = gens(at['x', (0.1, 0.2)])
for key, V in at.items():
    d = np.abs(gens(V) - ref) / np.maximum(np.abs(ref), 1)
    print('generators at P0, %s vs (x, +delta): max rel diff %.1e  (worst %s)' % (key, d.max(), S['names'][int(d.argmax())]))
dims = S['dims']
for key, V in at.items():
    if key == ('x', (0.1, 0.2)):
        continue
    for w in range(1, 6):
        o = NT + sum(dims[1:w])
        a, b = V[w, o:o + dims[w]], at['x', (0.1, 0.2)][w, o:o + dims[w]]
        print('  basis y_%d, %s vs (x,+delta): max diff %.1e' % (w, key, np.abs(a - b).max()))

# production file: start at P0 from the x-regularised start, + delta
V = at['x', (0.1, 0.2)]
yv = {0: np.array([1.0 + 0j])}
for w in range(1, 6):
    o = NT + sum(dims[1:w])
    yv[w] = V[w, o:o + dims[w]]
T = V[:, :NT].copy()
iisys.export(S, '%s/generators_bs.bin' % out, P0, yv, T)
np.save('%s/P0_values.npy' % out, gens(V))
print('wrote %s/generators_bs.bin' % out)
