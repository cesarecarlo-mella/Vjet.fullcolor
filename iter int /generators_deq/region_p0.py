"""the start values P0 of a production region A, defined as for the decay: the exact corner series continued into the
region, two start points near the region's corner (X, Y -> 0, the decay corner continued), and a quad transport over two
paths to the hub (0.3, 0.3) of the region's triangle.  (notes/evaluator_notes.tex)

  region z: corner c1 (x, y -> 0).  The functions are a double series there (charts 1 and 4 reversed, quad_data.double_series):
            evaluated at x = -X/Z, y = -Y/Z, both crossing zero from above: log -> log|.| + i pi, sqrt -> i sqrt|.|.
  region x: corner c2 (y, z -> 0).  No double series there (the letter y + z = 1 - x vanishes, and the non-transversal
            letters are seen): the series of chart 5 (y = t z^2, B = z), continued: y = -X/Z, z = -Y/Z, so
            t = y/z^2 < 0 with log t -> log|t| - i pi, and log B -> log|B| + i pi.  The start points have small |t|.
  region y: corner c3 (x, z -> 0): chart 3 (x = t z^2, B = z), x = -X/Z, z = -Y/Z, the same continuation.
The region's integer system (M'_j = sum_k c_kj M_k, out/region_<A>.json) carries the values to the hub along two paths;
the two results are compared and averaged.  Cross-check: the hub values of the direct i0 transport from the decay P0
(region_quad.py).  --write: out/region_<A>_iter_quad.bin (start at the hub), the head of out/region_<A>_bs.bin, and
out/region_<A>_P0_values_quad.m.
usage: region_p0.py <x|y|z> [--write]"""
import json, os, re, struct, sys
from fractions import Fraction
import numpy as np
import mpmath as mp
import gen_eval
import quad_data as Q

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'out')
DIST = os.path.join(HERE, '..', '..', 'blowup_series', 'families', 'ITER')
mp.mp.dps = 40
region = sys.argv[1]
HUB = (0.3, 0.3)
TOL = 1e-28
STARTS = {'z': [(0.05, 0.05), (0.04, 0.06)], 'x': [(0.001, 0.2), (0.0005, 0.15)], 'y': [(0.001, 0.2), (0.0005, 0.15)]}[region]


def series_at(S, top=None):
    """the 586 functions at the start point S = (X, Y) of the region, from the continued corner series.
    top = (dt, dB): drop the dt highest orders in the first variable and dB in the second (truncation check)"""
    X, Y = mp.mpf(S[0]), mp.mpf(S[1])
    Z = 1 - X - Y
    a, b = X / Z, Y / Z                                  # |small variables| at the corner
    if region == 'z':
        ser, point, ph = Q.double_series(os.path.join(DIST, 'dist_16_48')), (a, b), (mp.pi, mp.pi)
    else:
        chart = {'x': 5, 'y': 3}[region]
        ser = Q.chart_series(os.path.join(DIST, 'dist_16_48_c35'), chart)
        point, ph = (a / b ** 2, b), (-mp.pi, mp.pi)    # |t| = |SMALL| / BIG^2, |B| = |BIG|
    if top is None:
        return Q.series_values(ser, point, phases=ph)
    if region == 'z':                                    # double series: truncate in the total degree instead
        return Q.series_values(ser, point, maxdeg=32 - top[0] - top[1], phases=ph)
    idx = np.concatenate([g[0] for g in ser['groups'].values()])
    am, bm = int(idx[:, 1].max()), int(idx[:, 2].max())
    sub = {'groups': {ck: _keep(g, am - 2 * top[0], bm - 2 * top[1]) for ck, g in ser['groups'].items()}}
    return Q.series_values(sub, point, phases=ph)


def _keep(g, amax2, bmax2):
    idx, cf = g
    m = (idx[:, 1] <= amax2) & (idx[:, 2] <= bmax2)
    return idx[m], [c for c, k in zip(cf, m) if k]


def region_matrices():
    src = re.sub(r'\(\*.*?\*\)', '', open(os.path.join(OUT, 'deq.wl')).read(), flags=re.S)
    mats = re.findall(r'SparseArray\[\{(.*?)\}, \{(\d+), \d+\}\]', src, re.S)
    M = [{(int(i) - 1, int(j) - 1): Fraction(v) for i, j, v in re.findall(r'\{(\d+),(\d+)\}->(-?\d+(?:/\d+)?)', body)} for body, _ in mats]
    cmap = json.load(open(os.path.join(OUT, 'region_%s.json' % region)))['letter_map']
    Mr = [dict() for _ in range(20)]
    for k in range(1, 21):
        for j, c in cmap[str(k)].items():
            for ij, v in M[k - 1].items():
                Mr[int(j) - 1][ij] = Mr[int(j) - 1].get(ij, 0) + c * v
    Mr = [{ij: v for ij, v in m.items() if v} for m in Mr]
    assert all(v.denominator == 1 for m in Mr for v in m.values())
    return Mr


def write_region_data(fn, Mr, start, F):
    """ITERQ001 file: the region system, start (X, Y) with the region's root +Sqrt[X Y Z], values F there"""
    X, Y = mp.mpf(start[0]), mp.mpf(start[1])
    r0 = mp.sqrt(X * Y * (1 - X - Y))
    with open(fn, 'wb') as f:
        f.write(b'ITERQ001' + struct.pack('<2i', Q.NF, 20))
        for m in Mr:
            ent = sorted(m.items())
            f.write(struct.pack('<i', len(ent)))
            if ent:
                f.write(struct.pack('<%di' % len(ent), *[i for (i, j), v in ent]))
                f.write(struct.pack('<%di' % len(ent), *[j for (i, j), v in ent]))
                f.write(struct.pack('<%dd' % len(ent), *[float(v) for (i, j), v in ent]))
        f.write(struct.pack('<4d', *Q._dd(X), *Q._dd(Y)))
        f.write(struct.pack('<4d', *Q._dd(r0), 0.0, 0.0))
        for v in F:
            v = mp.mpc(v)
            f.write(struct.pack('<4d', *Q._dd(v.real), *Q._dd(v.imag)))


rel = lambda A, B, sc: float(max(abs(a - b) for a, b in zip(A, B)) / sc)
Mr = region_matrices()
print('region %s: hub %s; start points %s' % (region, HUB, STARTS))
used = []
res = []
for n, S in enumerate(STARTS):
    F = series_at(S)
    sc = max(abs(v) for v in F)
    tr = [rel(F, series_at(S, top), sc) for top in ((2, 0), (0, 4))]   # z: total degree <= 30 / <= 28
    tmp = os.path.join(OUT, 'region_%s_start_%d.bin' % (region, n))
    write_region_data(tmp, Mr, S, F)
    cand = sorted(((Q.path_margin(S, HUB, d, region)[0], d) for d in gen_eval.DEFORM), reverse=True)
    for m, d in cand:
        if d in used:
            continue
        try:
            (inf, FH), = Q.run([HUB], tmp, 'quad', TOL, d, kmax=12)
        except RuntimeError as e:
            print('   start %s, path %s: %s' % (S, d, e)); continue
        used.append(d); res.append(FH)
        what = 'total degree <= 30 / <= 28' if region == 'z' else 'top 2 orders in t / 4 in B dropped'
        print('   start %s: change of the series with %s: %.0e / %.0e; path %s (margin %.3f): %d steps, %.1f s'
              % (S, what, tr[0], tr[1], d, m, inf['steps'], inf['time']), flush=True)
        break
    os.remove(tmp)
sc = max(abs(v) for v in res[0])
print('   the two determinations at the hub agree to %.1e relative to max|F| = %.3g' % (rel(res[0], res[1], sc), float(sc)))
P0 = [(a + b) / 2 for a, b in zip(*res)]
(i0, Fi0), = Q.run([HUB], os.path.join(OUT, 'region_%s_iter_quad.bin' % region), 'quad', TOL, (0.1, 0.2), kmax=12)
print('   cross-check: vs the hub values in out/region_%s_iter_quad.bin (direct i0 transport from the decay P0): %.1e'
      % (region, rel(P0, Fi0, sc)))
if '--write' in sys.argv:
    write_region_data(os.path.join(OUT, 'region_%s_iter_quad.bin' % region), Mr, HUB, P0)
    with open(os.path.join(OUT, 'region_%s_P0_values_quad.m' % region), 'w') as f:
        f.write('(* region %s: the 586 functions of out/basis.wl, continued, at the hub (X, Y) = (3/10, 3/10) of the region (as\n'
                '   doubles), region root +Sqrt[X Y Z]; from the continued corner series at %s and %s, quad transport over two\n'
                '   paths (region_p0.py); agreement %.0e relative to max|F| = %.3g *)\n{\n'
                % (region, STARTS[0], STARTS[1], rel(res[0], res[1], sc), float(sc)))
        f.write(',\n'.join('%s -> %s + %s*I' % (l, mp.nstr(v.real, 30).replace('e', '*^'), mp.nstr(v.imag, 30).replace('e', '*^'))
                           for l, v in zip(Q._labels, P0)))
        f.write('\n}\n')
    import quad_p0_to_double as QD
    new_y, Tm = QD.head_values(P0, verbose=False)
    QD.patch_head(os.path.join(OUT, 'region_%s_bs.bin' % region), new_y, Tm, hub=HUB, write=True)
    print('   written: out/region_%s_iter_quad.bin, out/region_%s_P0_values_quad.m, head of out/region_%s_bs.bin' % ((region,) * 3))
