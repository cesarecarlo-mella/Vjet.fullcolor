"""quad start values in a production region A > 1 (region_build.py: (X, Y) its triangle, hub (0.3, 0.3)):

 1. the 586 functions of the exact system at the physical image of the hub, by a quad integration (cpp/iter_quad)
    from the decay quad values at P0 along a complex-deformed straight path on which the two variables that become
    negative cross zero with a positive imaginary part (the Feynman prescription; as region_build.py's check).  Two
    such paths, compared; the root must arrive at -Sqrt[X Y Z]/Z^2 (the continued root).
 2. out/region_<A>_iter_quad.bin: the region's integer system for iter_quad, M'_j = sum_k c_kj M_k
    (out/region_<A>.json), started at the hub with the region's own root r = +Sqrt[X Y Z].
 3. the head of out/region_<A>_bs.bin (the 503 system, double) from the same values (quad_p0_to_double).
usage: region_quad.py <x|y|z> [--write]"""
import itertools, json, os, re, struct, sys
from fractions import Fraction
import numpy as np
import mpmath as mp
import gen_eval
import quad_data as Q

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'out')
mp.mp.dps = 40
region = sys.argv[1]
HUB = (0.3, 0.3)
TOL_I0 = 1e-28
P0 = gen_eval.P0


def phys(X, Y):
    X, Y = mp.mpf(X), mp.mpf(Y)
    Z = 1 - X - Y
    return {'z': (-X / Z, -Y / Z), 'x': (1 / Z, -X / Z), 'y': (-X / Z, 1 / Z)}[region]


# 1. the hub values by the i0 path
xt, yt = phys(*HUB)
zt = 1 - xt - yt
z0 = 1 - P0[0] - P0[1]
neg = [v for v, a, b in (('x', P0[0], xt), ('y', P0[1], yt), ('z', z0, zt)) if a > 0 > b]
print('region %s: hub %s -> (x, y, z) = (%.6g, %.6g, %.6g); negative: %s' % (region, HUB, xt, yt, zt, neg))


def side_ok(d):
    """Im of each variable that becomes negative is > 0 on the path (Im x = 4 dx t(1-t) (x_t - x0), z = 1 - x - y)"""
    ix, iy = d[0] * float(xt - P0[0]), d[1] * float(yt - P0[1])
    im = dict(x=ix, y=iy, z=-(ix + iy))
    return all(im[v] > 0.02 for v in neg)


grid = [round(float(v), 2) for v in np.arange(-0.3, 0.31, 0.05) if abs(v) > 1e-9]
cand = sorted(((Q.path_margin(P0, (float(xt), float(yt)), (dx, dy))[0], (round(dx, 2), round(dy, 2)))
               for dx, dy in itertools.product(grid, grid) if side_ok((dx, dy))), reverse=True)
paths, res = [], []
rexp = -mp.sqrt(HUB[0] * HUB[1] * (1 - HUB[0] - HUB[1])) / (1 - HUB[0] - HUB[1]) ** 2
for m, d in cand:                       # by margin; a path whose step control fails at TOL_I0 is skipped
    if len(paths) == 2:
        break
    if paths and abs(d[0] - paths[0][0]) + abs(d[1] - paths[0][1]) < 0.15:
        continue
    try:
        (inf, F), = Q.run([(xt, yt)], os.path.join(OUT, 'iter_quad_P0.bin'), 'quad', TOL_I0, d, kmax=12)
    except RuntimeError as e:
        print('   path %s (margin %.3f): %s, skipped' % (d, m, str(e).strip()), flush=True)
        continue
    paths.append(d); res.append(F)
    print('   path %s (margin %.3f): %d steps, %.1f s; root at the end %s (expected -Sqrt[XYZ]/Z^2 = %.15g)'
          % (d, m, inf['steps'], inf['time'], inf['r'], rexp), flush=True)
    assert abs(inf['r'] - float(rexp)) < 1e-10, 'root continued to the wrong sign'
sc = max(abs(v) for v in res[0])
print('   the two paths agree to %.1e relative to max|F| = %.3g' % (float(max(abs(a - b) for a, b in zip(*res)) / sc), float(sc)))
F = [(a + b) / 2 for a, b in zip(*res)]

# 2. the region's integer system for iter_quad
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
print('   region system: letters %s, max |entry| %d' % ([j + 1 for j in range(20) if Mr[j]], max(abs(v) for m in Mr for v in m.values())))
r0 = mp.sqrt(mp.mpf(HUB[0]) * HUB[1] * (1 - mp.mpf(HUB[0]) - HUB[1]))
fq = os.path.join(OUT, 'region_%s_iter_quad.bin' % region)
with open(fq, 'wb') as f:
    f.write(b'ITERQ001' + struct.pack('<2i', Q.NF, 20))
    for m in Mr:
        ent = sorted(m.items())
        f.write(struct.pack('<i', len(ent)))
        if ent:
            f.write(struct.pack('<%di' % len(ent), *[i for (i, j), v in ent]))
            f.write(struct.pack('<%di' % len(ent), *[j for (i, j), v in ent]))
            f.write(struct.pack('<%dd' % len(ent), *[float(v) for (i, j), v in ent]))
    f.write(struct.pack('<4d', *Q._dd(mp.mpf(HUB[0])), *Q._dd(mp.mpf(HUB[1]))))
    f.write(struct.pack('<4d', *Q._dd(r0), 0.0, 0.0))
    for v in F:
        v = mp.mpc(v)
        f.write(struct.pack('<4d', *Q._dd(v.real), *Q._dd(v.imag)))
print('   written %s' % fq)

# check: from the hub with the region system to another region point = directly from P0 along an i0 path
Pt = (0.15, 0.6)
(i1, Fa), = Q.run([Pt], fq, 'quad', TOL_I0, gen_eval.choose_deformation(Pt, region=region), kmax=12)
xp, yp = phys(*Pt)
cchk = sorted(((Q.path_margin(P0, (float(xp), float(yp)), d)[0], d) for d in itertools.product(grid, grid)
               if all(dict(x=d[0] * float(xp - P0[0]), y=d[1] * float(yp - P0[1]),
                      z=-(d[0] * float(xp - P0[0]) + d[1] * float(yp - P0[1])))[v] > 0.02 for v in neg)), reverse=True)
for m, dchk in cchk[:8]:
    try:
        (i2, Fb), = Q.run([(xp, yp)], os.path.join(OUT, 'iter_quad_P0.bin'), 'quad', TOL_I0, dchk, kmax=12)
        break
    except RuntimeError:
        continue
print('   check at (X, Y) = %s: region system from the hub vs direct i0 path from P0 (%s): %.1e relative to max|F|'
      % (Pt, dchk, float(max(abs(a - b) for a, b in zip(Fa, Fb)) / max(abs(v) for v in Fb))))

# 3. the head of the double evaluator's region file
import quad_p0_to_double as QD
new_y, Tm = QD.head_values(F, verbose=False)
QD.patch_head(os.path.join(OUT, 'region_%s_bs.bin' % region), new_y, Tm, hub=HUB, write='--write' in sys.argv)
