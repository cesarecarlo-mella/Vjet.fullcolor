"""all 84 remainders (every colour layer, real and imaginary parts) at a point of a production region, three ways:
  A  the evaluator: the region's integer system (586 functions) from the region hub, in quad (quad_eval, region=A);
  B  directly from the decay: the decay system integrated in quad from P0 to the physical point along a complex path
     on which the two variables that become negative cross zero with a positive imaginary part (no region system);
  C  the earlier, independent continuation through the corner (iπ at the corner; the region heads *_edgeP0, double).
usage: region_continuation_check.py [X,Y]  (default 0.2,0.5; all three regions)"""
import itertools, os, sys
import numpy as np
import mpmath as mp
import gen_eval, quad_data as Q, quad_eval as QE, remainder_eval as R, remainder_quad as RQ
from ddnp import DD, CD

HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 40
P = tuple(float(a) for a in sys.argv[1].split(',')) if len(sys.argv) > 1 else (0.2, 0.5)
P0 = gen_eval.P0
grid = [round(float(v), 2) for v in np.arange(-0.3, 0.31, 0.05) if abs(v) > 1e-9]


def direct(region):
    """B: the 586 functions at the physical image of P, from the decay P0 along an i0 path (quad)"""
    X, Y = mp.mpf(P[0]), mp.mpf(P[1])
    Z = 1 - X - Y
    xt, yt = {'z': (-X / Z, -Y / Z), 'x': (1 / Z, -X / Z), 'y': (-X / Z, 1 / Z)}[region]
    zt = 1 - xt - yt
    neg = [v for v, a, b in (('x', P0[0], xt), ('y', P0[1], yt), ('z', 1 - P0[0] - P0[1], zt)) if a > 0 > b]
    def side_ok(d):
        ix, iy = d[0] * float(xt - P0[0]), d[1] * float(yt - P0[1])
        return all(dict(x=ix, y=iy, z=-(ix + iy))[v] > 0.02 for v in neg)
    cand = sorted(((Q.path_margin(P0, (float(xt), float(yt)), d)[0], d) for d in itertools.product(grid, grid) if side_ok(d)),
                  reverse=True)
    for m, d in cand[:10]:
        try:
            (inf, F), = Q.run([(xt, yt)], os.path.join(HERE, 'out', 'iter_quad_P0.bin'), 'quad', 1e-28, d, kmax=12)
            return F, d, inf
        except RuntimeError:
            continue
    raise RuntimeError('no path')


cols = [k for k in R.KEYS]
layers = sorted({k for _, k in cols}, key=lambda s: [int(v) for v in __import__('re').findall(r'-?\d+', s)], reverse=True)
for region in ('x', 'y', 'z'):
    A = QE.evaluate_remainders_quad([P], 1e-28, region=region)[0]
    F, d, inf = direct(region)
    Fcd = CD.from_mp(F)
    Fcd = CD(DD(Fcd.re.hi[None, :], Fcd.re.lo[None, :]), DD(Fcd.im.hi[None, :], Fcd.im.lo[None, :]))
    B = RQ.remainders_from_dd([P], RQ.generators_dd(Fcd), region=region)
    Tc, _ = gen_eval.evaluate_direct([P], 1e-12, data=os.path.join(HERE, 'out', 'region_%s_bs_edgeP0.bin' % region), region=region)
    C = R.remainders_from(R.region_to_phys(region, [P]), Tc, rsign=-1)[0]
    a = A.to_mp()[0]; b = B.to_mp()[0]; a_c = A.to_complex()[0]
    x, y = R.region_to_phys(region, [P])[0]
    print('region %s, (X, Y) = %s -> (x, y, z) = (%.4g, %.4g, %.4g); direct path from P0: deformation %s, %d steps'
          % (region, P, x, y, 1 - x - y, d, inf['steps']))
    print('   %-17s %-24s %-24s %-24s' % ('colour layer', 'A vs B: Re / Im (quad)', 'A vs C: Re / Im (double)', '|Im| / |value|'))
    for lay in layers:
        js = [j for j, (n, k) in enumerate(cols) if k == lay]
        rb = max(float(abs(mp.re(a[j] - b[j])) / abs(a[j])) for j in js); ib = max(float(abs(mp.im(a[j] - b[j])) / abs(a[j])) for j in js)
        rc = max(abs((a_c[j] - C[j]).real) / abs(a_c[j]) for j in js); ic = max(abs((a_c[j] - C[j]).imag) / abs(a_c[j]) for j in js)
        im = min(abs(a_c[j].imag) / abs(a_c[j]) for j in js)
        print('   %-17s %9.1e / %-9.1e     %9.1e / %-9.1e     min %.2f' % (lay, rb, ib, rc, ic, im))
