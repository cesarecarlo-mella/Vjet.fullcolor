#!/usr/bin/env python3
"""precision of the double evaluator on the LHC domain of the production regions, made larger (--closer, default 2) towards the
edges of the triangle (u, v): v = 1/A (inverse energy), u = -b/A, w = -c/A = 1 - u - v, pT^2 / M_Z^2 = u w / v.

    LHC domain:      ptmin < pT < ptmax,  sqrt(s_hat) <= sqrts       (u w >= k1 v,  u w <= k2 v,  v >= M_Z^2 / sqrts^2)
                     and (--etamax) the jet rapidity in the partonic frame |eta*| <= etamax: u / w = exp(2 eta*), so
                     exp(-2 etamax) <= u / w <= exp(2 etamax) (at the high energies where precision matters the boost
                     of the partonic system is small, and eta* is close to the lab rapidity)
    checked domain:  u w >= k1/c v,  u w <= c k2 v,  v >= M_Z^2 / sqrts^2 / c      (c times closer to the edges)

Points: v log-spaced from its smallest to its largest value; at each v, u log-spaced from the lower root of
u (1 - u - v) = k1/c v to the middle (1 - v)/2, and mirrored (w the same): dense towards both collinear edges.
Stage 1 (fast): every point in double with the error estimate (an upper bound, ~1000 times conservative).
Stage 2 (--quad): every point whose estimate exceeds the target in quad (each point on its own; a point quad cannot
do is reported), which gives the actual error of double.  Results: accuracy_study/pheno_<region>.npz.

usage: pheno_domain_check.py [--regions z,x,y] [--ptmin 10] [--ptmax 3000] [--sqrts 13600] [--target 1e-6]
                             [--nv 15] [--nu 12] [--quad]"""
import argparse, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'iter int ', 'generators_deq'))
import remainder_eval as R    # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('--regions', default='z,x,y')
ap.add_argument('--ptmin', type=float, default=10.0)
ap.add_argument('--ptmax', type=float, default=3000.0)
ap.add_argument('--sqrts', type=float, default=13600.0)
ap.add_argument('--target', type=float, default=1e-6)
ap.add_argument('--closer', type=float, default=2.0, help='c: how many times closer to the edges than the LHC domain')
ap.add_argument('--etamax', type=float, default=None, help='jet rapidity acceptance |eta*| <= etamax (none by default)')
ap.add_argument('--nv', type=int, default=15)
ap.add_argument('--nu', type=int, default=12)
ap.add_argument('--quad', action='store_true', help='stage 2: the actual error, in quad, where the estimate exceeds the target')
ap.add_argument('--per-bin', type=int, default=0, help='stage 2 only at the k worst points (by estimate) of each min(u, w) bin')
a = ap.parse_args()
MZ = 91.1876
c = a.closer
k1, k2 = (a.ptmin / MZ) ** 2 / c, (a.ptmax / MZ) ** 2 * c
vlo = (MZ / a.sqrts) ** 2 / c
vhi = 1 + 2 * k1 - 2 * np.sqrt(k1 + k1 ** 2)
P = []
for v in np.geomspace(vlo, vhi * (1 - 1e-6), a.nv):
    disc = np.sqrt(max((1 - v) ** 2 - 4 * k1 * v, 0))
    ulo = (1 - v - disc) / 2
    if a.etamax is not None:                                   # u / w >= exp(-2 etamax), w = 1 - u - v
        r = np.exp(-2 * a.etamax)
        ulo = max(ulo, r * (1 - v) / (1 + r))
    us = np.geomspace(ulo, (1 - v) / 2, a.nu)
    for u in np.concatenate([us, (1 - v - us)[:-1]]):
        if u * (1 - u - v) <= k2 * v:
            P.append((u, v))
P = np.array(P)
U, V = P[:, 0], P[:, 1]
W = 1 - U - V
PT = MZ * np.sqrt(U * W / V)
print('checked domain: %.2f GeV < pT < %.0f GeV, v >= %.1e (sqrt(s) <= %.0f GeV), v <= %.3f%s; %d points per region'
      % (a.ptmin / np.sqrt(c), a.ptmax * np.sqrt(c), vlo, MZ / np.sqrt(vlo), vhi,
         ', |eta*| <= %g' % a.etamax if a.etamax is not None else '', len(U)))
print('   closest to the edges: v >= %.1e, u >= %.1e, w >= %.1e, 1 - v >= %.3f' % (V.min(), U.min(), W.min(), 1 - V.max()), flush=True)
for A in a.regions.split(','):
    out = os.path.join(HERE, 'pheno_%s%s.npz' % (A, '_eta%g' % a.etamax if a.etamax is not None else ''))
    pts = [tuple(p) for p in np.stack([W, U], -1)]                # the region triangle (X, Y) = (w, u)
    if os.path.exists(out) and np.load(out)['u'].shape == U.shape and np.allclose(np.load(out)['u'], U):
        d = dict(np.load(out))
    else:
        t0 = time.time()
        Rd, info, ed = R.evaluate_remainders(pts, 1e-12, error=True, region=A)
        d = dict(u=U, v=V, values=Rd, est=ed / np.maximum(np.abs(Rd), 1e-300), actual=np.full(Rd.shape, np.nan))
        np.savez(out, **d)
        print('region %s: double %.0f s' % (A, time.time() - t0), flush=True)
    est = d['est']
    sus = np.where(est.max(1) > a.target)[0]
    mind = np.minimum(U, W)
    print('region %s: estimate (84 values x %d points): median %.1e, 99%% %.1e, worst %.1e; points with an estimate above '
          '%.0e: %d' % (A, len(pts), np.median(est), np.percentile(est, 99), est.max(), a.target, len(sus)))
    for lo, hi in ((1e-2, 1), (1e-3, 1e-2), (1e-4, 1e-3), (1e-5, 1e-4), (1e-6, 1e-5), (0, 1e-6)):
        m = (mind >= lo) & (mind < hi)
        if m.any():
            print('   min(u, w) in [%.0e, %.0e): %3d points, worst estimate per point: median %.1e, max %.1e'
                  % (lo, hi, m.sum(), np.median(est[m].max(1)), est[m].max()))
    if not a.quad:
        continue
    import quad_eval as QE
    BINS = ((1e-2, 1), (1e-3, 1e-2), (1e-4, 1e-3), (1e-5, 1e-4), (1e-6, 1e-5), (0, 1e-6))
    if a.per_bin:
        pick = []
        for lo, hi in BINS:
            b = [i for i in sus if lo <= mind[i] < hi]
            pick += sorted(b, key=lambda i: -est[i].max())[:a.per_bin]
        sus = np.array(sorted(pick), int)
    t0 = time.time()
    for n, i in enumerate(sus):
        if np.isfinite(d['actual'][i]).all():
            continue
        try:
            Rq, iq = QE.evaluate_remainders_quad([pts[i]], 1e-28, region=A)
            vq = Rq.to_complex()[0]
            d['actual'][i] = np.abs(d['values'][i] - vq) / np.maximum(np.abs(vq), 1e-300)
        except Exception as e:
            print('   quad failed at (u, v) = (%.2e, %.2e): %s' % (U[i], V[i], str(e)[:80]), flush=True)
        if n % 10 == 9:
            np.savez(out, **d)
            print('   quad: %d of %d points, %.0f s' % (n + 1, len(sus), time.time() - t0), flush=True)
    np.savez(out, **d)
    act = d['actual'][sus]
    ok = np.isfinite(act).all(1)
    worst = np.nanmax(act) if ok.any() else np.nan
    print('region %s: actual error of double where the estimate exceeds %.0e: worst %.1e; values above the target: %d of %d'
          % (A, a.target, worst, int((act[ok] > a.target).sum()), act[ok].size))
    for lo, hi in ((1e-2, 1), (1e-3, 1e-2), (1e-4, 1e-3), (1e-5, 1e-4), (1e-6, 1e-5), (0, 1e-6)):
        m = (mind[sus] >= lo) & (mind[sus] < hi) & ok
        if m.any():
            print('   min(u, w) in [%.0e, %.0e): actual worst %.1e' % (lo, hi, act[m].max()))
    print('region %s: %s' % (A, 'double better than %.0e everywhere' % a.target if worst <= a.target else
                             'double NOT better than %.0e everywhere' % a.target), flush=True)
