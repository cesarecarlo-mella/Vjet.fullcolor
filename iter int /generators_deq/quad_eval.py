"""the 84 remainders in full double-double precision (about 32 digits, 27 at worst from the start values):
cpp/iter_quad integrates the exact integer system (586 functions) from the quad values at P0 (out/iter_quad_P0.bin,
notes/evaluator_notes.tex) to each point, on the complex-deformed path of gen_eval.choose_deformation; the generators
(remainder_quad.generators_dd) and the remainders (remainder_quad.remainders_from_dd) follow in double-double.

    from quad_eval import evaluate_remainders_quad
    R, info, err = evaluate_remainders_quad(points, error=True)     # R: ddnp.CD (npts, 84), columns remainder_eval.KEYS

Points near a spurious curve: as in double (gen_eval.plan), but with a 5-level symmetric stencil, error O(h^10), at
h = rho / 1000 (error ~1e-26); closer than NEAR_Q * rho to letters 15, 16 (also poles of the rational functions).
error=True: a second evaluation on another deformed path (and half the stencil step); ERR_FACTOR_Q times the
difference, plus ERR_EPS_Q sum |terms| of the sum over rows: the start values at P0 are common to both paths (their
error is ~3e-28 relative to the largest function), and the cancellation in the sum amplifies it."""
import os, subprocess
import numpy as np
import mpmath as mp
import gen_eval
import remainder_quad as RQ
from ddnp import DD, CD

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_Q = os.path.join(HERE, 'out', 'iter_quad_P0.bin')
BIN_Q = os.path.join(HERE, 'cpp', 'iter_quad')
NF = 586
NEAR_Q_DE, NEAR_Q_POLE = 1e-12, 1e-4             # spurious curves: poles of the DE only / also of the rational functions
KMAX = 12                                          # extrapolation levels (16 amplifies round-off ~1e5: the dd floor on long paths)
HQ = 1e-3                                           # stencil step relative to rho
LEVELS = 5
ERR_FACTOR_Q = 10
ERR_EPS_Q = 1e-27


def _weights(L):
    """w_k, k = 1..L: sum_k w_k (f(P + k h n) + f(P - k h n)) / 2 = f(P) + O(h^(2L))"""
    with mp.workdps(40):
        A = mp.matrix([[mp.mpf(k) ** (2 * j) for k in range(1, L + 1)] for j in range(L)])
        w = mp.lu_solve(A, mp.matrix([1] + [0] * (L - 1)))
        return [DD.from_mp(w[k]) for k in range(L)]


W_Q = _weights(LEVELS)


def data_q(region='decay'):
    """the iter_quad data file of the decay triangle (start at P0) or of a production region (start at its hub)"""
    return DATA_Q if region == 'decay' else os.path.join(HERE, 'out', 'region_%s_iter_quad.bin' % region)


def spurious_near_q(P, region='decay'):
    x, y = P
    rho = min(x, y, 1 - x - y)
    best, e = None, 1e-7 * rho
    pr = gen_eval.profile(region)
    for k in pr['spurious']:
        near = NEAR_Q_POLE if pr['near'][k] > 1e-6 else NEAR_Q_DE
        f = pr['sing'][k]
        g = np.array([(f(x + e, y) - f(x - e, y)) / (2 * e), (f(x, y + e) - f(x, y - e)) / (2 * e)])
        d = abs(f(x, y)) / np.linalg.norm(g)
        if d < near * rho and (best is None or d < best[0]):
            best = (d, g / np.linalg.norm(g), k)
    return best


def plan_q(points, hfac=1.0, region='decay'):
    flat, pl = [], []
    for P in points:
        near = spurious_near_q(P, region)
        if near is None:
            pl.append((len(flat), 1, None)); flat.append(P)
            continue
        d, n, k = near
        h = hfac * HQ * min(P[0], P[1], 1 - P[0] - P[1])
        pl.append((len(flat), 2 * LEVELS, dict(letter=k, distance=d, h=h)))
        flat += [(P[0] + s * j * h * n[0], P[1] + s * j * h * n[1]) for j in range(1, LEVELS + 1) for s in (1, -1)]
    return flat, pl


def _run(points, tol, delta, data):
    inp = '\n'.join('%.17g %.17g' % (p[0], p[1]) for p in points)
    r = subprocess.run([BIN_Q, data, 'quad', repr(float(tol)), repr(float(delta[0])), repr(float(delta[1])), str(KMAX)],
                       input=inp, capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError('iter_quad: ' + r.stderr[-300:])
    lines = r.stdout.split('\n')
    heads = [l.split() for l in lines if l.startswith('#')]
    vals = np.array(' '.join(l for l in lines if l and not l.startswith('#')).split(), float).reshape(len(points), NF, 4)
    info = [dict(steps=int(h[3]), evals=int(h[4]), time=float(h[5]), deformation=tuple(delta)) for h in heads]
    return vals, info


def functions_q(points, tol=1e-28, second=False, data=None, region='decay'):
    """the 586 functions at the points: CD array (npts, 586), info"""
    data = data or data_q(region)
    groups = {}
    for i, P in enumerate(points):
        groups.setdefault(gen_eval.choose_deformation(P, second, region), []).append(i)
    V = np.zeros((len(points), NF, 4)); info = [None] * len(points)
    for d, idx in groups.items():
        v, inf = _run([points[i] for i in idx], tol, d, data)
        V[idx] = v
        for i, f in zip(idx, inf):
            info[i] = f
    return CD(DD(V[..., 0], V[..., 1]), DD(V[..., 2], V[..., 3])), info


def _collapse(R, infos, pl):
    """CD (nflat, 84) -> CD (npts, 84), the stencils combined"""
    rows_re_h, rows_re_l, rows_im_h, rows_im_l, inf = [], [], [], [], []
    for i0, c, av in pl:
        if c == 1:
            v = R[i0]
            inf.append(infos[i0])
        else:
            v = None
            for j in range(LEVELS):
                s = (R[i0 + 2 * j] + R[i0 + 2 * j + 1]) * (DD.of(0.5) * W_Q[j])
                v = s if v is None else v + s
            inf.append(dict(infos[i0], time=sum(i['time'] for i in infos[i0:i0 + c]), averaged=av))
        rows_re_h.append(v.re.hi); rows_re_l.append(v.re.lo); rows_im_h.append(v.im.hi); rows_im_l.append(v.im.lo)
    return CD(DD(np.array(rows_re_h), np.array(rows_re_l)), DD(np.array(rows_im_h), np.array(rows_im_l))), inf


def _once(points, tol, second, hfac, absolute=False, region='decay'):
    """-> (R, info) and, with absolute, sum |terms| of the sum over rows (npts, 84) as well"""
    import remainder_eval
    flat, pl = plan_q(points, hfac, region)
    F, infos = functions_q(flat, tol, second, region=region)
    T = RQ.generators_dd(F)
    R, info = _collapse(RQ.remainders_from_dd(flat, T, region=region), infos, pl)
    if not absolute:
        return R, info
    if region == 'decay':
        phys, rsign = flat, 1
    else:
        phys, rsign = remainder_eval.region_to_phys(region, flat), -1
    S, _ = _collapse(CD(DD(remainder_eval.remainders_from(phys, T.to_complex(), rsign, absolute=True))), infos, pl)
    return R, info, np.abs(S.re.to_float())


def evaluate_remainders_quad(points, tol=1e-28, error=False, region='decay'):
    """-> (R, info) or (R, info, err): R a ddnp.CD array (npts, 84); err (npts, 84) float.  points: (x, y) of the decay
    triangle, or (X, Y) of the triangle of a production region (remainder_eval.REGION_MAP)"""
    pts = [tuple(map(float, p)) for p in points]
    if not error:
        return _once(pts, tol, False, 1.0, region=region)
    R, info, S = _once(pts, tol, False, 1.0, absolute=True, region=region)
    R2, info2 = _once(pts, tol, True, 0.5, region=region)
    d = (R - R2)
    err = ERR_FACTOR_Q * np.hypot(d.re.to_float(), d.im.to_float()) + ERR_EPS_Q * S
    for a, b in zip(info, info2):
        a['time_second'] = b['time']
    return R, info, err
