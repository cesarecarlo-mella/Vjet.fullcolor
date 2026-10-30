#!/usr/bin/env python3
"""3D plots of one helicity coefficient, all 14 colour layers, in a production region A > 1, in the variables (u, v):

    v = 1/A (the inverse energy),  u = -b/A,  w = -c/A = 1 - u - v       (u, v, w > 0)
    region z: u = -y/z, w = -x/z;   region y: u = -z/y, w = -x/y;   region x: u = -z/x, w = -y/x

(the package's region triangle (X, Y) is (w, u), Z = v).

Default: the pheno domain, the value itself.  A cut on the transverse momentum of the outgoing parton,
pT^2 / M_Z^2 = u w / v >= kappa = (ptmin / M_Z)^2 (q^2 = M_Z^2), removes threshold (v -> 1, where u + w -> 0) and the
collinear edges u -> 0, w -> 0, which all have pT -> 0; v reaches 1 + 2 kappa - 2 Sqrt[kappa + kappa^2].  At fixed v,
u runs between the roots of u (1 - u - v) = kappa v (a grid rectangular in v and the fraction of that interval).

--log: the larger range, the (u, v) triangle up to near threshold and the collinear edges, with a logarithmic vertical
axis: the height is log10|value|, the colour its sign.  v from vlo to vhi, and at each v, u from a fraction sedge to
1 - sedge of its interval (0, 1 - v).  Every point is evaluated with its error estimate; where the double error exceeds
1e-4 of the value the point is redone in quad, and a point still above 1e-4 is left out.

    plot_region_layers.py [--region z] [--coef alpha1] [--ptmin 20] [--vmin 0.02] [--n 30]
    plot_region_layers.py --log [--vlo 0.01] [--vhi 0.95] [--sedge 0.02] [--n 25]

Writes region_<A>_<coef>_uv_{re,im}.png, resp. region_<A>_<coef>_log_{re,im}.png, next to this script (values cached
in region_<A>_<coef>_{uv,log}.npz).  Evaluator: the repository's (iter int /generators_deq), double and quad."""
import argparse, os, re, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'iter int ', 'generators_deq'))
import remainder_eval as R    # noqa: E402

MZ = 91.1876
ap = argparse.ArgumentParser()
ap.add_argument('--region', default='z', choices='xyz')
ap.add_argument('--coef', default='alpha1', choices=R.NAMES)
ap.add_argument('--ptmin', type=float, default=20.0, help='pT cut of the outgoing parton in GeV (q^2 = M_Z^2)')
ap.add_argument('--vmin', type=float, default=0.02, help='smallest v = 1/A (high energy)')
ap.add_argument('--n', type=int, default=None, help='grid points per direction (30; --log: 25)')
ap.add_argument('--log', action='store_true', help='the large range in log variables (see above)')
ap.add_argument('--vlo', type=float, default=0.01, help='--log: smallest v')
ap.add_argument('--vhi', type=float, default=0.95, help='--log: largest v')
ap.add_argument('--sedge', type=float, default=0.02, help='--log: u from sedge to 1 - sedge of its interval (0, 1 - v)')
a = ap.parse_args()
A = a.region
n = a.n or (25 if a.log else 30)
_pw = lambda k: tuple(int(v) for v in re.fullmatch(r'Nc\^(-?\d+)\*Nf\^(\d+)\*NfZ\^(\d+)', k).groups())
LAYERS = sorted([k for c, k in R.KEYS if c == a.coef], key=lambda k: (-_pw(k)[0], _pw(k)[1], _pw(k)[2]))   # Nc^3 first
cols = [R.KEYS.index((a.coef, k)) for k in LAYERS]
defs = {'z': ('u = -y/z', 'w = -x/z'), 'y': ('u = -z/y', 'w = -x/y'), 'x': ('u = -z/x', 'w = -y/x')}[A]
same = lambda c, x: c.shape == x.shape and np.allclose(c, x)


def label(k):
    c, f, z = _pw(k)
    return '$N_c^{%d}%s%s$' % (c, {0: '', 1: '\\,N_f', 2: '\\,N_f^2'}[f], '\\,N_{f,Z}' if z else '')


def evaluate(U, V, cache, careful=False):
    """the layers at (u, v); careful: with error estimates, quad where double is worse than 1e-4, NaN if still worse"""
    if os.path.exists(cache):
        d = np.load(cache)
        if same(d['u'], U) and same(d['v'], V):
            return d['values'][:, [list(d['layers']).index(k) for k in LAYERS]], d.get('how')
    pts = [tuple(p) for p in np.stack([1 - U - V, U], -1)]          # the region triangle (X, Y) = (w, u)
    if not careful:
        Val, how = R.evaluate_remainders(pts, 1e-12, region=A)[0][:, cols], None
    else:
        Rd, info, ed = R.evaluate_remainders(pts, 1e-12, error=True, region=A)
        Val, err = Rd[:, cols], ed[:, cols]
        how = np.zeros(len(pts), int)                             # 0 double, 1 quad, 2 left out
        bad = np.where((err / np.maximum(np.abs(Val), 1e-300)).max(1) > 1e-4)[0]
        print('%d of %d points: double error above 1e-4, redone in quad' % (len(bad), len(pts)), flush=True)
        if len(bad):
            import quad_eval
            Rq, iq, eq = quad_eval.evaluate_remainders_quad([pts[i] for i in bad], 1e-28, error=True, region=A)
            Vq, eq = Rq.to_complex()[:, cols], eq[:, cols]
            Val[bad] = Vq; how[bad] = 1
            still = bad[(eq / np.maximum(np.abs(Vq), 1e-300)).max(1) > 1e-4]
            Val[still] = np.nan; how[still] = 2
            print('   %d still above 1e-4 in quad: left out' % len(still), flush=True)
    np.savez(cache, u=U, v=V, values=Val, layers=np.array(LAYERS), how=how if how is not None else np.zeros(0))
    return Val, how


def figure(panels, txt, out):
    fig = plt.figure(figsize=(19, 17))
    for l, k in enumerate(LAYERS):
        ax = fig.add_subplot(4, 4, l + 1, projection='3d')
        panels(ax, l)
        ax.set_title(label(k), fontsize=13)
        ax.tick_params(labelsize=7)
    ax = fig.add_subplot(4, 4, 16); ax.axis('off')
    ax.text(0.0, 0.95, '\n'.join(txt), va='top', fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=100)
    plt.close(fig)
    print('written', out)


head = lambda part: ['%s, %s part, region %s > 1' % (a.coef, 'real' if part == 're' else 'imaginary', A), '',
                     'v = 1/%s (inverse energy),  %s,  %s = 1 - u - v' % (A, defs[0], defs[1]),
                     'v -> 1: threshold;  v -> 0: high energy', '']

if not a.log:
    kappa = (a.ptmin / MZ) ** 2
    vmax = 1 + 2 * kappa - 2 * np.sqrt(kappa + kappa ** 2)
    VV, SS = np.meshgrid(np.linspace(a.vmin, vmax * (1 - 1e-9), n), np.linspace(0.0, 1.0, n), indexing='ij')
    disc = np.sqrt(np.maximum((1 - VV) ** 2 - 4 * kappa * VV, 0.0))
    UU = (1 - VV - disc) / 2 + SS * disc                       # between the two roots of u (1 - u - v) = kappa v
    Val, _ = evaluate(UU.ravel(), VV.ravel(), os.path.join(HERE, 'region_%s_%s_uv.npz' % (A, a.coef)))
    for part, fn in (('re', np.real), ('im', np.imag)):
        def panels(ax, l):
            H = fn(Val[:, l]).reshape(n, n)
            m = max(np.abs(H).max(), 1e-300)
            ax.plot_surface(UU, VV, H, facecolors=plt.cm.coolwarm(TwoSlopeNorm(vcenter=0.0, vmin=-m, vmax=m)(H)),
                            linewidth=0, antialiased=True, shade=False)
            ax.set_xlim(0, 1); ax.set_ylim(0, vmax)
            ax.set_xlabel('u', fontsize=9); ax.set_ylabel('v', fontsize=9)
            ax.set_zlabel('%s %s' % ('Re' if part == 're' else 'Im', a.coef), fontsize=8)
        figure(panels, head(part) + ['domain: $p_T > %g$ GeV, i.e. $u\\,w/v \\geq (p_{T,min}/M_Z)^2 = %.4f$' % (a.ptmin, kappa),
                                     '    (removes threshold and the collinear edges)', '    $%g \\leq v \\leq %.3f$' % (a.vmin, vmax),
                                     'height: the value itself ($q^2 = M_Z^2 = 1$)', '',
                                     '%d points, double evaluator' % VV.size, '(start values: continued corner series + quad transport)'],
               os.path.join(HERE, 'region_%s_%s_uv_%s.png' % (A, a.coef, part)))
else:
    VV, SS = np.meshgrid(np.linspace(a.vlo, a.vhi, n), np.linspace(a.sedge, 1 - a.sedge, n), indexing='ij')
    UU = SS * (1 - VV)
    Val, how = evaluate(UU.ravel(), VV.ravel(), os.path.join(HERE, 'region_%s_%s_log.npz' % (A, a.coef)), careful=True)
    nq, nout = int((how == 1).sum()), int((how == 2).sum())
    for part, fn in (('re', np.real), ('im', np.imag)):
        def panels(ax, l):
            H = fn(Val[:, l]).reshape(n, n)
            L = np.log10(np.maximum(np.abs(H), 1e-300))
            col = np.where(H >= 0, 0.85, 0.15)                     # sign: red positive, blue negative
            fc = plt.cm.coolwarm(col)
            fc[np.isnan(H)] = (0, 0, 0, 0)
            ax.plot_surface(UU, VV, np.where(np.isnan(H), np.nan, L), facecolors=fc, linewidth=0, antialiased=True, shade=True)
            ax.set_xlim(0, 1); ax.set_ylim(0, 1)
            ax.set_xlabel('u', fontsize=9); ax.set_ylabel('v', fontsize=9)
            ax.set_zlabel('log$_{10}$ |%s %s|' % ('Re' if part == 're' else 'Im', a.coef), fontsize=8)
        figure(panels, head(part) + ['shown: %g <= v <= %g; u from %g to %g of (0, 1 - v)' % (a.vlo, a.vhi, a.sedge, 1 - a.sedge),
                                     'height: log$_{10}$ of |value| ($q^2 = 1$); colour: red > 0, blue < 0', '',
                                     '%d points: double, %d redone in quad (double error > 1e-4), %d left out' % (VV.size, nq, nout)],
               os.path.join(HERE, 'region_%s_%s_log_%s.png' % (A, a.coef, part)))
