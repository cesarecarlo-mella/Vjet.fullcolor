#!/usr/bin/env python3
"""
democratic.py -- fixed geometric split of the production triangle among the five series, symmetric
under u <-> w, and the resulting accuracy (he_corner.py outputs, one json per corner).

    threshold region        : threshold series in (w, u), Bernoulli in both; either a band v >= v0 ('band'),
                              or v >= lambda max(u, w) ('bary', lambda < 1; lambda = 1: the third where v is largest),
                              or v >= a max(u, w) + b min(u, w) ('cut', with the s split s0 below)
    otherwise, u <= w       : corner (u, v):  s  = u/(u+v)  <= s0  -> chart at s  = 0,  else at s  = 1
    otherwise, w <  u       : corner (w, v):  s~ = w/(w+v)  <= s0  -> chart at s~ = 0,  else at s~ = 1
                              (s0 = 1/2 for 'band' and 'bary')
    (high-energy charts with Bernoulli in p = u + v, resp. p~ = w + v)

Every chart is used inside its convergence domain: the s = 0 chart needs u < w (p (1+s) < 1), the s = 1
chart v < w (p (2-s) < 1); outside the threshold region v < lambda w (or v < v0 <= 1/3) with u <= w gives w > v.
In the (s, p) square the 'bary' cut is p = lambda/(1 - s + lambda) for u <= w and the vertical line
s = 1/(1 + lambda) for u > w.

Without a rule the band and bary families are scanned and the one with the best worst case is kept.

Metric 'part' (default): digits of the real and imaginary parts separately, each relative to itself (the
minimum of the two). Metric 'scaled': -log10 |S - E| / max(|E|, median |E| over the triangle), i.e. relative
where the function is large and absolute (in units of its typical size) where it dips; needs the 'err'
field of he_corner.py.

    python3 democratic.py <tri json with CORNER> <sp json with CORNER> <out.png> <title> [rule: band:v0, bary:lambda, cut:a,b,s0, or scan] [part|scaled]
"""
import json
import sys
from fractions import Fraction

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np

SERIES = [('th', 'bern', 'threshold (w, u), Bernoulli'),
          ('he2', 'bern_rho', '(u,v) corner, s = 0'),
          ('he2b', 'bern_rho', '(u,v) corner, s = 1'),
          ('he3', 'bern_rho', '(w,v) corner, s~ = 0'),
          ('he3b', 'bern_rho', '(w,v) corner, s~ = 1')]
COLORS = ['#4C72B0', '#DD8452', '#C44E52', '#55A868', '#8172B3']
RULES = [('band', x) for x in (0.15, 0.2, 0.25, 0.3)] + [('bary', x) for x in (0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8)]


def region(u, v, rule):
    """series index (into SERIES) for arrays u, v; w = 1 - u - v"""
    w = 1 - u - v
    kind, par = rule[0], rule[1:]
    s0 = 0.5
    if kind == 'band':
        th = v >= par[0]
    elif kind == 'bary':
        th = v >= par[0] * np.maximum(u, w)
    else:
        a, b, s0 = par
        th = v >= a * np.maximum(u, w) + b * np.minimum(u, w)
    return np.where(th, 0, np.where(u <= w, np.where(u <= s0 * (u + v), 1, 2), np.where(w <= s0 * (w + v), 3, 4)))


def rule_str(rule):
    return '%s %s' % (rule[0], ', '.join('%g' % x for x in rule[1:]))


def load(pattern, metric='part', scale=None):
    """digits of the five series at the grid points; scale = median |E| (from the first file if not given)"""
    D, pts = [], None
    for corner, rep, _ in SERIES:
        d = json.load(open(pattern.replace('CORNER', corner)))
        p = [(float(Fraction(r['u'])), float(Fraction(r['v']))) for r in d['rows']]
        if pts is not None and p != pts:
            sys.exit('different lattices')
        pts = p
        if metric == 'part':
            D.append([min([x for x in r[rep] if x is not None], default=np.nan) for r in d['rows']])
        else:
            A = np.array([abs(complex(*r['exact'])) for r in d['rows']])
            scale = np.median(A) if scale is None else scale
            err = np.array([r['err'][rep] for r in d['rows']])
            D.append(np.minimum(16, -np.log10(err / np.maximum(A, scale) + 1e-300)))
    return np.array(pts), np.array(D, float), d['N'], d['K'], scale


def apply(pts, D, rule):
    reg = region(pts[:, 0], pts[:, 1], rule)
    return reg, D[reg, np.arange(len(pts))]


def stats(x):
    return 'median %4.1f | >=4: %3d%% | >=6: %3d%% | >=8: %3d%% | >=12: %3d%% | worst %5.1f' % (
        np.nanmedian(x), 100 * np.mean(x >= 4), 100 * np.mean(x >= 6), 100 * np.mean(x >= 8), 100 * np.mean(x >= 12),
        np.nanmin(x))


def boundaries(ax, x, y, R, col):
    """region boundaries of a fine region map R (NaN outside the domain): each indicator at level 1/2"""
    for k in range(len(SERIES)):
        ind = np.where(np.isnan(R), np.nan, (R == k).astype(float))
        if np.nansum(ind) > 0:
            ax.contour(x, y, ind, levels=[0.5], colors=col, linewidths=0.9)


def main():
    tri_pat, sp_pat, out, title = sys.argv[1:5]
    metric = sys.argv[6] if len(sys.argv) > 6 else 'part'
    pts, D, N, _, scale = load(tri_pat, metric)
    print('triangle lattice, %d points, order %d, metric %s%s'
          % (len(pts), N, metric, ', median |f| = %.3g' % scale if scale else ''))
    if len(sys.argv) > 5 and sys.argv[5] != 'scan':
        kind, par = sys.argv[5].split(':')
        rule = (kind, *map(float, par.split(',')))
    else:
        best = None
        for r in RULES:
            reg, dg = apply(pts, D, r)
            print('  %-10s: %s' % (rule_str(r), stats(dg)))
            key = (np.nanmin(dg), np.nanmedian(dg))
            if best is None or key > best[0]:
                best = (key, r)
        rule = best[1]
    reg, dg = apply(pts, D, rule)
    print('rule: %s' % rule_str(rule))
    for k, (_, _, lab) in enumerate(SERIES):
        sel = reg == k
        if sel.any():
            print('  %-30s %3d points: %s' % (lab, sel.sum(), stats(dg[sel])))
    print('  %-30s %3d points: %s' % ('ALL', len(dg), stats(dg)))
    spts, SD, _, K, _ = load(sp_pat, metric, scale)
    sreg, sdg = apply(spts, SD, rule)
    print('(s, p) grid, %d points: %s' % (len(sdg), stats(sdg)))

    fig, axes = plt.subplots(1, 4, figsize=(19, 4.6), gridspec_kw=dict(wspace=0.28))
    cmap, norm = plt.get_cmap('viridis'), plt.Normalize(0, 16)
    rcmap, rnorm = ListedColormap(COLORS), BoundaryNorm(np.arange(-0.5, 5.5), 5)

    # exact split on a fine grid, triangle (u, v) and square (s, p)
    f = (np.arange(600) + 0.5) / 600                  # cell centres: p = 0 (u = v = 0) is not a point of the split
    fu, fv = np.meshgrid(f, f)
    RT = np.where(fu + fv <= 1, region(fu, fv, rule), np.nan)
    fs, fp = np.meshgrid(f, f)
    RS = region(fs * fp, (1 - fs) * fp, rule).astype(float)

    tri = mtri.Triangulation(pts[:, 0], pts[:, 1])
    im = axes[0].tripcolor(tri, np.clip(dg, 0, 16), cmap=cmap, norm=norm, shading='gouraud')
    axes[0].set_title('digits, fixed split (triangle)\nmin %.1f, median %.1f' % (np.nanmin(dg), np.nanmedian(dg)), fontsize=9)
    axes[1].imshow(RT, origin='lower', extent=(0, 1, 0, 1), cmap=rcmap, norm=rnorm, interpolation='nearest')
    axes[1].set_title('regions (triangle), %s' % rule_str(rule), fontsize=9)
    for ax, col in ((axes[0], 'w'), (axes[1], 'k')):
        boundaries(ax, f, f, RT, col)
        ax.plot([0, 1, 0, 0], [0, 0, 1, 0], c='k', lw=0.8)
        ax.set_aspect('equal'); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xlabel('u'); ax.set_ylabel('v')

    s = spts[:, 0] / (spts[:, 0] + spts[:, 1]); p = spts[:, 0] + spts[:, 1]
    G = np.full((K - 1, K - 1), np.nan)
    for si, pi, x in zip(np.rint(s * K).astype(int), np.rint(p * K).astype(int), sdg):
        G[pi - 1, si - 1] = x
    axes[2].imshow(np.clip(G, 0, 16), origin='lower', extent=(0.5 / K, 1 - 0.5 / K, 0.5 / K, 1 - 0.5 / K),
                   cmap=cmap, norm=norm, aspect='equal')
    axes[2].set_title('digits, fixed split ((s, p) square)\nmin %.1f, median %.1f' % (np.nanmin(sdg), np.nanmedian(sdg)), fontsize=9)
    axes[3].imshow(RS, origin='lower', extent=(0, 1, 0, 1), cmap=rcmap, norm=rnorm, interpolation='nearest')
    axes[3].set_title('regions ((s, p) square)', fontsize=9)
    for ax, col in ((axes[2], 'w'), (axes[3], 'k')):
        boundaries(ax, f, f, RS, col)
        ax.set_aspect('equal'); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_xlabel('s = u/(u+v)'); ax.set_ylabel('p = u+v')

    fig.colorbar(im, ax=axes.tolist(), fraction=0.012, pad=0.015,
                 label='correct digits (min of Re, Im)' if metric == 'part' else 'digits: |error| / max(|f|, median |f|)')
    fig.legend([plt.Rectangle((0, 0), 1, 1, color=c) for c in COLORS], [l for _, _, l in SERIES],
               loc='lower center', ncol=5, fontsize=9, frameon=False, bbox_to_anchor=(0.45, -0.06))
    fig.suptitle(title, fontsize=11)
    fig.savefig(out, dpi=130, bbox_inches='tight')
    fig.savefig(out.replace('.png', '.pdf'), bbox_inches='tight')
    print('wrote', out)


if __name__ == '__main__':
    main()
