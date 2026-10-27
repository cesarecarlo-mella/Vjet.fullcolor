#!/usr/bin/env python3
"""
plot_convergence.py -- digit maps of the five production-region series (he_corner.py outputs).

Top row: the (s, p) square of the corner u, v -> 0 (u = s p, v = (1-s) p, w = 1 - p), on a uniform
(s, p) grid; one panel per series, then the best of the five and which series wins.  The dashed lines
are the chart boundaries: s = 1/2 (u = v), p (1 + s) = 1 (u = w), p (2 - s) = 1 (v = w).
Bottom row: best of the five and which series wins, on the (u, v) triangle.

In the CC region {z > 1, 1 - z < y < 0}:  s = y/(y - 1), p = (1 - y)/z, and the region is 0 < s, p < 1.

    python3 plot_convergence.py <sp json path with CORNER> <triangle json path with CORNER> <out.png> [title]
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

SERIES = [('th', 'bern', 'threshold\n(Bernoulli w, u)'),
          ('he2', 'bern_rho', '(u,v) corner, s = 0\n(Bernoulli p)'),
          ('he2b', 'bern_rho', '(u,v) corner, s = 1\n(Bernoulli p)'),
          ('he3', 'bern_rho', '(w,v) corner, s~ = 0\n(Bernoulli p~)'),
          ('he3b', 'bern_rho', '(w,v) corner, s~ = 1\n(Bernoulli p~)')]
COLORS = ['#4C72B0', '#DD8452', '#C44E52', '#55A868', '#8172B3']


def load(path, rep):
    d = json.load(open(path))
    pts = np.array([[float(Fraction(r['u'])), float(Fraction(r['v']))] for r in d['rows']])
    dg = np.array([min([x for x in r[rep] if x is not None], default=np.nan) for r in d['rows']], float)
    return pts, dg, d


def lines(ax):
    s = np.linspace(0, 1, 200)
    ax.axvline(0.5, ls='--', lw=0.8, c='k')
    ax.plot(s, 1 / (1 + s), ls='--', lw=0.8, c='k')
    ax.plot(s, 1 / (2 - s), ls='--', lw=0.8, c='k')
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)


def main():
    sp_prefix, tri_prefix, out = sys.argv[1], sys.argv[2], sys.argv[3]
    title = sys.argv[4] if len(sys.argv) > 4 else ''
    fig, axes = plt.subplots(2, 7, figsize=(24, 7.2), gridspec_kw=dict(wspace=0.32, hspace=0.42))
    cmap = plt.get_cmap('viridis')
    norm = plt.Normalize(0, 16)
    # ---- top: (s, p) square
    D = []
    for k, (corner, rep, lab) in enumerate(SERIES):
        pts, dg, d = load(sp_prefix.replace('CORNER', corner), rep)
        K = d['K']
        u, v = pts[:, 0], pts[:, 1]
        p = u + v; s = u / p
        grid = np.full((K - 1, K - 1), np.nan)
        for (si, pi, x) in zip(np.rint(s * K).astype(int), np.rint(p * K).astype(int), dg):
            grid[pi - 1, si - 1] = x
        D.append(grid)
        ax = axes[0, k]
        ax.imshow(np.clip(grid, 0, 16), origin='lower', extent=(0.5 / K, 1 - 0.5 / K, 0.5 / K, 1 - 0.5 / K),
                  cmap=cmap, norm=norm, aspect='equal')
        lines(ax)
        ax.set_title(lab, fontsize=9)
        ax.set_xlabel('s = u/(u+v)'); ax.set_ylabel('p = u+v')
    D = np.array(D)
    best = np.nanmax(D, axis=0); which = np.nanargmax(np.nan_to_num(D, nan=-99), axis=0)
    ax = axes[0, 5]
    im = ax.imshow(np.clip(best, 0, 16), origin='lower', extent=(0.5 / K, 1 - 0.5 / K, 0.5 / K, 1 - 0.5 / K), cmap=cmap, norm=norm, aspect='equal')
    cs = ax.contour(np.arange(1, K) / K, np.arange(1, K) / K, best,
                    levels=[4, 8, 12], colors='w', linewidths=0.7)
    ax.clabel(cs, fmt='%d', fontsize=7)
    lines(ax); ax.set_title('best of the five\n(min %.1f, median %.1f digits)' % (np.nanmin(best), np.nanmedian(best)), fontsize=9)
    ax.set_xlabel('s'); ax.set_ylabel('p')
    fig.colorbar(im, ax=axes[0, :6], fraction=0.012, pad=0.01, label='correct digits (min of Re, Im)')
    ax = axes[0, 6]
    ax.imshow(which, origin='lower', extent=(0.5 / K, 1 - 0.5 / K, 0.5 / K, 1 - 0.5 / K), cmap=ListedColormap(COLORS),
              norm=BoundaryNorm(np.arange(-0.5, 5.5), 5), aspect='equal')
    lines(ax); ax.set_title('which series wins', fontsize=9); ax.set_xlabel('s'); ax.set_ylabel('p')
    # ---- bottom: triangle
    T = []
    for k, (corner, rep, lab) in enumerate(SERIES):
        pts, dg, d = load(tri_prefix.replace('CORNER', corner), rep)
        T.append(dg)
    T = np.array(T)
    tb = np.nanmax(T, axis=0); tw = np.nanargmax(np.nan_to_num(T, nan=-99), axis=0)
    tri = mtri.Triangulation(pts[:, 0], pts[:, 1])
    for ax, data, ttl, cm, nm in ((axes[1, 5], np.clip(tb, 0, 16), 'best of the five (triangle)', cmap, norm),
                                  (axes[1, 6], tw, 'which series wins (triangle)', ListedColormap(COLORS), BoundaryNorm(np.arange(-0.5, 5.5), 5))):
        ax.tripcolor(tri, data, cmap=cm, norm=nm, shading='gouraud' if ax is axes[1, 5] else 'flat')
        ax.plot([0, 1, 0, 0], [0, 0, 1, 0], c='k', lw=0.8)
        ax.plot([0, 1 / 3], [0, 1 / 3], ls='--', lw=0.8, c='k')          # u = v
        ax.plot([0, 1 / 2], [1, 0], ls='--', lw=0.8, c='k')              # u = w
        ax.plot([1, 0], [0, 1 / 2], ls='--', lw=0.8, c='k')              # v = w
        ax.set_aspect('equal'); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_xlabel('u'); ax.set_ylabel('v'); ax.set_title(ttl, fontsize=9)
    for k, (corner, rep, lab) in enumerate(SERIES):
        ax = axes[1, k]
        ax.tripcolor(tri, np.clip(T[k], 0, 16), cmap=cmap, norm=norm, shading='gouraud')
        ax.plot([0, 1, 0, 0], [0, 0, 1, 0], c='k', lw=0.8)
        ax.set_aspect('equal'); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_xlabel('u'); ax.set_ylabel('v'); ax.set_title(lab.replace('\n', ' '), fontsize=8)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in COLORS]
    fig.legend(handles, [l.replace('\n', ' ') for _, _, l in SERIES], loc='lower center', ncol=5, fontsize=9, frameon=False)
    fig.suptitle(title, fontsize=12)
    fig.savefig(out, dpi=130, bbox_inches='tight')
    fig.savefig(out.replace('.png', '.pdf'), bbox_inches='tight')
    print('wrote', out)


if __name__ == '__main__':
    main()
