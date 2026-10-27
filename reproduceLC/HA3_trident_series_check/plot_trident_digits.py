#!/usr/bin/env python3
"""Digit maps for C3F0Z0 HA3 with Trident functions replaced by their corner series (c1, c2, c3),
versus the exact GPL value; Re and Im separately.  Style: plot_w6_paper.py.

usage: python3 plot_trident_digits.py <resdir> <outdir> <ha_series_digits_order18.json>
"""
import json
import os
import sys

import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plot_w6_paper as P
P.CLIP = True
mp.mp.dps = 30
RES, OUT, HASER = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(OUT, exist_ok=True)
CHOP = 1e-4


def digits(s, e):
    s, e = mp.mpf(s), mp.mpf(e)
    if abs(e) < CHOP:
        return np.nan
    d = abs(s - e)
    if d < 1e-16 * abs(e):
        return 16.0
    return float(min(16, -mp.log10(d / abs(e))))


def load(corner, N):
    R = json.load(open(os.path.join(RES, 'res_%s_N%d.json' % (corner, N))))['rows']
    xy = np.array([[r['x'], r['y']] for r in R])
    D = np.array([[digits(r['series'][p], r['exact'][p]) for p in (0, 1)] for r in R])
    C = np.array([r['cancel'] for r in R])
    return xy, D, C


def mk(xy, v):
    return np.c_[xy, 1 - xy.sum(1), np.ones(len(xy)), v]


def haser():
    R = json.load(open(HASER))['rows']
    return np.array([r['native'][2] for r in R], float)    # HA3, (npts, 2)


def figure(N, with_ref):
    cols = []
    for c in ('c1', 'c2', 'c3'):
        xy, D, _ = load(c, N)
        cols.append((D, '%sonly' % c, {'c1': r'corner 1, $(x,y)\to0$', 'c2': r'corner 2, $(y,z)\to0$',
                                         'c3': r'corner 3, $(x,z)\to0$'}[c]))
    with np.errstate(all='ignore'):
        best = np.fmax(np.fmax(cols[0][0], cols[1][0]), cols[2][0])
    cols.append((best, 'coverage', 'best of the three corners'))
    if with_ref:
        cols.append((haser(), 'c1only', 'HA3 series (GPLs), corner 1'))
    n = len(cols)
    fig, ax = plt.subplots(2, n, figsize=(1.55 * n + 0.6, 3.3), gridspec_kw=dict(wspace=0.10, hspace=0.18))
    for c, (D, mode, title) in enumerate(cols):
        for r in (0, 1):
            P.panel(ax[r, c], mk(xy, D[:, r]), 4, mode, show_x=(r == 1), show_y=(c == 0), soft=False)
        ax[0, c].set_title(title, color=P.INK, pad=3, fontsize=7.5)
    for r, part in enumerate(('real part', 'imaginary part')):
        ax[r, 0].annotate(part, xy=(0, 0.5), xycoords='axes fraction', xytext=(-30, 0), textcoords='offset points',
                          rotation=90, ha='center', va='center', color=P.INK, fontsize=9)
    P.colorbar(fig, ax)
    fig.suptitle(r'$N_c^3$, HA3: Trident functions $\to$ corner series (order %d) vs exact' % N,
                 color=P.INK, fontsize=9, y=1.0)
    f = os.path.join(OUT, 'Nc3_HA3_trident_series_digits_N%d.pdf' % N)
    fig.savefig(f, bbox_inches='tight'); fig.savefig(f[:-4] + '.png', dpi=200, bbox_inches='tight'); plt.close(fig)
    print('wrote', f)
    for D, mode, title in cols:
        for p, pn in ((0, 'Re'), (1, 'Im')):
            v = D[:, p]
            print('N=%d %-32s %s: median %.1f, >=8: %3.0f%%, >=12: %3.0f%%, max %.1f, worst %.1f, chopped %d'
                  % (N, title, pn, np.nanmedian(v), 100 * np.mean(v >= 8), 100 * np.mean(v >= 12),
                     np.nanmax(v), np.nanmin(v), int(np.isnan(v).sum())))


if __name__ == '__main__':
    Ns = [int(a) for a in (sys.argv[4].split(',') if len(sys.argv) > 4 else ['18'])]
    for N in Ns:
        figure(N, with_ref=(N == 18))
    for c in ('c1', 'c2', 'c3'):
        _, _, C = load(c, Ns[0])
        print('%s cancellation factor sum|R T| / |sum|: median Re %.1e, Im %.1e; max Re %.1e, Im %.1e'
              % (c, np.median(C[:, 0]), np.median(C[:, 1]), C[:, 0].max(), C[:, 1].max()))
