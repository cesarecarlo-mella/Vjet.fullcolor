#!/usr/bin/env python3
"""Digit maps: fully expanded HA3 (C3F0Z0) around c1, c2, c3 (native and Bernoulli) vs exact GPLs.
usage: python3 plot_expanded_digits.py <dir with eval_cK.json> <outdir>"""
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
DIR, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)


def digits(s, e):
    s, e = mp.mpf(s), mp.mpf(e)
    if abs(e) < 1e-4:
        return np.nan
    d = abs(s - e)
    return 16.0 if d < 1e-16 * abs(e) else float(min(16, -mp.log10(d / abs(e))))


data, order = {}, {}
for c in ('c1', 'c2', 'c3'):
    d = json.load(open(os.path.join(DIR, 'eval_%s.json' % c)))
    order[c] = d['order']
    R = d['rows']
    xy = np.array([[r['x'], r['y']] for r in R])
    data[c] = {rep: np.array([[digits(r[rep][p], r['exact'][p]) for p in (0, 1)] for r in R])
               for rep in ('native', 'bernoulli')}


def mk(v):
    return np.c_[xy, 1 - xy.sum(1), np.ones(len(xy)), v]


with np.errstate(all='ignore'):
    for rep in ('native', 'bernoulli'):
        data.setdefault('best', {})[rep] = np.fmax(np.fmax(data['c1'][rep], data['c2'][rep]), data['c3'][rep])

cols = [('c1', 'c1only', r'corner 1, $(x,y)\to0$' + '\n' + r'$N_{\rm tot}=%d$' % order['c1']),
        ('c2', 'c2only', r'corner 2, $(y,z)\to0$' + '\n' + r'$N_{\rm tot}=%d$' % order['c2']),
        ('c3', 'c3only', r'corner 3, $(x,z)\to0$' + '\n' + r'$N_{\rm tot}=%d$' % order['c3']),
        ('best', 'coverage', 'best of the three corners\n ')]
rows = [('native', 0, 'native, Re'), ('native', 1, 'native, Im'),
        ('bernoulli', 0, 'Bernoulli, Re'), ('bernoulli', 1, 'Bernoulli, Im')]
fig, ax = plt.subplots(4, 4, figsize=(7.0, 6.6), gridspec_kw=dict(wspace=0.10, hspace=0.16))
for j, (c, mode, title) in enumerate(cols):
    for i, (rep, p, lab) in enumerate(rows):
        P.panel(ax[i, j], mk(data[c][rep][:, p]), 4, mode, show_x=(i == 3), show_y=(j == 0), soft=False)
        if j == 0:
            ax[i, 0].annotate(lab, xy=(0, 0.5), xycoords='axes fraction', xytext=(-30, 0), textcoords='offset points',
                              rotation=90, ha='center', va='center', color=P.INK, fontsize=8.5)
    ax[0, j].set_title(title, color=P.INK, pad=3, fontsize=7.5)
P.colorbar(fig, ax)
fig.suptitle(r'$N_c^3$, HA3: full expansion (Trident series and rational functions) vs exact GPLs',
             color=P.INK, fontsize=9, y=0.95)
f = os.path.join(OUT, 'Nc3_HA3_expanded_3corners_digits.pdf')
fig.savefig(f, bbox_inches='tight'); fig.savefig(f[:-4] + '.png', dpi=200, bbox_inches='tight'); plt.close(fig)
print('wrote', f)
for c, _, _ in cols:
    for rep, p, lab in rows:
        v = data[c][rep][:, p]
        print('%-4s %-14s median %4.1f   >=8: %3.0f%%   >=12: %3.0f%%   max %4.1f   worst %5.1f   chopped %d'
              % (c, lab, np.nanmedian(v), 100 * np.mean(v >= 8), 100 * np.mean(v >= 12), np.nanmax(v),
                 np.nanmin(v), int(np.isnan(v).sum())))
