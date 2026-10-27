#!/usr/bin/env python3
"""Digit maps: Nc^3 HA1..6, corner-1 expansion vs exact GPLs (real and imaginary part separately)."""
import json, sys, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plot_w6_paper as P
P.CLIP = True
res, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
HA = ['$\\alpha_1$', '$\\alpha_2$', '$\\alpha_3$', '$\\gamma_1$', '$\\gamma_2$', '$\\gamma_3$']

def load(N):
    d = json.load(open(os.path.join(res, 'digits_order%d.json' % N)))
    R = d['rows']
    xy = np.array([[r['x'], r['y']] for r in R])
    nat = np.array([r['native'] for r in R], float)      # (npts, 6, 2)
    ber = np.array([r['bernoulli'] for r in R], float)
    return xy, nat, ber

def mk(xy, v):
    return np.c_[xy, 1 - xy.sum(1), np.ones(len(xy)), v]

def grid6(N, rep):
    xy, nat, ber = load(N)
    D = nat if rep == 'native' else ber
    fig, ax = plt.subplots(2, 6, figsize=(8.6, 3.2), gridspec_kw=dict(wspace=0.10, hspace=0.18))
    for c in range(6):
        for r, part in enumerate(('Re', 'Im')):
            P.panel(ax[r, c], mk(xy, D[:, c, r]), 4, 'c1only', show_x=(r == 1), show_y=(c == 0), soft=False)
            if r == 0: ax[r, c].set_title('HA%d = %s' % (c + 1, HA[c]), color=P.INK, pad=3, fontsize=9)
        for r, part in enumerate(('real part', 'imaginary part')):
            if c == 0:
                ax[r, 0].annotate(part, xy=(0, 0.5), xycoords='axes fraction', xytext=(-30, 0), textcoords='offset points',
                                  rotation=90, ha='center', va='center', color=P.INK, fontsize=9)
    P.colorbar(fig, ax)
    f = os.path.join(out, 'Nc3_HA_digits_%s_N%d.pdf' % (rep, N))
    fig.savefig(f, bbox_inches='tight'); fig.savefig(f[:-4] + '.png', dpi=200, bbox_inches='tight'); plt.close(fig)
    print('wrote', f)

def summary(Ns):
    fig, ax = plt.subplots(len(Ns), 4, figsize=(6.3, 1.6 * len(Ns) + 0.5), squeeze=False, gridspec_kw=dict(wspace=0.10, hspace=0.18))
    for r, N in enumerate(Ns):
        xy, nat, ber = load(N)
        cols = [(nat, 0, 'native, Re'), (nat, 1, 'native, Im'), (ber, 0, 'Bernoulli, Re'), (ber, 1, 'Bernoulli, Im')]
        for c, (D, p, t) in enumerate(cols):
            with np.errstate(all='ignore'):
                v = np.nanmin(D[:, :, p], axis=1)
            P.panel(ax[r, c], mk(xy, v), 4, 'c1only', show_x=(r == len(Ns) - 1), show_y=(c == 0), soft=False)
            if r == 0: ax[r, c].set_title(t, color=P.INK, pad=3, fontsize=9)
        ax[r, 0].annotate('$N=%d$' % N, xy=(0, 0.5), xycoords='axes fraction', xytext=(-30, 0), textcoords='offset points',
                          rotation=90, ha='center', va='center', color=P.INK, fontsize=9)
    P.colorbar(fig, ax)
    f = os.path.join(out, 'Nc3_HA_digits_min_summary.pdf')
    fig.savefig(f, bbox_inches='tight'); fig.savefig(f[:-4] + '.png', dpi=200, bbox_inches='tight'); plt.close(fig)
    print('wrote', f)

for N in (15, 18):
    for rep in ('native', 'bernoulli'):
        grid6(N, rep)
summary([15, 18])
# table of fractions
for N in (15, 18):
    xy, nat, ber = load(N)
    for name, D in (('native', nat), ('Bernoulli', ber)):
        for p, pn in ((0, 'Re'), (1, 'Im')):
            with np.errstate(all='ignore'):
                v = np.nanmin(D[:, :, p], axis=1)
            print('N=%d %-9s %s  min-over-HA digits: median %.1f, >=8: %3.0f%%, >=12: %3.0f%%, worst %.1f   chopped cells %d'
                  % (N, name, pn, np.nanmedian(v), 100 * np.mean(v >= 8), 100 * np.mean(v >= 12), np.nanmin(v), int(np.isnan(D[:, :, p]).sum())))
