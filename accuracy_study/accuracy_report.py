#!/usr/bin/env python3
"""summary of an accuracy study (accuracy_job.py chunks in DIR), per colour layer:
  double:  the relative error estimate err_d / |v| (median, 99th percentile);
  quad (if present), the reference v_q with its own estimate err_q:
    true error       |v_d - v_q| / |v|  (median, worst);
    covered          |v_d - v_q| <= err_d               (the double estimate bounds the actual error);
    contained        |v_d - v_q| + err_q <= err_d       (the quad interval lies inside the double one);
    estimate / true  median of err_d / |v_d - v_q|      (how conservative the estimate is).
Also the same by distance d = min(x, y, z) to the edges, and the values not covered.
usage: accuracy_report.py DIR"""
import glob, os, sys
import numpy as np

D = sys.argv[1]
files = sorted(glob.glob(os.path.join(D, 'chunk_*.npz')), key=lambda f: int(f.split('_')[-1][:-4]))
C = [np.load(f) for f in files]
cat = lambda k: np.concatenate([c[k] for c in C])
P, vd, ed = cat('points'), cat('vd'), cat('ed')
layers, coefs = list(C[0]['layers']), list(C[0]['coefficients'])
quad = all('vq' in c for c in C)
npts = len(P)
print('%s: %d points (%d chunks), package %s; double %.2f s/pt%s' % (D, npts, len(C), C[0]['version'],
      np.mean([c['td'] for c in C]), ', quad %.1f s/pt' % np.mean([c['tq'] for c in C]) if quad else ''))
if quad:
    vq, vql, eq = cat('vq'), cat('vq_lo'), cat('eq')
    scale = np.maximum(np.abs(vq), 1e-300)
    true = np.abs((vd - vq) - vql)
    cov, con = true <= ed, true + eq <= ed
else:
    scale = np.maximum(np.abs(vd), 1e-300)

pc = lambda a, q: np.percentile(a, q)
hdr = '%-18s %13s %13s' % ('layer', 'est. median', 'est. 99%')
if quad:
    hdr += ' %12s %12s %9s %10s %11s %12s' % ('true median', 'true worst', 'covered', 'contained', 'est./true', 'quad est.max')
print('\nrelative values, over %d points x 6 coefficients per layer' % npts)
print(hdr)
for l, name in enumerate(layers + ['all']):
    s = (slice(None), slice(None), l) if name != 'all' else (slice(None),) * 3
    rel_e = (ed[s] / scale[s]).ravel()
    row = '%-18s %13.1e %13.1e' % (name, np.median(rel_e), pc(rel_e, 99))
    if quad:
        rt = (true[s] / scale[s]).ravel()
        ratio = (ed[s] / np.maximum(true[s], 1e-300)).ravel()
        row += ' %12.1e %12.1e %8.2f%% %9.2f%% %11.1f %12.1e' % (np.median(rt), rt.max(), 100 * cov[s].mean(),
                                                               100 * con[s].mean(), np.median(ratio), (eq[s] / scale[s]).max())
    print(row)

d = P.min(1)
d = np.minimum(d, 1 - P.sum(1))
print('\nby distance d = min(x, y, z) to the edges')
for lo, hi in ((0.1, 1), (0.03, 0.1), (0.01, 0.03), (0.003, 0.01), (0, 0.003)):
    m = (d >= lo) & (d < hi)
    if not m.any():
        continue
    rel_e = (ed[m] / scale[m]).ravel()
    row = '  %5.3f <= d < %-5g %5d points: estimate median %.0e, 99%% %.0e' % (lo, hi, m.sum(), np.median(rel_e), pc(rel_e, 99))
    if quad:
        rt = (true[m] / scale[m]).ravel()
        row += '; true median %.0e, worst %.0e; covered %.2f%%, contained %.2f%%' % (np.median(rt), rt.max(), 100 * cov[m].mean(), 100 * con[m].mean())
    print(row)

if quad:
    miss = np.argwhere(~cov)
    print('\nvalues whose actual error exceeds the double estimate: %d of %d' % (len(miss), cov.size))
    for p, i, l in miss[:30]:
        print('  (x, y) = (%.4f, %.4f)  %-6s %-18s  |v| %.2e  actual %.1e  estimate %.1e  (ratio %.1f)'
              % (P[p, 0], P[p, 1], coefs[i], layers[l], scale[p, i, l], true[p, i, l], ed[p, i, l], true[p, i, l] / ed[p, i, l]))
