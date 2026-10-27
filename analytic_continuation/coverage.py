#!/usr/bin/env python3
"""
coverage.py -- coverage of the production triangle by a set of corner expansions (he_corner.py
outputs on the same lattice): at every point the best of the chosen representations is kept.

    python3 coverage.py <json>:<rep> [<json>:<rep> ...]
    e.g.  coverage.py th.json:bern he2.json:square he3.json:square

Prints, for every input and for the combination, the min(Re, Im) digits statistics over the
triangle, and the map of the combination (rows v = 1 -> 0, i.e. threshold corner at the top,
columns u = 0 -> 1; 0-9 digits, '-' below 1).
"""
import json
import sys
from fractions import Fraction

import numpy as np

sets = []
for arg in sys.argv[1:]:
    path, rep = arg.rsplit(':', 1)
    d = json.load(open(path))
    pts = [(Fraction(r['u']), Fraction(r['v'])) for r in d['rows']]
    dg = []
    for r in d['rows']:
        vals = [x for x in r[rep] if x is not None]
        dg.append(min(vals) if vals else np.nan)
    sets.append(('%s %s:%s' % (d['name'], d.get('corner', ''), rep), d['K'], pts, np.array(dg, float)))
K, pts = sets[0][1], sets[0][2]
assert all(s[2] == pts for s in sets), 'different lattices'


def stats(label, dg):
    print("  %-40s median %4.1f | >=4: %3d%% | >=8: %3d%% | >=12: %3d%% | worst %5.1f"
          % (label, np.nanmedian(dg), 100 * np.mean(dg >= 4), 100 * np.mean(dg >= 8), 100 * np.mean(dg >= 12), np.nanmin(dg)))


best = np.fmax.reduce([s[3] for s in sets])
which = np.argmax(np.nan_to_num(np.array([s[3] for s in sets]), nan=-99), axis=0)
print("coverage over %d points of the triangle (K = %d):" % (len(pts), K))
for s in sets:
    stats(s[0], s[3])
stats('BEST OF THE ABOVE', best)
M = {(int(u * K), int(v * K)): (best[k], which[k]) for k, (u, v) in enumerate(pts)}
print("\nmap of the best digits (rows v = %d/%d .. 1/%d, columns u = 1/%d ..):" % (K - 1, K, K, K))
for j in range(K - 1, 0, -1):
    line = ''.join(('-' if M[(i, j)][0] < 1 else str(min(9, int(M[(i, j)][0])))) for i in range(1, K) if (i, j) in M)
    print('  v=%5.2f |%s' % (j / K, line))
print("\nwhich corner is best (index in the argument list):")
for j in range(K - 1, 0, -1):
    print('  v=%5.2f |%s' % (j / K, ''.join(str(M[(i, j)][1]) for i in range(1, K) if (i, j) in M)))
