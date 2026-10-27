#!/usr/bin/env python3
"""
validate_tridents.py -- exact comparison of the trident series of
trident_rotation.py (blow-up variables) with the published physical-variable
series (trident_c<c>_w<w>.json, order 20), in the window the blow-up order
determines: SMALL <= ORDt, BIG <= ORDy - 2 ORDt.

    python3 validate_tridents.py --tridents ../trident/test_5_15 --ref ../../data/series_trident_all_weights \
        --ordt 5 --ordy 15
"""
import argparse
import json
import os
import sys
import time
from types import SimpleNamespace

from gmpy2 import mpq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import phys
import pipe_common as pc

ap = argparse.ArgumentParser()
ap.add_argument('--tridents', required=True)
ap.add_argument('--ref', required=True)
ap.add_argument('--ordt', type=int, required=True)
ap.add_argument('--ordy', type=int, required=True)
ap.add_argument('-c', '--corners', default='1,2,3')
ap.add_argument('-w', '--weights', default='0,1,2,3,4,5,6')
a = ap.parse_args()
X, Y = a.ordt, a.ordy - 2 * a.ordt
CONST = ('Pi', 'Zeta[3]', 'Zeta[5]')

allok = True
for c in [int(s) for s in a.corners.split(',')]:
    sm, bg = phys.MAP[c]
    for w in [int(s) for s in a.weights.split(',')]:
        t0 = time.time()
        mine = phys.slots_to_phys(SimpleNamespace(c=c), pc.load_pickle(
            os.path.join(a.tridents, 'trident_c%d_w%d.pkl' % (c, w))), X, Y)
        d = json.load(open(os.path.join(a.ref, 'trident_c%d_w%d.json' % (c, w))))
        u, v = d['variables']
        tokens, half = d['tokens'], [p == 'sqrt_uv' for p in d['function_prefactors']]
        ref = {}
        for j, row in enumerate(d['rows']):
            h = 1 if half[j] else 0
            for tid, n, den in row:
                pu, pv, lu, lv, npi, z3, z5 = tokens[tid]
                e = {u: (2 * pu + h, lu), v: (2 * pv + h, lv)}
                (ex2, lx), (ey2, ly) = e[sm], e[bg]
                if ex2 > 2 * X or ey2 > 2 * Y:
                    continue
                ck = tuple(sorted((at, 2 * k) for at, k in zip(CONST, (npi, z3, z5)) if k))
                ref[(j, ex2, ey2, lx, ly, ck)] = mpq(int(n), int(den))
        del d
        onlyr, onlym, diff = phys.compare(mine, ref)
        ok = not (onlyr or onlym or diff)
        allok &= ok
        print("c%d w%d (%s<=%d, %s<=%d): ref %d, mine %d | missing %d extra %d different %d  %s (%.0fs)"
              % (c, w, sm, X, bg, Y, len(ref), len(mine), len(onlyr), len(onlym), len(diff),
                 "OK" if ok else "MISMATCH", time.time() - t0), flush=True)
        if not ok:
            bad = sorted({k[0] + 1 for k in onlyr + onlym + diff})
            print("     tridents with differences (%d): %s" % (len(bad), bad[:40]), flush=True)
print("ALL IDENTICAL" if allok else "DIFFERENCES FOUND")
sys.exit(0 if allok else 1)
