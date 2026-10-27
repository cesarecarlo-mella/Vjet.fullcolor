#!/usr/bin/env python3
"""compare.py -- exact component-wise comparison of boundary-constant files (lists of Pi/Zeta/I constants).

    python3 compare.py MINE_PATTERN REF_PATTERN W1 W2      (patterns with {w}, e.g. ../families/H/inputs/bdy_RZ_yv_w{w}.m)
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'code'))
import mparse as mp

mine, ref, w1, w2 = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
ok = True
for w in range(w1, w2 + 1):
    a, b = mp.parse_file(mine.format(w=w)), mp.parse_file(ref.format(w=w))
    if len(a) != len(b):
        print('weight %d: lengths differ %d %d' % (w, len(a), len(b))); ok = False; continue
    bad = [i + 1 for i, (x, y) in enumerate(zip(a, b)) if mp.padd(x, y, -1)]
    ok &= not bad
    print('weight %d: %d of %d components differ%s' % (w, len(bad), len(a), (', e.g. %s' % bad[:5]) if bad else ''))
print('IDENTICAL' if ok else 'DIFFERENT')
