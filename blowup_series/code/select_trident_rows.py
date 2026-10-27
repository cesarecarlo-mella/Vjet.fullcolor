#!/usr/bin/env python3
"""
select_trident_rows.py -- choose 1318 independent rows of the trident matrix
m and store them exactly, so that the rotation family -> trident needs only a
small file (run once, locally; needs python-flint).

The ancillary trident_basis/ gives the crossed family elements through the
1318 trident functions,

    BASISLC[F, crossing, i] = phase * sum_j m[(F, crossing, i), j] Trident[j],

with m a 23520 x 1318 matrix of real rationals.  Every crossing is available
at every corner (standard + mirror charts), so any 1318 independent rows give
the tridents:  Trident = m_sel^-1 . (BASISLC_sel / phase).
Rows are taken greedily, cheapest first (fewest nonzeros, then fewest digits):
931 tridents are single family elements, the rest small combinations.

    python3 select_trident_rows.py <trident_basis dir> ../trident/rotation_rows.json
"""
import argparse
import json
import os
import re
import sys
import time

import numpy as np
from flint import nmod_mat

ap = argparse.ArgumentParser()
ap.add_argument('basis', help='trident_basis/ directory of the ancillary files')
ap.add_argument('out', help='output JSON')
a = ap.parse_args()
P = 2**61 - 1
NT = 1318
t0 = time.time()


def read(name, mode='r'):
    with open(os.path.join(a.basis, name), mode) as f:
        return f.read()


labels = re.findall(r'BASISLC\[(\w+),\s*"(\w*)",\s*(\d+)\]', read('family_basis_labels.wl'))
phases = {(f, c, int(i)): ph for f, c, i, ph in
          re.findall(r'BASISLC\[(\w+),\s*"(\w*)",\s*(\d+)\]\s*->\s*(-?I|-?1)', read('family_basis_i_factors.wl'))}
pos = np.array(re.findall(rb'\{(\d+),(\d+)\}', read('family_to_trident_positions.wl', 'rb')), dtype=np.int64) - 1
vals = read('family_to_trident_values.wl', 'rb').replace(b'\\\n', b'').strip().lstrip(b'{').rstrip(b'}').split(b',')
vals = [v.strip().decode() for v in vals]
assert len(vals) == len(pos), "positions and values differ in length"
print("read %d labels, %d entries, %d phases (%.0fs)" % (len(labels), len(vals), len(phases), time.time() - t0), flush=True)


def mod(v):
    n, _, d = v.partition('/')
    return int(n) % P * pow(int(d or 1) % P, -1, P) % P


rows, cols = pos[:, 0], pos[:, 1]
nnz = np.bincount(rows, minlength=len(labels))
digits = np.bincount(rows, weights=np.array([len(v) for v in vals]), minlength=len(labels))
order = np.lexsort((digits, nnz))                 # cheapest rows first
rank_of = np.empty(len(labels), dtype=np.int64)
rank_of[order] = np.arange(len(labels))
AT = np.zeros((NT, len(labels)), dtype=np.uint64)
AT[cols, rank_of[rows]] = [mod(v) for v in vals]
R, rk = nmod_mat(AT.tolist(), P).rref()
if rk != NT:
    sys.exit("rank %d < %d: the rows do not determine all tridents" % (rk, NT))
piv, r = [], 0
for j in range(len(labels)):
    if r < rk and int(R[r, j]) != 0:
        piv.append(j); r += 1
sel = [int(i) for i in order[piv]]
print("selected %d rows (%.0fs)" % (len(sel), time.time() - t0), flush=True)

where = {}
for k, r in enumerate(rows):
    where.setdefault(int(r), []).append(k)
out = []
for r in sel:
    f, c, i = labels[r]
    out.append({"family": f, "crossing": c, "index": int(i), "phase": phases.get((f, c, int(i)), "1"),
                "entries": [[int(cols[k]) + 1, vals[k]] for k in where[r]]})
with open(a.out, 'w') as fh:
    json.dump({"description": "1318 independent rows of the trident matrix m (family_to_trident): "
                              "BASISLC[family, crossing, index] = phase * sum entries[j][1] * Trident[entries[j][0]]",
               "source": os.path.abspath(a.basis), "rows": out}, fh)
print("wrote %s (%.1f MB, %.0fs)" % (a.out, os.path.getsize(a.out) / 1e6, time.time() - t0))
