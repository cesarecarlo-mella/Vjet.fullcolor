#!/bin/bash
# The series of all 84 finite remainders at order N = 20 (cluster): the ITER family at 24/72 in the charts 1, 3, 5,
# then the assembly into remainder_series_N20/<coefficient>/<layer>/series{c1,c2,c3}.m, seriesBer{c1,c2,c3}.m.
# Resumable: rerun the same line and finished charts / layers are skipped.
#
# One line, from the repository root (optional argument: cores, default 32):
#   nohup bash blowup_series/run_iter_N20.sh 32 > blowup_series/iter_N20.log 2>&1 &
# Progress:
#   tail -5 blowup_series/iter_N20.log
set -e
cd "$(dirname "$0")"
NC=${1:-32}
PY=venv/bin/python
mkdir -p families/ITER/dist_24_72
echo "[$(date +%T)] step 1: ITER family at 24/72, charts 1, 3, 5, $NC cores"
$PY code/run_array.py -w 6 -n "$NC" -c 1,3,5 --ordt 24 --ordy 72 --workdir families/ITER/dist_24_72 \
    --inputs families/ITER/inputs --letters letters --clean-residues
echo "[$(date +%T)] step 2: the 84 remainders at the three corners"
$PY code/assemble_remainders.py --iter families/ITER/dist_24_72 --inputs families/ITER/inputs \
    --generators "../iter int /generators_deq/out/generators.wl" --zgqq ../Zgqq_finite_remainders \
    --out ../remainder_series_N20 -n $(( NC / 2 > 1 ? NC / 2 : 1 ))
echo "[$(date +%T)] ALL DONE: $(find ../remainder_series_N20 -name 'series*.m' | wc -l) files in remainder_series_N20 (504 expected)"
