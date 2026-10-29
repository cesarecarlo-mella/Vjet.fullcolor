#!/bin/bash
# Accuracy study (cluster): N points uniform in the decay triangle, double (two paths) and full quad (two paths), with
# the package built from this checkout; then the report.  Resumable: finished chunks are skipped.
#
# One line, from the repository root (arguments: points, parallel jobs; default 1000 32):
#   nohup bash accuracy_study/run_cluster.sh 1000 32 > accuracy_study/run.log 2>&1 &
# Progress:
#   ls accuracy_study/results_1000/*.npz | wc -l ; tail -3 accuracy_study/run.log
# Cost: about 10 s per point (quad, both paths), i.e. ~3 CPU hours for 1000 points.
set -e
cd "$(dirname "$0")/.."
N=${1:-1000}; NC=${2:-32}; CH=$(( N / 10 > 0 ? N / 10 : 1 ))
PY=${PY:-blowup_series/venv/bin/python}
$PY -c "import numpy, mpmath" || { echo "need numpy and mpmath in $PY"; exit 1; }
$PY package/build_dist.py
PKG=dist/Zgqq_evaluator-$(cat package/VERSION)
make -C "$PKG" ${CXX:+CXX=$CXX} > /dev/null
OUT=accuracy_study/results_$N
mkdir -p "$OUT"
echo "[$(date +%T)] $N points in $CH chunks, $NC at a time"
seq 0 $(( CH - 1 )) | xargs -P "$NC" -I{} sh -c "[ -f $OUT/chunk_{}.npz ] || $PY accuracy_study/accuracy_job.py --package $PKG --npts $N --chunk {} --nchunks $CH --quad --out $OUT > $OUT/chunk_{}.log 2>&1"
echo "[$(date +%T)] all chunks done"
$PY accuracy_study/accuracy_report.py "$OUT" | tee "$OUT/report.txt"
