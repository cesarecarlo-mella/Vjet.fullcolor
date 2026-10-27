#!/bin/sh
# GPL form of every HA file built only from letters 1..6 (3 jobs in parallel, low priority)
cd "$(dirname "$0")"
SRC=../../remainder_denominator_bracket
for L in C1F2Z0 C9F2Z0 C1F1Z1 C9F1Z1 C2F1Z0 C3F0Z0 C8F1Z0; do
  for i in 1 2 3 4 5 6; do
    echo "nice -n 10 python3 ii2g.py $SRC/${L}_HA$i.wl ${L}_HA${i}_G.m 'HAG[\"$L\", $i]'"
  done
done | xargs -P 3 -I{} sh -c '{}' > run_all.log 2>&1
echo ALLDONE >> run_all.log
