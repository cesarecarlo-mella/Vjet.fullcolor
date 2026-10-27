#!/bin/sh
cd "$(dirname "$0")"
SRC=../../remainder_denominator_bracket
for f in C2F1Z0_HA4 C2F1Z0_HA5 C2F1Z0_HA6 C3F0Z0_HA1 C3F0Z0_HA2 C3F0Z0_HA3 C3F0Z0_HA4 C3F0Z0_HA5 C3F0Z0_HA6 C8F1Z0_HA1 C8F1Z0_HA2 C8F1Z0_HA3 C8F1Z0_HA4 C8F1Z0_HA5 C8F1Z0_HA6; do
  L=${f%_HA*}; i=${f#*_HA}
  [ -s ${f}_G.m ] || echo "nice -n 10 python3 ii2g.py $SRC/$f.wl ${f}_G.m 'HAG[\"$L\", $i]'"
done | xargs -P 3 -I{} sh -c '{}' >> run_all.log 2>&1
echo ALLDONE >> run_all.log
