#!/bin/sh
# After makeNS.wl and extendFibr2.wl: singlet conversion, then collect/shuffle + exact checks, one file at a time,
# low priority; then combine. Log: run_all.log
cd "$(dirname "$0")"
LCX=../../LC_Z_convert_to_II/XY_POLES
echo "$(date '+%H:%M:%S') convertPS"
nice -n 10 wolframscript -file convertPS.wl > convertPS.log 2>&1
grep -hE "colour layers|not in fibrXY|WARNING|wrote|total" convertPS.log
for i in 1 2 3 4 5 6; do
  # non-singlet: poles redone (new Delta term); eps^0 unchanged -> reuse yesterday's, re-checked against the new file
  echo "$(date '+%H:%M:%S') NS $i"
  nice -n 10 wolframscript -file polesII_polesonly.wl fc_poles_new_$i.m > FC_POLES/polesII_poles_ns_$i.log 2>&1
  grep -hE "wrote|WARNING|nonlinear" FC_POLES/polesII_poles_ns_$i.log
  for s in collected shuffled; do cp $LCX/FC_POLES/fc_poles_${i}_eps0_$s.m FC_POLES/fc_poles_new_${i}_eps0_$s.m; done
  nice -n 10 wolframscript -file checkPolesOnly.wl fc_poles_new_$i.m poles > FC_POLES/checkPoles_ns_$i.log 2>&1; cat FC_POLES/checkPoles_ns_$i.log
  nice -n 10 wolframscript -file checkPolesOnly.wl fc_poles_new_$i.m eps0 > FC_POLES/checkEps0_ns_$i.log 2>&1; cat FC_POLES/checkEps0_ns_$i.log
  # singlet
  echo "$(date '+%H:%M:%S') PS $i"
  for part in poles eps0; do
    nice -n 10 wolframscript -file polesII_$( [ $part = poles ] && echo polesonly || echo eps0 ).wl fc_ps_$i.m > FC_POLES/polesII_${part}_ps_$i.log 2>&1
    grep -hE "wrote|WARNING|nonlinear" FC_POLES/polesII_${part}_ps_$i.log
    nice -n 10 wolframscript -file checkPolesOnly.wl fc_ps_$i.m $part > FC_POLES/check_${part}_ps_$i.log 2>&1; cat FC_POLES/check_${part}_ps_$i.log
  done
done
echo "$(date '+%H:%M:%S') combine"
nice -n 10 wolframscript -file FC_SHUFFLED/combine.wl
echo "$(date '+%H:%M:%S') ALLDONE"
