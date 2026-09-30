#!/bin/bash
# Exhaustive S23 census, n = $1..$2, 6 gentreeg res/mod shards per n.
cd "$(dirname "$0")"
for n in $(seq $1 $2); do
  t0=$(date +%s)
  for r in 0 1 2 3 4 5; do
    ( /opt/homebrew/bin/gentreeg -p -q $n $r/6 | ./s23_census $n > raw/n${n}_s${r}of6.txt ) &
  done
  wait
  echo "n=$n done in $(( $(date +%s) - t0 ))s: $(grep -h STATS raw/n${n}_s*of6.txt | awk '{for(i=1;i<=NF;i++){split($i,a,"="); s[a[1]]+=a[2]}} END{printf "trees=%d S23viol=%d zero=%d PVviol=%d", s["trees"], s["r1_viol"], s["r1_zero"], s["pv_viol"]}')" >> s23_progress.log
done
echo ALLDONE >> s23_progress.log
