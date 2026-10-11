#!/bin/sh
cd "$(dirname "$0")"
until grep -q ALLDONE ems.log 2>/dev/null; do sleep 20; done
for E in 5 20 60 150; do for o in 1 2 3; do
  case $o in 1) IP=10.85;; 2) IP=16.4;; 3) IP=27.0;; esac
  v=dd_o${o}_E$E; python3 mkd.py $v $o $E && gfortran -O3 -march=native -fopenmp -ffixed-line-length-132 -w $v.for -o $v
  t0=$(date +%s); j=0
  for th in $(python3 nodes.py $E $IP); do j=$((j+1)); d=r_$v/n$j; mkdir -p $d; echo $th > $d/ets.in
    (cd $d && echo "1 5" | OMP_NUM_THREADS=4 ../../$v > out.txt 2>&1); done
  echo "$v done $(( $(date +%s)-t0 )) s"
done; done
echo ALLDONE
