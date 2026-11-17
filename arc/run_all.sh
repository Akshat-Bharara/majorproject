#!/bin/bash
# Reproduces every experiment of NOTES.md / Section 7 of the paper.
# Correctness first; the timing runs start only if every check passes.
set -e
cd "$(dirname "$0")"
make -s all
mkdir -p data/small results
rm -f data/small/*.txt
bin/arc_tests 300 1 data results > results/tests_stdout.txt || { echo "CORRECTNESS FAILURES"; exit 1; }
python3 -I pybrute.py data results/small_cpp.txt > results/pybrute.txt || { echo "PYTHON MISMATCH"; exit 1; }
python3 -I check_hypotheses.py data/hypothesis_counterexamples.json > results/check_hypotheses.txt || exit 1
bin/arc_bench feas > results/feasibility.txt
if [ "$1" != "--no-timing" ]; then
  echo "n,k,mode,algo,rep,time_ms,decisions,r2,ok" > results/scaling.csv
  for seed in 1 2 3; do
    for n in 100 400 1000 10000; do
      for k in 4 16 64 $((n/10)); do
        bin/arc_bench run $n $k narrow 1 $seed | sed "s/,0,\([0-9.e+-]*,[0-9]*,[0-9]*,[01]\)$/,$seed,\1/" >> results/scaling.csv
      done
    done
    for n in 1000 10000 100000 1000000; do
      for k in 4 16 64 $((n/10)); do
        bin/arc_bench run $n $k wide 1 $seed | sed "s/,0,\([0-9.e+-]*,[0-9]*,[0-9]*,[01]\)$/,$seed,\1/" >> results/scaling.csv
      done
    done
  done
fi
echo done > results/.done
