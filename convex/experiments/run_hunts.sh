#!/bin/sh
# Large randomized searches for counterexamples (4 parallel seeds each).
#   V5: alternating greedy seeded with a tight pair     (n <= 22, k <= 9)
#   C1: tight-pair conjecture, exact DP-based check     (n <= 30, k <= 12)
# Logs go to convex/results/hunts/; summarise with hunt_macros.py.
set -e
cd "$(dirname "$0")/.."
mkdir -p build results/hunts
g++ -O2 -std=c++17 -o build/v5check experiments/v5check.cpp convex_dispersion.cpp geometry.cpp
for s in 1 2 3 4; do
  python3 experiments/v5_hunt.py $s 1500 22 9 | ./build/v5check > results/hunts/v5_$s.log &
done
wait
for s in 21 22 23 24; do
  python3 experiments/v5_hunt.py $s 500 30 12 | ./build/v5check --c1 > results/hunts/c1_$s.log &
done
wait
python3 experiments/hunt_macros.py
