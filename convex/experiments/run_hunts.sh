#!/bin/sh
# Randomized searches for counterexamples, as reported in the paper.
#   v5_*   V5: alternating greedy seeded with a tight pair     (n <= 22, k <= 9)
#   c1_*   C1: tight-pair conjecture, exact DP-based check
#   v4_*   V4: tight-pair seed, every left/right pattern (DFS)
#   *p_*   perturbations of all known counterexamples
# Logs go to convex/results/hunts/; hunt_macros.py writes the paper macros.
# (c1_21..24 used 500 instances with n <= 30 and were stopped at a time limit;
#  their last progress line is what is counted.)
set -e
cd "$(dirname "$0")/.."
mkdir -p build results/hunts
g++ -O2 -std=c++17 -o build/v5check experiments/v5check.cpp convex_dispersion.cpp geometry.cpp
H=results/hunts
for s in 1 2 3 4; do
  python3 experiments/v5_hunt.py $s 1500 22 9 | ./build/v5check > $H/v5_$s.log &
done
wait
for s in 21 22 23 24; do
  python3 experiments/v5_hunt.py $s 500 30 12 | ./build/v5check --c1 > $H/c1_$s.log &
done
wait
for s in 31 32; do python3 experiments/v5_hunt.py $s 3000 22 10 | ./build/v5check --v4 > $H/v4_$s.log & done
for s in 41 42; do python3 experiments/v5_hunt.py $s 3000 20 10 | ./build/v5check --c1 > $H/c1_$s.log & done
wait
for s in 51 52; do python3 experiments/v5_hunt.py $s 4000 26 12 twochain | ./build/v5check --v4 > $H/v4_$s.log & done
for s in 61 62; do python3 experiments/v5_hunt.py $s 4000 24 12 twochain | ./build/v5check --c1 > $H/c1_$s.log & done
wait
for s in 71 72; do python3 experiments/perturb_hunt.py $s 20000 | ./build/v5check --v4 > $H/v4p_$s.log & done
for s in 81 82; do python3 experiments/perturb_hunt.py $s 20000 | ./build/v5check --c1 > $H/c1p_$s.log & done
wait
for s in 91 92; do python3 experiments/v5_hunt.py $s 3000 22 10 | ./build/v5check --v4pm > $H/v4pm_$s.log & done
for s in 93 94; do python3 experiments/perturb_hunt.py $s 20000 | ./build/v5check --v4pm > $H/v4pmp_$s.log & done
wait
python3 experiments/hunt_macros.py
# The C1 counterexample (results/counterexamples/C1_n8_k4.txt) was constructed,
# not found by the hunts above:
#   python3 experiments/cap_lemma.py 4 60        # quadrilateral whose caps are all bad
#   python3 experiments/build_from_cap.py "<Q>"  # one blocking point per cap, then exact check
