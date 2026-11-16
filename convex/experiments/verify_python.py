#!/usr/bin/env python3
"""Independent re-check of the `verify` split.

    python3 experiments/verify_python.py results/verify.csv [path-to-fmm-binary]

For every row of the benchmark CSV this script
  * recomputes the optimum by brute force over all k-subsets (pure Python,
    exact integers) and compares it with the value reported by the C++ code;
  * re-validates the returned sets (size, distinctness, minimum distance);
  * optionally runs the general-dimension clique/matrix-product solver from
    the repository root (the "FMM" solver) and compares its optimum.
"""
import csv
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exactcheck import d2, optimum, strictly_convex_ccw  # noqa: E402
from itertools import combinations  # noqa: E402


def read(path):
    with open(path) as f:
        n, k = map(int, f.readline().split())
        return [tuple(map(int, f.readline().split())) for _ in range(n)], k


def main():
    rows = list(csv.DictReader(open(sys.argv[1])))
    fmm = sys.argv[2] if len(sys.argv) > 2 else None
    base = os.path.join(os.path.dirname(os.path.abspath(sys.argv[1])), "..")
    bad = fmm_ok = 0
    for row in rows:
        P, k = read(os.path.join(base, row["file"]))
        opt = optimum(P, k)
        ok = strictly_convex_ccw(P) and int(row["r2"]) == opt
        for col in ("set_capped", "set_certified"):
            S = list(map(int, row[col].split()))
            ok = ok and len(set(S)) == k and min(d2(P, a, b) for a, b in combinations(S, 2)) == opt
        if fmm:
            inp = f"{len(P)} 2 {k}\n" + "".join(f"{x} {y}\n" for x, y in P)
            out = subprocess.run([fmm], input=inp, capture_output=True, text=True).stdout
            val = int(out.split()[1])
            ok = ok and val == opt
            fmm_ok += val == opt
        if not ok:
            bad += 1
            print("MISMATCH", row["file"])
    print(f"{len(rows) - bad}/{len(rows)} instances agree with Python brute force"
          + (f"; FMM solver agrees on {fmm_ok}/{len(rows)}" if fmm else ""))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
