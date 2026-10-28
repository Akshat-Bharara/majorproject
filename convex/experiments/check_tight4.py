#!/usr/bin/env python3
"""Cross-check the O(n^2 log n) k = 4 algorithm (--algo=tight4) against the
capped DP on random convex instances of all families and on every k = 4
counterexample file. usage: check_tight4.py SEED COUNT NMAX
"""

import glob
import os
import random
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from greedy_stress import family  # noqa: E402
from structure import gen         # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(os.path.dirname(HERE), "build", "convex_dispersion")


def solve(text, algo):
    out = subprocess.run([EXE, f"--algo={algo}"], input=text, capture_output=True,
                         text=True, check=True).stdout.split()
    return int(out[1]), [int(x) for x in out[3:7]]


def main():
    seed, count, nmax = map(int, sys.argv[1:4])
    rng = random.Random(seed)
    texts = []
    for path in glob.glob(os.path.join(os.path.dirname(HERE), "results", "counterexamples", "*k4.txt")):
        texts.append(open(path).read())
    for _ in range(count):
        # ("thin" and "clusters" cannot produce large hulls quickly)
        kind = rng.choice(["twochain", "valtr", "triangle", "ellipse", "disk"])
        n = rng.randint(8, nmax)
        pts = gen(rng, n, kind) if kind in ("ellipse", "thin", "disk") else family(rng, n, kind)
        texts.append(f"{len(pts)} 2 4\n" + "".join(f"{x} {y}\n" for x, y in pts))
    bad = 0
    for t in texts:
        a, _ = solve(t, "tight4")
        b, _ = solve(t, "capped")
        if a != b:
            bad += 1
            print("MISMATCH", a, b)
            print(t)
    print(f"{len(texts)} instances, {bad} mismatches")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
