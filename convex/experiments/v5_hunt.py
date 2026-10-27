#!/usr/bin/env python3
"""Generate many random convex instances and stream them to v5check.

usage: v5_hunt.py SEED COUNT NMAX KMAX | build/v5check
"""

import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from greedy_stress import family  # noqa: E402
from structure import gen         # noqa: E402

seed, count, nmax, kmax = map(int, sys.argv[1:5])
rng = random.Random(seed)
out = sys.stdout
for _ in range(count):
    kind = rng.choice(["twochain", "twochain", "valtr", "triangle", "clusters",
                       "ellipse", "thin", "disk"])
    n = rng.randint(8, nmax)
    pts = gen(rng, n, kind) if kind in ("ellipse", "thin", "disk") else family(rng, n, kind)
    k = rng.randint(3, min(kmax, len(pts) - 1))
    out.write(f"{len(pts)} 2 {k}\n")
    for x, y in pts:
        out.write(f"{x} {y}\n")
