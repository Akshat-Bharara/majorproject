#!/usr/bin/env python3
"""Search for C1 counterexamples with larger k by lifting a known one.

Takes a C1 counterexample (P, k) and emits random instances P + {extra points}
with parameter k + 1 (extra points sampled around P, kept only if every point
stays a strictly convex hull vertex). Pipe into `build/v5check --c1`.

usage: lift_search.py FILE SEED COUNT [EXTRA]
"""

import math
import random
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from structure import hull   # noqa: E402


def main():
    path, seed, count = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    extra_max = int(sys.argv[4]) if len(sys.argv) > 4 else 3
    tok = open(path).read().split()
    n, k = int(tok[0]), int(tok[2])
    P = [(int(tok[3 + 2 * i]), int(tok[4 + 2 * i])) for i in range(n)]
    cx = sum(x for x, _ in P) / n
    cy = sum(y for _, y in P) / n
    R = max(math.dist((cx, cy), p) for p in P)
    rng = random.Random(seed)
    out = sys.stdout
    made = 0
    while made < count:
        Q = list(P)
        for _ in range(rng.randint(1, extra_max)):
            a = rng.uniform(0, 2 * math.pi)
            rad = R * rng.uniform(0.9, 2.5)
            Q.append((round(cx + rad * math.cos(a)), round(cy + rad * math.sin(a))))
        H = hull(Q)
        if len(H) != len(Q):
            continue
        made += 1
        out.write(f"{len(H)} 2 {k + 1}\n" + "".join(f"{x} {y}\n" for x, y in H))


if __name__ == "__main__":
    main()
