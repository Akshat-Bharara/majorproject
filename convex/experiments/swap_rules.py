#!/usr/bin/env python3
"""Which s makes the swap S - succ(s) + a(s) succeed? Test selection rules.

For every r-far k-set S of small random instances, and every rule that picks
one s in S, count how often the swap at that s fails (or s is already tight,
which counts as success).
"""

import collections
import itertools
import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from structure import Inst, gen          # noqa: E402
from greedy_stress import family         # noqa: E402


def angle(P, a, b, c):  # angle at b
    ax, ay = P[a][0] - P[b][0], P[a][1] - P[b][1]
    cx, cy = P[c][0] - P[b][0], P[c][1] - P[b][1]
    return math.atan2(abs(ax * cy - ay * cx), ax * cx + ay * cy)


def main():
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
    trials = int(sys.argv[2]) if len(sys.argv) > 2 else 150
    fails = collections.Counter()
    sets = 0
    for _ in range(trials):
        kind = rng.choice(["twochain", "valtr", "triangle", "clusters", "ellipse", "thin"])
        n = rng.randint(7, 12)
        P = gen(rng, n, kind) if kind in ("ellipse", "thin") else family(rng, n, kind)
        I = Inst(P); n = I.n; P = I.P; D = I.D
        k = rng.randint(3, min(6, n - 1))
        for r2 in sorted({D[a][b] for a, b in itertools.combinations(range(n), 2)}):
            A = [I.a(i, r2, 1) for i in range(n)]
            for S in I.far_sets(range(n), k, r2):
                S = list(S); sets += 1
                succ = {S[i]: S[(i + 1) % k] for i in range(k)}
                pred = {S[(i + 1) % k]: S[i] for i in range(k)}

                def ok(s):
                    t = succ[s]
                    if A[s] == t:
                        return True
                    T = [x for x in S if x != t] + [A[s]]
                    return all(D[a][b] >= r2 for a, b in itertools.combinations(T, 2))

                rules = {
                    "shortest_succ_edge": min(S, key=lambda s: D[s][succ[s]]),
                    "longest_succ_edge": max(S, key=lambda s: D[s][succ[s]]),
                    "closest_pair_first": min(S, key=lambda s: min(D[s][x] for x in S if x != s)),
                    "min_interior_angle_at_s": min(S, key=lambda s: angle(P, pred[s], s, succ[s])),
                    "max_interior_angle_at_s": max(S, key=lambda s: angle(P, pred[s], s, succ[s])),
                    "min_angle_at_succ": min(S, key=lambda s: angle(P, s, succ[s], succ[succ[s]])),
                    "max_angle_at_succ": max(S, key=lambda s: angle(P, s, succ[s], succ[succ[s]])),
                    "max_gap_to_a": max(S, key=lambda s: (succ[s] - A[s]) % n),
                    "min_D_a_succ": min(S, key=lambda s: D[A[s]][succ[s]]),
                }
                for name, s in rules.items():
                    if not ok(s):
                        fails[name] += 1
    print("sets", sets)
    for name in sorted(fails, key=fails.get):
        print(f"{name:26s} fails {fails[name]}")


if __name__ == "__main__":
    main()
