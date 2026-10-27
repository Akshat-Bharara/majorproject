#!/usr/bin/env python3
"""Tight-pair conjecture C1 and greedy variants seeded by a tight pair.

C1  every feasible (P, k, r) has a feasible k-set containing q and a(q) for
    some q (a consecutive pair of the solution that is "tight": a(q) is the
    first point ccw from q at distance >= r).
V4  seed {q, a(q)}; extend greedily (first far point beyond a frontier) with
    every left/right pattern of length k-2.
V5  seed {q, a(q)}; alternate the two frontiers (both starting sides).
"""

import itertools
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from structure import Inst          # noqa: E402
from greedy_conj import extend      # noqa: E402
from greedy_stress import family, rstar  # noqa: E402


def seeded(I, q, k, r2, pattern):
    a = I.a(q, r2, 1)
    if a is None:
        return False
    chosen, x, y = [q, a], a, q   # ccw frontier a, cw frontier q
    t = 0
    while len(chosen) < k:
        cl = extend(I, chosen, x, y, 1, r2)
        cr = extend(I, chosen, y, x, -1, r2)
        side = pattern(t, cl, cr)
        if side == "L" and cl is not None:
            chosen.append(cl); x = cl
        elif side == "R" and cr is not None:
            chosen.append(cr); y = cr
        else:
            return False
        t += 1
    return True


def v4(I, q, k, r2):
    return any(seeded(I, q, k, r2, lambda t, cl, cr, b=bits: b[t])
               for bits in itertools.product("LR", repeat=k - 2))


def v5(I, q, k, r2):
    return (seeded(I, q, k, r2, lambda t, cl, cr: "LR"[t % 2]) or
            seeded(I, q, k, r2, lambda t, cl, cr: "RL"[t % 2]))


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    trials = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    rng = random.Random(seed)
    cnt = {"tests": 0, "C1": 0, "V4": 0, "V5": 0}
    for trial in range(trials):
        kind = rng.choice(["twochain", "valtr", "triangle", "clusters"])
        I = Inst(family(rng, rng.randint(8, 14), kind))
        n = I.n
        k = rng.randint(3, min(7, n - 1))
        rs = rstar(I, k)
        A = [I.a(i, 0, 1) for i in range(n)]
        for r2 in sorted({I.D[a][b] for a, b in itertools.combinations(range(n), 2)}):
            if r2 > rs:
                break
            cnt["tests"] += 1
            A = [I.a(i, r2, 1) for i in range(n)]
            sets = list(I.far_sets(range(n), k, r2))
            if not any(A[q] in s for s in sets for q in s):
                cnt["C1"] += 1
                print("C1 FAIL", kind, n, k, r2, I.P); sys.stdout.flush()
            if not any(v4(I, q, k, r2) for q in range(n)):
                cnt["V4"] += 1
                if cnt["V4"] <= 2:
                    print("V4 FAIL", kind, n, k, r2, I.P); sys.stdout.flush()
            if not any(v5(I, q, k, r2) for q in range(n)):
                cnt["V5"] += 1
    print(f"seed {seed}: {cnt}")


if __name__ == "__main__":
    main()
