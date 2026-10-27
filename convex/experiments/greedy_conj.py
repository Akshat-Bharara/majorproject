#!/usr/bin/env python3
"""Two-sided greedy conjectures for k-dispersion in convex position.

Start from a seed q and the set {q, a(q), b(q)} (conjecture C1'). Keep a ccw
frontier x and a cw frontier y; extending a side means adding the first point
beyond the frontier (towards the other frontier) that is r-far from every
chosen point. Variants:
  V1  alternate ccw, cw, ccw, ...
  V2  extend the side whose candidate is closer in index (ties: ccw)
  V3  try every left/right pattern (2^(k-3) per seed)
A variant "succeeds" on (P, k, r) if some seed q yields k points. Compared with
the exhaustive answer; reports failures (variant says no, brute says yes).
"""

import collections
import itertools
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from structure import Inst, gen  # noqa: E402


def extend(I, chosen, frontier, other, step, r2):
    """First point strictly between frontier and other (moving by step) far from all chosen."""
    n = I.n
    j = (frontier + step) % n
    while j != other:
        if all(I.D[j][c] >= r2 for c in chosen):
            return j
        j = (j + step) % n
    return None


def run_pattern(I, q, k, r2, pattern):
    """pattern: function(step_index, cand_ccw, cand_cw) -> 'L' or 'R' or None."""
    a, b = I.a(q, r2, 1), I.a(q, r2, -1)
    if a is None:
        return False
    chosen = [q, a] if a == b else [q, a, b]
    if a != b and I.D[a][b] < r2:
        return False
    if len(chosen) >= k:
        return True
    x, y = a, b  # ccw frontier, cw frontier
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


def v1(I, q, k, r2):
    return run_pattern(I, q, k, r2, lambda t, cl, cr: "L" if t % 2 == 0 else "R") or \
           run_pattern(I, q, k, r2, lambda t, cl, cr: "R" if t % 2 == 0 else "L")


def v2(I, q, k, r2):
    def pat(t, cl, cr):
        if cl is None:
            return "R"
        if cr is None:
            return "L"
        dl = (cl - q) % I.n
        dr = (q - cr) % I.n
        return "L" if dl <= dr else "R"
    return run_pattern(I, q, k, r2, pat)


def v3(I, q, k, r2):
    for bits in itertools.product("LR", repeat=max(0, k - 3)):
        if run_pattern(I, q, k, r2, lambda t, cl, cr, b=bits: b[t] if t < len(b) else None):
            return True
    return False


def main():
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 21)
    trials = int(sys.argv[2]) if len(sys.argv) > 2 else 800
    st = collections.Counter()
    first_fail = {}
    for trial in range(trials):
        kind = rng.choice(["ellipse", "disk", "thin"])
        n = rng.randint(7, 13)
        I = Inst(gen(rng, n, kind))
        k = rng.randint(3, min(7, n - 1))
        for r2 in sorted({I.D[a][b] for a, b in itertools.combinations(range(n), 2)}):
            truth = I.feasible(range(n), k, r2)
            for name, f in (("V1", v1), ("V2", v2), ("V3", v3)):
                got = any(f(I, q, k, r2) for q in range(n))
                st[(name, k, "tests")] += 1
                if got and not truth:
                    st[(name, k, "UNSOUND")] += 1
                if truth and not got:
                    st[(name, k, "fail")] += 1
                    first_fail.setdefault((name, k), (kind, I.P, r2))
    for key in sorted(st):
        print(key, st[key])
    for key, v in sorted(first_fail.items()):
        print("first failure", key, v[0], "n=", len(v[1]), "r2=", v[2])


if __name__ == "__main__":
    main()
