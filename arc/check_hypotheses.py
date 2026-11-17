#!/usr/bin/env python3
"""Independent check of the counterexamples to the weakened hypotheses.

    python3 -I check_hypotheses.py data/hypothesis_counterexamples.json

H1-k6 : an r-far k-set S exists with Far_r(s) arc-shaped for every s in S,
        yet no greedy chain is r-far.
H2u-k4: the arc property holds at r = r*, yet binary search with the
        UNCHECKED chain test (no wrap + closing distance) does not return r*.
H3-k4 : (no hypothesis) no greedy chain is r*-far.
Pstar : the arc property fails at r^2 = 10240 and no greedy chain is r-far.
"""
import json
import sys
from itertools import combinations


def d2(P, i, j):
    return (P[i][0] - P[j][0]) ** 2 + (P[i][1] - P[j][1]) ** 2


def convex(P):
    n = len(P)
    cr = lambda a, b, c: (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    return all(cr(P[i], P[(i + 1) % n], P[j]) > 0 for i in range(n) for j in range(n) if j not in (i, (i + 1) % n))


def arcs(P, q, r2):
    n = len(P)
    far = [d2(P, q, (q + t) % n) >= r2 for t in range(1, n)]
    return sum(1 for t in range(n - 1) if far[t] and (t == 0 or not far[t - 1]))


def chain(P, q, k, r2):
    n = len(P)
    out, x, total = [q], q, 0
    for _ in range(k - 1):
        step = next((t for t in range(1, n) if d2(P, x, (x + t) % n) >= r2), None)
        if step is None or total + step >= n:
            return None
        total += step
        x = (x + step) % n
        out.append(x)
    return out


def greedy(P, k, r2):
    return any(c and all(d2(P, a, b) >= r2 for a, b in combinations(c, 2)) for c in (chain(P, q, k, r2) for q in range(len(P))))


def unchecked(P, k, r2):
    return any(c and d2(P, c[0], c[-1]) >= r2 for c in (chain(P, q, k, r2) for q in range(len(P))))


def far_sets(P, k, r2):
    return [S for S in combinations(range(len(P)), k) if all(d2(P, a, b) >= r2 for a, b in combinations(S, 2))]


def main():
    bad = 0
    for inst in json.load(open(sys.argv[1]))["instances"]:
        P, k, r2, name = [tuple(p) for p in inst["points"]], inst["k"], inst["r2"], inst["name"]
        n = len(P)
        cand = sorted({d2(P, i, j) for i in range(n) for j in range(i + 1, n)})
        opt = max(r for r in cand if far_sets(P, k, r))
        ok = convex(P)
        if name.startswith("H1"):
            good = [S for S in far_sets(P, k, r2) if all(arcs(P, s, r2) <= 1 for s in S)]
            ok = ok and bool(good) and not greedy(P, k, r2)
            info = f"{len(good)} r-far k-set(s) with arc-shaped Far sets, greedy fails"
        elif name.startswith("H2u"):
            lo, hi = 0, len(cand) - 1
            while lo < hi:
                mid = (lo + hi + 1) // 2
                if unchecked(P, k, cand[mid]):
                    lo = mid
                else:
                    hi = mid - 1
            ok = ok and r2 == opt and all(arcs(P, q, opt) <= 1 for q in range(n)) and cand[lo] != opt
            info = f"r*^2 = {opt}, arc property at r*, unchecked binary search returns {cand[lo]}"
        elif name.startswith("H3"):
            ok = ok and r2 == opt and not greedy(P, k, opt)
            info = f"r*^2 = {opt}, no greedy chain is r*-far"
        else:
            beta = sum(max(0, arcs(P, q, r2) - 1) for q in range(n))
            ok = ok and r2 == opt and beta > 0 and not greedy(P, k, r2)
            info = f"r*^2 = {opt}, beta = {beta}, no greedy chain is r-far"
        bad += not ok
        print(("ok   " if ok else "FAIL ") + f"{name:7s} n={n} k={k}: {info}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
