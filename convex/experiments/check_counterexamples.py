#!/usr/bin/env python3
"""Re-verify every explicit instance of Section 6 by exhaustive search.

    python3 experiments/check_counterexamples.py [data/counterexamples.json]

Uses only exact Python integers and brute force over all k-subsets; it shares
no code with the C++ solvers.  Also checks the G_k family of Proposition
"G_k fails for every k >= 4" for k = 4..12.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exactcheck import far_sets, properties, strictly_convex_ccw, arc, optimum  # noqa: E402


def lens_family(k):
    """The family of Proposition 6.2: u alone on the left parabola, k-1 points on
    the right one (levels 0, +-1, +-2, ...); u and v = level 0 are the only close pair."""
    levels = [0] + [j for t in range(1, k) for j in (t, -t)][: k - 2]
    H = 3 * (k + 1)
    right = sorted(((H * H - 9 * j * j, 3 * j * H * H) for j in levels), key=lambda p: p[1])
    return [(-H * H, 0)] + right, (3 * H * H) ** 2


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "data", "counterexamples.json")
    bad = 0
    for inst in json.load(open(path))["instances"]:
        P = [tuple(p) for p in inst["points"]]
        k, r2 = inst["k"], inst["r2"]
        ok = strictly_convex_ccw(P)
        got = properties(P, k, r2)
        if "optimum_r2" in inst["expect"]:
            got["optimum_r2"] = optimum(P, k)
        for key, want in inst["expect"].items():
            ok = ok and got.get(key) == want
        bad += not ok
        shown = {key: got.get(key) for key in inst["expect"]}
        print(f"{'ok  ' if ok else 'FAIL'} {inst['name']:16s} n={len(P):2d} k={k} r2={r2}: {shown}")
    for k in range(4, 13):
        P, r2 = lens_family(k)
        ok = (len(P) == k and strictly_convex_ccw(P) and not far_sets(P, k, r2)
              and any(far_sets(P, k - 1, r2, within=arc(P, r2, i)) for i in range(k)))
        bad += not ok
        print(f"{'ok  ' if ok else 'FAIL'} lens family k={k}: no r-far k-set, but some arc A(i) holds an r-far (k-1)-set")
    print("all verified" if not bad else f"{bad} FAILED")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
