#!/usr/bin/env python3
"""Turn a counterexample to Lemma B' into a candidate counterexample to the
swap lemma W+ (and test C1 on it).

Given Q (k vertices, ccw), r = min pairwise distance of Q. For every edge
(s_i, s_{i+1}) pick a point w_i in the cap R(s_i, s_{i+1}) that is within
distance < r of some other vertex, with margin, by random sampling. Then
P = Q + {w_i}, rounded to integers, must be in strictly convex position; if
so, S = Q is an r-far k-set in which every swap S - s_{i+1} + a(s_i) fails
(a(s_i) = w_i, the only point between s_i and s_{i+1}). Finally C1 and W+
are checked exactly by brute force on P.
"""

import itertools
import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cap_lemma import in_R, cap_boundary_samples  # noqa: E402
from structure import Inst, hull       # noqa: E402

SCALE = int(sys.argv[2]) if len(sys.argv) > 2 else 10 ** 6
MARGIN = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-3


def bad_point(rng, Q, i, r, tries=200000):
    """A point strictly inside the cap R(s,t) within r(1-MARGIN) of a vertex:
    start from the boundary minimiser and nudge it inwards at random."""
    k = len(Q)
    others = [Q[j] for j in range(k) if j not in (i, (i + 1) % k)]
    B = cap_boundary_samples(Q, i, r, m=2000)
    if B:
        d0, w0 = min((min(math.dist(u, w) for u in others), w) for w in B)
        best = None
        for scale in (1e-2, 3e-3, 1e-3, 3e-4, 1e-4):
            for _ in range(4000):
                w = (w0[0] + rng.gauss(0, scale * r), w0[1] + rng.gauss(0, scale * r))
                if not in_R(w, Q, i, r * (1 + MARGIN)):
                    continue
                d = min(math.dist(u, w) for u in others)
                if d < r * (1 - MARGIN) and (best is None or d < best[0]):
                    best = (d, w)
            if best:
                return best
    s, t = Q[i], Q[(i + 1) % k]
    span = max(math.dist(Q[a], Q[b]) for a in range(k) for b in range(k))
    cx, cy = (s[0] + t[0]) / 2, (s[1] + t[1]) / 2
    best = None
    for _ in range(tries):
        w = (cx + rng.uniform(-span, span), cy + rng.uniform(-span, span))
        if not in_R(w, Q, i, r * (1 + MARGIN)):
            continue
        d = min(math.dist(Q[j], w) for j in range(k) if j not in (i, (i + 1) % k))
        if d < r * (1 - MARGIN) and (best is None or d < best[0]):
            best = (d, w)
    return best


def main():
    Q = eval(sys.argv[1])
    rng = random.Random(1)
    k = len(Q)
    r = min(math.dist(Q[a], Q[b]) for a in range(k) for b in range(a + 1, k))
    W = []
    for i in range(k):
        b = bad_point(rng, Q, i, r)
        if b is None:
            print("edge", i, "has no bad cap point with margin")
            return
        W.append(b[1])
        print(f"edge {i}: bad w at distance {b[0] / r:.4f} r")
    pts = []
    for i in range(k):
        pts.append(Q[i]); pts.append(W[i])
    P = [(round(SCALE * x), round(SCALE * y)) for x, y in pts]
    H = hull(P)
    print("convex position:", len(H) == len(P), f"({len(H)} of {len(P)} on hull)")
    if len(H) != len(P):
        return
    I = Inst(H)
    idxQ = [H.index((round(SCALE * x), round(SCALE * y))) for x, y in Q]
    r2 = min(I.D[a][b] for a, b in itertools.combinations(idxQ, 2))
    A = [I.a(i, r2, 1) for i in range(I.n)]
    S = sorted(idxQ)
    swaps = []
    for j in range(k):
        s, t = S[j], S[(j + 1) % k]
        T = [x for x in S if x != t] + [A[s]]
        swaps.append(A[s] == t or all(I.D[a][b] >= r2 for a, b in itertools.combinations(T, 2)))
    print("S = Q is r-far:", all(I.D[a][b] >= r2 for a, b in itertools.combinations(S, 2)))
    print("swap succeeds at each s:", swaps)
    sets = list(I.far_sets(range(I.n), k, r2))
    c1 = any(A[q] in st for st in sets for q in st)
    print(f"r-far {k}-sets: {len(sets)}; C1 holds: {c1}")
    print(f"{len(H)} 2 {k}")
    for x, y in H:
        print(x, y)


if __name__ == "__main__":
    main()
