#!/usr/bin/env python3
"""Lemma B (isolated edge): in a convex polygon Q whose vertices are pairwise
>= r apart, some edge e has dist(u, e) >= r for every vertex u not on e.

The hardest case is r = min pairwise distance of Q. We test random convex
polygons and run a hill climber that minimises the "best isolation"
  iso(Q) = max over edges e of  min over u not on e of dist(u, e) / r,
looking for iso(Q) < 1 (a counterexample).
"""

import math
import random
import sys


def hull(pts):
    pts = sorted(set(pts))
    if len(pts) < 3:
        return pts
    def cr(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    lo, up = [], []
    for p in pts:
        while len(lo) >= 2 and cr(lo[-2], lo[-1], p) <= 0:
            lo.pop()
        lo.append(p)
    for p in reversed(pts):
        while len(up) >= 2 and cr(up[-2], up[-1], p) <= 0:
            up.pop()
        up.append(p)
    return lo[:-1] + up[:-1]


def seg_dist(p, a, b):
    ax, ay = b[0] - a[0], b[1] - a[1]
    t = ((p[0] - a[0]) * ax + (p[1] - a[1]) * ay) / (ax * ax + ay * ay)
    t = max(0.0, min(1.0, t))
    return math.hypot(p[0] - a[0] - t * ax, p[1] - a[1] - t * ay)


def iso(Q):
    k = len(Q)
    r = min(math.dist(Q[i], Q[j]) for i in range(k) for j in range(i + 1, k))
    best = 0.0
    for i in range(k):
        a, b = Q[i], Q[(i + 1) % k]
        m = min(seg_dist(Q[j], a, b) for j in range(k) if j != i and j != (i + 1) % k)
        best = max(best, m / r)
    return best


def random_poly(rng, k):
    kind = rng.random()
    if kind < 0.4:   # thin two-chain lens
        h = rng.uniform(0.01, 0.5)
        pts = [(x, h * (1 - x * x) * rng.choice([1, -1])) for x in
               (rng.uniform(-1, 1) for _ in range(k))]
    elif kind < 0.7:  # random points on an ellipse
        b = rng.uniform(0.05, 1)
        pts = [(math.cos(t), b * math.sin(t)) for t in (rng.uniform(0, 6.3) for _ in range(k))]
    else:             # hull of random points
        pts = [(rng.random(), rng.random()) for _ in range(6 * k)]
    Q = hull(pts)
    return Q if len(Q) >= 3 else random_poly(rng, k)


def main():
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
    restarts = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    worst = (9, None)
    for _ in range(restarts):
        Q = random_poly(rng, rng.randint(3, 12))
        cur = iso(Q)
        step = 0.05
        for it in range(1500):   # hill climbing on vertex positions
            i = rng.randrange(len(Q))
            R = list(Q)
            R[i] = (R[i][0] + rng.gauss(0, step), R[i][1] + rng.gauss(0, step))
            R = hull(R)
            if len(R) != len(Q):
                continue
            v = iso(R)
            if v <= cur:
                Q, cur = R, v
            if it % 300 == 299:
                step *= 0.5
        if cur < worst[0]:
            worst = (cur, Q)
        if cur < 1 - 1e-9:
            print("COUNTEREXAMPLE iso =", cur, Q)
            break
    print("smallest isolation found:", worst[0], "k =", len(worst[1]))
    print(worst[1])


if __name__ == "__main__":
    main()
