#!/usr/bin/env python3
"""Lemma B' (isolated cap): numerical test and adversarial search.

Q = (q_0..q_{k-1}) convex, ccw, pairwise distances >= r (we take r = min
pairwise distance, the hardest case). For the edge (s, t) = (q_i, q_{i+1}) the
cap region is
  R(s,t) = { w : w strictly right of s->t (outside Q),
                 orient(q_{i-1}, s, w) > 0, orient(w, t, q_{i+2}) > 0
                 (Q with w inserted between s and t stays convex),
                 |s w| >= r }.
Claim: some edge has dist(u, R(s,t)) >= r for every vertex u != s, t.
dist(u, R) is computed by dense sampling of the boundary of R (the closest
point of R to a vertex u lies on that boundary), refined around the minimum.
"""

import math
import random
import sys

sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from isolated_edge import hull, random_poly  # noqa: E402


def orient(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def line_inter(p1, p2, p3, p4):
    d = (p2[0] - p1[0]) * (p4[1] - p3[1]) - (p2[1] - p1[1]) * (p4[0] - p3[0])
    if abs(d) < 1e-15:
        return None
    t = ((p3[0] - p1[0]) * (p4[1] - p3[1]) - (p3[1] - p1[1]) * (p4[0] - p3[0])) / d
    return (p1[0] + t * (p2[0] - p1[0]), p1[1] + t * (p2[1] - p1[1]))


def in_R(w, Q, i, r):
    k = len(Q)
    s, t = Q[i], Q[(i + 1) % k]
    sp, tn = Q[(i - 1) % k], Q[(i + 2) % k]
    return (orient(s, t, w) < 0 and orient(sp, s, w) > 0 and orient(w, t, tn) > 0
            and math.dist(s, w) >= r)


def cap_boundary_samples(Q, i, r, m=400):
    """Points on the boundary of R(s,t) (clipped to a big box)."""
    k = len(Q)
    s, t = Q[i], Q[(i + 1) % k]
    sp, tn = Q[(i - 1) % k], Q[(i + 2) % k]
    span = 10 * max(math.dist(Q[a], Q[b]) for a in range(k) for b in range(k))
    # apex of the wedge (may not exist / be far)
    X = line_inter(sp, s, t, tn)
    # direction of the two wedge rays beyond s and t
    ds = (s[0] - sp[0], s[1] - sp[1]); ls = math.hypot(*ds); ds = (ds[0] / ls, ds[1] / ls)
    dt = (t[0] - tn[0], t[1] - tn[1]); lt = math.hypot(*dt); dt = (dt[0] / lt, dt[1] / lt)
    pts = []
    def seg(a, b):
        for j in range(m + 1):
            u = j / m
            pts.append((a[0] + u * (b[0] - a[0]), a[1] + u * (b[1] - a[1])))
    if X is not None and orient(s, t, X) < 0:     # bounded triangle s, X, t
        seg(s, X); seg(X, t); seg(s, t)
    else:                                         # unbounded: rays from s and t
        seg(s, (s[0] + span * ds[0], s[1] + span * ds[1]))
        seg(t, (t[0] + span * dt[0], t[1] + span * dt[1]))
        seg(s, t)
    for j in range(4 * m + 1):                    # circle |sw| = r
        a = 2 * math.pi * j / (4 * m)
        pts.append((s[0] + r * math.cos(a), s[1] + r * math.sin(a)))
    eps = 1e-12
    return [p for p in pts if in_R_closed(p, Q, i, r, eps)]


def in_R_closed(w, Q, i, r, eps):
    k = len(Q)
    s, t = Q[i], Q[(i + 1) % k]
    sp, tn = Q[(i - 1) % k], Q[(i + 2) % k]
    sc = max(1.0, max(abs(c) for c in s))
    return (orient(s, t, w) <= eps * sc and orient(sp, s, w) >= -eps * sc
            and orient(w, t, tn) >= -eps * sc and math.dist(s, w) >= r * (1 - 1e-12))


def cap_iso(Q):
    k = len(Q)
    r = min(math.dist(Q[a], Q[b]) for a in range(k) for b in range(a + 1, k))
    best = 0.0
    for i in range(k):
        B = cap_boundary_samples(Q, i, r)
        if not B:
            return float("inf")   # empty cap: this edge trivially works
        m = min(math.dist(Q[j], w) for j in range(k) if j not in (i, (i + 1) % k) for w in B)
        best = max(best, m / r)
    return best


def main():
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
    restarts = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    fixed_k = int(sys.argv[3]) if len(sys.argv) > 3 else None
    worst = (9.0, None)
    for _ in range(restarts):
        Q = random_poly(rng, fixed_k or rng.randint(3, 9))
        while fixed_k and len(Q) != fixed_k:
            Q = random_poly(rng, fixed_k)
        cur = cap_iso(Q)
        step = 0.05
        for it in range(400):
            i = rng.randrange(len(Q))
            R = list(Q)
            R[i] = (R[i][0] + rng.gauss(0, step), R[i][1] + rng.gauss(0, step))
            R = hull(R)
            if len(R) != len(Q):
                continue
            v = cap_iso(R)
            if v <= cur:
                Q, cur = R, v
            if it % 100 == 99:
                step *= 0.5
        if cur < worst[0]:
            worst = (cur, Q)
        if cur < 1 - 1e-6:
            print("COUNTEREXAMPLE cap_iso =", cur, "k =", len(Q), Q)
            sys.stdout.flush()
            if fixed_k:
                break
    print("smallest cap isolation found:", worst[0], "k =", len(worst[1]))
    print(worst[1])


if __name__ == "__main__":
    main()
