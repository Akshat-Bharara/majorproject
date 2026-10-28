#!/usr/bin/env python3
"""Two-sided cap lemma (B'+-): adversarial search.

For the edge (s,t) of Q there are two caps:
  R_s(s,t): w beyond st, Q+w convex, |s w| >= r   (ccw swap: w plays a(s))
  R_t(s,t): w beyond st, Q+w convex, |t w| >= r   (cw swap:  w plays b(t))
An edge is *safe* if, in one of its caps, every w is at distance >= r from
all other vertices. B'+- claims some edge is safe; by the quadrilateral
inequality this holds whenever an edge of Q joins mutual... nearest
neighbours in either direction. A counterexample is searched for by hill
climbing two_iso(Q) = max_e max(ccw_iso(e), cw_iso(e)) below 1.
"""

import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from isolated_edge import hull, random_poly   # noqa: E402
from cap_lemma import orient, line_inter       # noqa: E402


def cap_samples(Q, i, r, anchor, m=300):
    k = len(Q)
    s, t = Q[i], Q[(i + 1) % k]
    sp, tn = Q[(i - 1) % k], Q[(i + 2) % k]
    c = s if anchor == "s" else t
    span = 10 * max(math.dist(Q[a], Q[b]) for a in range(k) for b in range(k))
    X = line_inter(sp, s, t, tn)
    ds = (s[0] - sp[0], s[1] - sp[1]); l = math.hypot(*ds); ds = (ds[0] / l, ds[1] / l)
    dt = (t[0] - tn[0], t[1] - tn[1]); l = math.hypot(*dt); dt = (dt[0] / l, dt[1] / l)
    pts = []

    def seg(a, b):
        for j in range(m + 1):
            u = j / m
            pts.append((a[0] + u * (b[0] - a[0]), a[1] + u * (b[1] - a[1])))
    if X is not None and orient(s, t, X) < 0:
        seg(s, X); seg(X, t); seg(s, t)
    else:
        seg(s, (s[0] + span * ds[0], s[1] + span * ds[1]))
        seg(t, (t[0] + span * dt[0], t[1] + span * dt[1]))
        seg(s, t)
    for j in range(4 * m + 1):
        a = 2 * math.pi * j / (4 * m)
        pts.append((c[0] + r * math.cos(a), c[1] + r * math.sin(a)))
    sc = max(1.0, max(abs(v) for v in s))
    eps = 1e-12 * sc
    return [w for w in pts
            if orient(s, t, w) <= eps and orient(sp, s, w) >= -eps
            and orient(w, t, tn) >= -eps and math.dist(c, w) >= r * (1 - 1e-12)]


def two_iso(Q):
    k = len(Q)
    r = min(math.dist(Q[a], Q[b]) for a in range(k) for b in range(a + 1, k))
    best = 0.0
    for i in range(k):
        others = [Q[j] for j in range(k) if j not in (i, (i + 1) % k)]
        for anchor in ("s", "t"):
            B = cap_samples(Q, i, r, anchor)
            if not B:
                return float("inf")
            best = max(best, min(math.dist(u, w) for u in others for w in B) / r)
    return best


def main():
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
    restarts = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    fixed_k = int(sys.argv[3]) if len(sys.argv) > 3 else None
    worst = (9.0, None)
    for _ in range(restarts):
        Q = random_poly(rng, fixed_k or rng.randint(4, 8))
        while fixed_k and len(Q) != fixed_k:
            Q = random_poly(rng, fixed_k)
        cur = two_iso(Q)
        step = 0.05
        for it in range(400):
            i = rng.randrange(len(Q))
            R = list(Q)
            R[i] = (R[i][0] + rng.gauss(0, step), R[i][1] + rng.gauss(0, step))
            R = hull(R)
            if len(R) != len(Q):
                continue
            v = two_iso(R)
            if v <= cur:
                Q, cur = R, v
            if it % 100 == 99:
                step *= 0.5
        if cur < worst[0]:
            worst = (cur, Q)
        if cur < 1 - 1e-6:
            print("COUNTEREXAMPLE two_iso =", cur, "k =", len(Q), Q)
            sys.stdout.flush()
            break
    print("smallest two-sided isolation found:", worst[0], "k =", len(worst[1]))
    print(worst[1])


if __name__ == "__main__":
    main()
