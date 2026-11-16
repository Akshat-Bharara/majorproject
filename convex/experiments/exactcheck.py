#!/usr/bin/env python3
"""Independent exhaustive checks in exact integer arithmetic (pure Python).

Shared helpers for check_counterexamples.py and verify_python.py; nothing
here uses the C++ code or the dynamic program.
"""
from itertools import combinations


def cross(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def strictly_convex_ccw(P):
    n = len(P)
    if n < 3:
        return False
    # every point strictly to the left of every directed hull edge
    return all(cross(P[i], P[(i + 1) % n], P[j]) > 0
               for i in range(n) for j in range(n) if j != i and j != (i + 1) % n)


def d2(P, i, j):
    return (P[i][0] - P[j][0]) ** 2 + (P[i][1] - P[j][1]) ** 2


def far_sets(P, k, r2, within=None):
    """All k-subsets (of `within`, default all points) with pairwise squared distance >= r2."""
    idx = range(len(P)) if within is None else within
    return [S for S in combinations(idx, k)
            if all(d2(P, a, b) >= r2 for a, b in combinations(S, 2))]


def optimum(P, k):
    return max(min(d2(P, a, b) for a, b in combinations(S, 2)) for S in combinations(range(len(P)), k))


def first_far(P, r2, q, direction):
    """a(q) (direction=+1) or b(q) (direction=-1): first point at distance >= r; None if none."""
    n = len(P)
    for t in range(1, n):
        v = (q + direction * t) % n
        if d2(P, q, v) >= r2:
            return v
    return None


def arc(P, r2, i):
    """A(i): the points from a(i) counterclockwise to b(i)."""
    n = len(P)
    a, b = first_far(P, r2, i, +1), first_far(P, r2, i, -1)
    if a is None:
        return []
    out, v = [], a
    while True:
        out.append(v)
        if v == b:
            return out
        v = (v + 1) % n


class Seeded:
    """Seeded greedies; offsets are counted counterclockwise from the seed q."""

    def __init__(self, P, k, r2):
        self.P, self.k, self.r2, self.n = P, k, r2, len(P)

    def ok(self, v, S):
        return all(d2(self.P, u, v) >= self.r2 for u in S)

    def seed(self, q, seed3):
        n = self.n
        a, b = first_far(self.P, self.r2, q, +1), first_far(self.P, self.r2, q, -1)
        if a is None:
            return None
        ta, tb = (a - q) % n, (q - b) % n
        if not seed3:
            return [q, a], ta, n
        if n - tb <= ta or d2(self.P, a, b) < self.r2:
            return None
        return [q, a, b], ta, n - tb

    def step(self, q, ccw, S, lo, hi):
        rng = range(lo + 1, hi) if ccw else range(hi - 1, lo, -1)
        for t in rng:
            v = (q + t) % self.n
            if self.ok(v, S):
                return (S + [v], t, hi) if ccw else (S + [v], lo, t)
        return None

    def alternating(self, q, seed3, first_ccw):
        st = self.seed(q, seed3)
        if st is None:
            return False
        side = first_ccw
        while len(st[0]) < self.k:
            st = self.step(q, side, *st)
            if st is None:
                return False
            side = not side
        return True

    def tree(self, q, seed3):
        st = self.seed(q, seed3)
        if st is None:
            return False

        def rec(S, lo, hi):
            if len(S) >= self.k:
                return True
            for ccw in (True, False):
                nx = self.step(q, ccw, S, lo, hi)
                if nx is not None and rec(*nx):
                    return True
            return False

        return rec(*st)


def properties(P, k, r2):
    """Truth values of the properties of Section 6 at (P, k, r)."""
    n = len(P)
    sets = far_sets(P, k, r2)
    feasible = bool(sets)
    a = [first_far(P, r2, q, +1) for q in range(n)]
    b = [first_far(P, r2, q, -1) for q in range(n)]
    out = {"feasible": feasible, "num_far_sets": len(sets)}
    out["G_rhs"] = any(far_sets(P, k - 1, r2, within=arc(P, r2, i)) for i in range(n))
    out["G"] = feasible == out["G_rhs"]
    if feasible:
        out["C1"] = any(q in S and a[q] in S for S in sets for q in range(n))
        out["C1p"] = any(q in S and a[q] in S and b[q] in S and a[q] != b[q] for S in sets for q in range(n))
        sd = Seeded(P, k, r2)
        out["V1"] = any(sd.alternating(q, True, s) for q in range(n) for s in (True, False))
        out["V3"] = any(sd.tree(q, True) for q in range(n))
        out["V5"] = any(sd.alternating(q, False, s) for q in range(n) for s in (True, False))
        out["V4"] = any(sd.tree(q, False) for q in range(n))
    return out
