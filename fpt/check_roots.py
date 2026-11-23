#!/usr/bin/env python3
"""Independent check (exact integers) of the instances on which natural root
rules for the two-sided search of Singireddy-Basappa-Mitchell fail.

    python3 -I check_roots.py data/root_rule_instances.json

For each instance: the points are in strictly convex ccw position; the
two-sided search (SBM, Sec. 3) succeeds from SOME root (Claim 2 holds), but
  closest-pair-*      : it fails from an endpoint of a closest pair of T;
  smallest-angle-*    : the witness-preserving search fails from the vertex of
                        T with the smallest interior angle;
  global-smallest-*   : ... from the vertex with the smallest angle over all
                        r-far k-sets.
"""
import json, math, sys
from itertools import combinations
from fractions import Fraction


def d2(P, a, b): return (P[a][0]-P[b][0])**2 + (P[a][1]-P[b][1])**2


def convex(P):
    n = len(P); cr = lambda a, b, c: (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])
    return all(cr(P[i], P[(i+1) % n], P[j]) > 0 for i in range(n) for j in range(n) if j not in (i, (i+1) % n))


def sbm(P, k, r2, root):
    n = len(P)
    def rec(C, lo, hi):
        if len(C) == k: return True
        for ccw in (True, False):
            rng = range(lo+1, hi) if ccw else range(hi-1, lo, -1)
            for t in rng:
                v = (root+t) % n
                if all(d2(P, u, v) >= r2 for u in C):
                    if rec(C+[v], t if ccw else lo, hi if ccw else t): return True
                    break
        return False
    return rec([root], 0, n)


def witness_search(P, k, r2, T, root):
    """two-sided search that keeps a witness: only steps x for which
    (W minus the next remaining point on that side) plus x is r-far."""
    n = len(P); off = lambda v: (v-root) % n
    rem0 = sorted([v for v in T if v != root], key=off)
    def rec(C, rem, lo, hi):
        if len(C) == k: return True
        for ccw in (True, False):
            rng = range(lo+1, hi) if ccw else range(hi-1, lo, -1)
            for t in rng:
                v = (root+t) % n
                if all(d2(P, u, v) >= r2 for u in C):
                    r2m = rem[1:] if ccw else rem[:-1]
                    if all(d2(P, s, v) >= r2 for s in r2m) and rec(C+[v], r2m, t if ccw else lo, hi if ccw else t):
                        return True
                    break
        return False
    return rec([root], rem0, 0, n)


def cos_angle(P, a, p, b):  # exact comparison key: cos with sign, as a Fraction of squares
    ax, ay = P[a][0]-P[p][0], P[a][1]-P[p][1]; bx, by = P[b][0]-P[p][0], P[b][1]-P[p][1]
    dot = ax*bx + ay*by
    return (1 if dot >= 0 else -1) * Fraction(dot*dot, (ax*ax+ay*ay)*(bx*bx+by*by))  # monotone in cos


def angles(P, T):
    k = len(T); return [cos_angle(P, T[(i-1) % k], T[i], T[(i+1) % k]) for i in range(k)]  # larger = smaller angle


def main():
    bad = 0
    for inst in json.load(open(sys.argv[1]))["instances"]:
        P = [tuple(p) for p in inst["points"]]; k, r2, name = inst["k"], inst["r2"], inst["name"]; n = len(P)
        sets = [S for S in combinations(range(n), k) if all(d2(P, a, b) >= r2 for a, b in combinations(S, 2))]
        claim2 = any(sbm(P, k, r2, v) for v in range(n))
        ok = convex(P) and claim2 and bool(sets)
        if name.startswith("closest-pair"):
            T = inst["T"]; ok = ok and tuple(T) in sets
            dmin = min(d2(P, a, b) for a, b in combinations(T, 2))
            ends = sorted({v for a, b in combinations(T, 2) if d2(P, a, b) == dmin for v in (a, b)})
            fails = [v for v in ends if not sbm(P, k, r2, v)]
            ok = ok and bool(fails)
            info = f"closest-pair endpoints {ends}; two-sided search fails from {fails}"
        elif name.startswith("smallest-angle"):
            T = inst["T"]; ok = ok and tuple(T) in sets
            A = angles(P, T); u = T[max(range(k), key=lambda i: A[i])]
            f = not witness_search(P, k, r2, T, u)
            ok = ok and f
            info = f"smallest-angle vertex {u}: witness-preserving search fails = {f}; succeeds from some vertex of T = {any(witness_search(P, k, r2, T, v) for v in T)}"
        elif name.startswith("stuck-vertex"):
            T = inst["T"]; ok = ok and tuple(T) in sets
            good = [v for v in T if witness_search(P, k, r2, T, v)]
            ok = ok and T[0] not in good and bool(good)
            info = f"vertex {T[0]} of T is not good (stuck after 0 steps, 3 points remaining); good vertices {good}"
        else:
            best = None
            for T in sets:
                A = angles(P, list(T))
                for i in range(k):
                    if best is None or A[i] > best[0]: best = (A[i], list(T), T[i])
            _, T, u = best
            f = not witness_search(P, k, r2, T, u)
            ok = ok and f
            info = f"globally smallest angle at {u} in T={T}: witness-preserving search fails = {f}"
        bad += not ok
        print(("ok   " if ok else "FAIL ") + f"{name}: n={n}, k={k}, {len(sets)} r-far k-sets, Claim 2 holds: {claim2}; {info}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
