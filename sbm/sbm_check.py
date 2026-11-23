#!/usr/bin/env python3
"""Exact re-implementation of the decision procedure of Singireddy, Basappa,
Mitchell (DAM 386 (2026), Sec. 3, Claim 2 / Lemma 3): two-way search tree.

Place d1 at a vertex v. Keep a clockwise and a counterclockwise frontier
(the last disk placed on each side, both = v initially). The next disk goes
to the FIRST vertex beyond one frontier (walking away from v, not passing the
other frontier) whose distance to EVERY placed centre is >= 2r. Branch on the
side. YES iff some root-leaf path of some root reaches k disks.
Distances are squared; D2 = (2r)^2.
"""
import sys
from itertools import combinations


def d2(P, i, j):
    return (P[i][0] - P[j][0]) ** 2 + (P[i][1] - P[j][1]) ** 2


def sbm_decide(P, k, D2, trace=None):
    n = len(P)

    def rec(S, lo, hi):          # free vertices: offsets lo+1 .. hi-1 (ccw from root)
        if len(S) == k:
            return S
        q = S[0]
        for ccw in (True, False):
            rng = range(lo + 1, hi) if ccw else range(hi - 1, lo, -1)
            for t in rng:
                v = (q + t) % n
                if all(d2(P, u, v) >= D2 for u in S):
                    res = rec(S + [v], t if ccw else lo, hi if ccw else t)
                    if res:
                        return res
                    break      # only the FIRST feasible vertex on this side
        return None

    for v in range(n):
        res = rec([v], 0, n)
        if res:
            return sorted(res)
    return None


def brute(P, k):
    best, arg = -1, None
    for S in combinations(range(len(P)), k):
        m = min(d2(P, a, b) for a, b in combinations(S, 2))
        if m > best:
            best, arg = m, S
    return best, arg


def sbm_optimize(P, k):
    """Algorithm 1 as described in the text: repeatedly test the (upper)
    median of the remaining candidates; YES keeps the elements >= median,
    NO keeps those < median."""
    X = sorted({d2(P, i, j) for i in range(len(P)) for j in range(i + 1, len(P))})
    while len(X) >= 2:
        m = X[len(X) // 2]
        X = [e for e in X if e >= m] if sbm_decide(P, k, m) else [e for e in X if e < m]
    return X[0]


if __name__ == "__main__":
    Pstar = [(16, 48), (-83, 26), (-114, 18), (-30, -39), (-16, -48), (83, -26), (114, -18), (30, 39)]
    k = 4
    opt, S = brute(Pstar, k)
    print("P* : exhaustive optimum (2 r_max)^2 =", opt, "set", S)
    print("     SBM decision at the optimum      :", sbm_decide(Pstar, k, opt))
    print("     SBM Algorithm 1 returns (2r)^2   =", sbm_optimize(Pstar, k))
