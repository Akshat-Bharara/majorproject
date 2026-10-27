#!/usr/bin/env python3
"""Shrink counterexamples to the greedy conjectures by deleting points.

usage: minimize.py {V1|V5} FILE      (FILE: "n 2 k" + points)

An instance is a counterexample for variant V if for some r <= r* (feasible)
no seed makes the greedy succeed. r* comes from the exact C++ solver and is
re-checked here by exhaustive search on the final instance. Prints the
minimal instance and the witnessing r^2.
"""

import itertools
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from structure import Inst, hull       # noqa: E402
from greedy_conj import v1              # noqa: E402
from tight_pair import v5               # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(os.path.dirname(HERE), "build", "convex_dispersion")


def rstar(P, k):
    txt = f"{len(P)} 2 {k}\n" + "".join(f"{x} {y}\n" for x, y in P)
    out = subprocess.run([EXE], input=txt, capture_output=True, text=True, check=True).stdout
    return int(out.split()[1])


def witness(P, k, var):
    P = hull(P)
    if len(P) <= k:
        return None
    I = Inst(P)
    rs = rstar(P, k)
    f = v1 if var == "V1" else v5
    for r2 in sorted({I.D[a][b] for a, b in itertools.combinations(range(I.n), 2)}):
        if r2 > rs:
            return None
        if not any(f(I, q, k, r2) for q in range(I.n)):
            return r2
    return None


def main():
    var, path = sys.argv[1], sys.argv[2]
    tok = open(path).read().split()
    n, k = int(tok[0]), int(tok[2])
    P = [(int(tok[3 + 2 * i]), int(tok[4 + 2 * i])) for i in range(n)]
    assert witness(P, k, var) is not None, "not a counterexample"
    changed = True
    while changed:
        changed = False
        for i in range(len(P)):
            Q = P[:i] + P[i + 1:]
            if witness(Q, k, var) is not None:
                P, changed = hull(Q), True
                break
    r2 = witness(P, k, var)
    I = Inst(P)
    # independent exhaustive confirmation that k points pairwise >= r exist
    sets = list(I.far_sets(range(I.n), k, r2))
    print(f"{var} counterexample: n={len(P)} k={k} r^2={r2} feasible_sets={len(sets)}")
    print(f"{len(P)} 2 {k}")
    for x, y in P:
        print(x, y)
    print("example far set:", [P[i] for i in sets[0]])


if __name__ == "__main__":
    main()
