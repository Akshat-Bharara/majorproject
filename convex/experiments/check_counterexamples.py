#!/usr/bin/env python3
"""Independently verify the counterexamples in results/counterexamples/.

Pure Python, exhaustive, exact integer arithmetic; shares no code with the
C++ solver. For each file it finds every squared radius r^2 at which the
named property fails and writes paper/generated/counterexamples.tex.

  G4  "k points pairwise >= r exist  <=>  some arc A(i) holds k-1 such points"
  V1  alternating greedy seeded with {q, a(q), b(q)} finds k points
  V5  alternating greedy seeded with {q, a(q)} finds k points
"""

import glob
import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from structure import Inst, hull   # noqa: E402
from greedy_conj import v1         # noqa: E402
from tight_pair import v5          # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def load(path):
    tok = open(path).read().split()
    n, k = int(tok[0]), int(tok[2])
    P = [(int(tok[3 + 2 * i]), int(tok[4 + 2 * i])) for i in range(n)]
    H = hull(P)
    assert len(H) == n, "not in strictly convex position"
    return H, k


def failures(prop, I, k):
    out = []
    for r2 in sorted({I.D[a][b] for a, b in itertools.combinations(range(I.n), 2)}):
        feas = I.feasible(range(I.n), k, r2)
        if prop == "G4":
            rhs = any(I.feasible(I.arc(i, r2), k - 1, r2) for i in range(I.n))
            if feas != rhs:
                out.append((r2, feas, rhs))
        else:
            f = v1 if prop == "V1" else v5
            if feas and not any(f(I, q, k, r2) for q in range(I.n)):
                out.append((r2, feas, False))
    return out


def main():
    lines = []
    ok = True
    for path in sorted(glob.glob(os.path.join(ROOT, "results", "counterexamples", "*.txt"))):
        prop = os.path.basename(path).split("_")[0]
        P, k = load(path)
        I = Inst(P)
        fl = failures(prop, I, k)
        ok &= bool(fl)
        r2, feas, rhs = fl[0]
        sets = list(I.far_sets(range(I.n), k, r2))
        print(f"{os.path.basename(path)}: property {prop} fails at {len(fl)} radii; "
              f"first r^2={r2} (k-set exists: {feas}, property says: {rhs}); "
              f"far k-sets at that radius: {len(sets)}")
        pts = ", ".join(f"({x},{y})" for x, y in P)
        sol = ", ".join(str(i) for i in sets[0]) if sets else "none"
        lines.append(f"{prop} & {len(P)} & {k} & {r2} & {sol} \\\\\n"
                     f"\\multicolumn{{5}}{{@{{}}p{{0.95\\textwidth}}@{{}}}}{{\\footnotesize\\ttfamily {pts}}} \\\\[2pt]\n")
    os.makedirs(os.path.join(ROOT, "paper", "generated"), exist_ok=True)
    with open(os.path.join(ROOT, "paper", "generated", "counterexamples.tex"), "w") as fh:
        fh.write("".join(lines))
    if not ok:
        sys.exit("some file is not a counterexample")


if __name__ == "__main__":
    main()
