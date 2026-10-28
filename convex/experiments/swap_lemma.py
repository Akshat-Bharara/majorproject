#!/usr/bin/env python3
"""Local swap properties that would imply the tight-pair conjecture C1.

For an r-far k-set S (ccw order) and s in S with successor t:
  swap+(s):  S - t + a(s)        (a(s) lies in (s, t])
  swap-(t):  S - s + b(t)        (b(t) lies in [s, t)); pair (b(t), t) is
             cw-tight, which is C1 for the mirrored instance
Properties, for EVERY r-far k-set S:
  W+   some swap+(s) is r-far                (implies C1)
  W+-  some swap+(s) or swap-(t) is r-far
Exhaustive on small instances; prints the first failure of each.
"""

import collections
import itertools
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from structure import Inst, gen          # noqa: E402
from greedy_stress import family         # noqa: E402


def far_set(I, S, r2):
    return all(I.D[a][b] >= r2 for a, b in itertools.combinations(S, 2))


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    trials = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    rng = random.Random(seed)
    st = collections.Counter()
    shown = set()
    for _ in range(trials):
        kind = rng.choice(["twochain", "valtr", "triangle", "clusters", "ellipse", "thin"])
        n = rng.randint(7, 12)
        P = gen(rng, n, kind) if kind in ("ellipse", "thin") else family(rng, n, kind)
        I = Inst(P)
        n = I.n
        k = rng.randint(3, min(6, n - 1))
        for r2 in sorted({I.D[a][b] for a, b in itertools.combinations(range(n), 2)}):
            A = [I.a(i, r2, 1) for i in range(n)]
            B = [I.a(i, r2, -1) for i in range(n)]
            for S in I.far_sets(range(n), k, r2):
                S = list(S)
                st[("sets", k)] += 1
                plus = minus = False
                for idx in range(k):
                    s, t = S[idx], S[(idx + 1) % k]
                    if A[s] == t:
                        plus = True
                        break
                    if far_set(I, [x for x in S if x != t] + [A[s]], r2):
                        plus = True
                        break
                if not plus:
                    st[("W+fail", k)] += 1
                    for idx in range(k):
                        s, t = S[idx], S[(idx + 1) % k]
                        if B[t] == s or far_set(I, [x for x in S if x != s] + [B[t]], r2):
                            minus = True
                            break
                    if not minus:
                        st[("W+-fail", k)] += 1
                    key = ("W+" if minus else "W+-", k)
                    if key not in shown:
                        shown.add(key)
                        print("FAIL", key, kind, "n", n, "r2", r2, "S", S, P)
    for key in sorted(st):
        print(key, st[key])


if __name__ == "__main__":
    main()
