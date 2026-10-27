#!/usr/bin/env python3
"""Empirical tests of structural conjectures for k-dispersion in convex position.

All arithmetic is exact (integer coordinates, squared distances).

For a threshold r (squared: r2) and a point p_i of the ccw cyclic sequence P:
  a(i)  first point counterclockwise from p_i with |p_i p| >= r
  b(i)  first point clockwise        from p_i with |p_i p| >= r
  A(i)  the ccw arc P[a(i), b(i)]

G_k  (generalised Lemma 16 of Tkachenko-Wang; true for k = 3):
      P has k points pairwise >= r  <=>  some arc A(i) has k-1 points pairwise >= r.
C1   if a feasible k-set exists, some feasible k-set contains a pair q, a(q).
C1'  ... some feasible k-set contains q, a(q) and b(q).

Writes convex/results/structure.csv and convex/paper/generated/structure.tex,
and the smallest counterexample to G_4 found to results/counterexamples/.
"""

import collections
import csv
import itertools
import math
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def cross(o, a, b):
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def hull(pts):
    pts = sorted(set(pts))
    lo, up = [], []
    for p in pts:
        while len(lo) >= 2 and cross(lo[-2], lo[-1], p) <= 0:
            lo.pop()
        lo.append(p)
    for p in reversed(pts):
        while len(up) >= 2 and cross(up[-2], up[-1], p) <= 0:
            up.pop()
        up.append(p)
    return lo[:-1] + up[:-1]


def gen(rng, n, kind):
    S = 1000
    while True:
        if kind == "ellipse":
            b = rng.choice([1, 0.5, 0.2, 0.05])
            raw = [(math.cos(t), b * math.sin(t))
                   for t in (rng.uniform(0, 2 * math.pi) for _ in range(n))]
        elif kind == "disk":
            raw = []
            while len(raw) < 30 * n:
                x, y = rng.uniform(-1, 1), rng.uniform(-1, 1)
                if x * x + y * y <= 1:
                    raw.append((x, y))
        else:  # thin
            raw = [(rng.uniform(-1, 1), rng.uniform(-0.05, 0.05)) for _ in range(40 * n)]
        h = hull([(round(S * x), round(S * y)) for x, y in raw])
        if len(h) >= n:
            keep = sorted(rng.sample(range(len(h)), n))
            return [h[i] for i in keep]


def d2(p, q):
    return (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2


class Inst:
    def __init__(self, P):
        self.P, self.n = P, len(P)
        self.D = [[d2(a, b) for b in P] for a in P]

    def far_sets(self, idx, k, r2):
        D = self.D
        for c in itertools.combinations(idx, k):
            if all(D[a][b] >= r2 for a, b in itertools.combinations(c, 2)):
                yield c

    def feasible(self, idx, k, r2):
        if k <= 1:
            return len(idx) >= k
        return next(self.far_sets(idx, k, r2), None) is not None

    def a(self, i, r2, step=1):
        for t in range(1, self.n):
            j = (i + step * t) % self.n
            if self.D[i][j] >= r2:
                return j
        return None

    def arc(self, i, r2):
        ai, bi = self.a(i, r2, 1), self.a(i, r2, -1)
        if ai is None:
            return []
        out, j = [ai], ai
        while j != bi:
            j = (j + 1) % self.n
            out.append(j)
        return out


def main():
    rng = random.Random(11)
    stats = collections.Counter()
    smallest = None
    for trial in range(1800):
        kind = rng.choice(["ellipse", "disk", "thin"])
        n = rng.randint(6, 11)
        I = Inst(gen(rng, n, kind))
        k = rng.randint(3, min(6, n - 1))
        for r2 in sorted({I.D[a][b] for a, b in itertools.combinations(range(n), 2)}):
            lhs = I.feasible(range(n), k, r2)
            rhs = any(I.feasible(I.arc(i, r2), k - 1, r2) for i in range(n))
            stats[("G", k, "tests")] += 1
            if lhs != rhs:
                stats[("G", k, "fail")] += 1
                if k == 4 and (smallest is None or n < len(smallest[0])):
                    smallest = (I.P, r2)
            if not lhs:
                continue
            A = [I.a(i, r2, 1) for i in range(n)]
            B = [I.a(i, r2, -1) for i in range(n)]
            sets = list(I.far_sets(range(n), k, r2))
            stats[("C1", k, "tests")] += 1
            stats[("C1p", k, "tests")] += 1
            if not any(A[q] in s for s in sets for q in s):
                stats[("C1", k, "fail")] += 1
            if not any(A[q] in s and B[q] in s for s in sets for q in s):
                stats[("C1p", k, "fail")] += 1

    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    rows = []
    for name in ("G", "C1", "C1p"):
        for k in range(3, 7):
            t = stats[(name, k, "tests")]
            if t:
                rows.append({"conjecture": name, "k": k, "tests": t,
                             "failures": stats[(name, k, "fail")]})
    with open(os.path.join(ROOT, "results", "structure.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    os.makedirs(os.path.join(ROOT, "paper", "generated"), exist_ok=True)
    label = {"G": "$G_k$ (Lemma 16 generalised)", "C1": "$C_1$ (a tight pair $q,a(q)$)",
             "C1p": "$C_1'$ ($q,a(q),b(q)$ together)"}
    with open(os.path.join(ROOT, "paper", "generated", "structure.tex"), "w") as fh:
        for r in rows:
            fh.write(f"{label[r['conjecture']]} & {r['k']} & {r['tests']} & {r['failures']} \\\\\n")
    if smallest:
        P, r2 = smallest
        os.makedirs(os.path.join(ROOT, "results", "counterexamples"), exist_ok=True)
        with open(os.path.join(ROOT, "results", "counterexamples",
                               f"G4_n{len(P)}_k4.txt"), "w") as fh:
            fh.write(f"{len(P)} 2 4\n")
            for x, y in P:
                fh.write(f"{x} {y}\n")
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
