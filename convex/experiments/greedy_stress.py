#!/usr/bin/env python3
"""Stress test of conjecture V1 (alternating two-sided greedy) on harder inputs.

Truth: r* by exhaustive search (feasible(r) <=> r <= r*). For every candidate
r <= r* we require that some seed q and one of the two alternation orders
produce k points. Families include adversarial two-chain shapes in which the
cyclic order is very different from geometric proximity.
"""

import itertools
import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from structure import Inst, hull  # noqa: E402
from greedy_conj import v1        # noqa: E402


def family(rng, n, kind):
    S = 10 ** 5
    while True:
        if kind == "twochain":       # thin lens: two nearly flat chains
            h = rng.uniform(0.001, 0.08)
            xs = sorted(rng.uniform(-1, 1) for _ in range(n))
            raw = [(x, h * (1 - x * x) * (1 if rng.random() < .5 else -1)) for x in xs]
        elif kind == "valtr":        # uniformly random convex polygon
            xs = sorted(rng.random() for _ in range(n)); ys = sorted(rng.random() for _ in range(n))
            def chain(v):
                a, b = [v[0]], [v[0]]
                for x in v[1:-1]:
                    (a if rng.random() < .5 else b).append(x)
                a.append(v[-1]); b.append(v[-1])
                return [a[i + 1] - a[i] for i in range(len(a) - 1)] + \
                       [b[i] - b[i + 1] for i in range(len(b) - 1)]
            dx, dy = chain(xs), chain(ys); rng.shuffle(dy)
            vec = sorted(zip(dx, dy), key=lambda v: math.atan2(v[1], v[0]))
            x = y = 0; raw = []
            for vx, vy in vec:
                raw.append((x, y)); x += vx; y += vy
        elif kind == "triangle":     # points on the slightly bulged sides of a thin triangle
            raw = []
            A, B, C = (-1, 0), (1, 0), (rng.uniform(-1, 1), rng.uniform(0.02, 0.4))
            for _ in range(n):
                P, Q = rng.choice([(A, C), (C, B), (B, A)])  # clockwise sides
                t = rng.random()
                nx, ny = P[1] - Q[1], Q[0] - P[0]            # outward normal (cw order)
                bulge = 0.02 * t * (1 - t)
                raw.append((P[0] + t * (Q[0] - P[0]) + bulge * nx,
                            P[1] + t * (Q[1] - P[1]) + bulge * ny))
        elif kind == "clusters":     # clusters of points around 3-5 arcs of a circle
            c = rng.randint(3, 5)
            cen = [rng.uniform(0, 2 * math.pi) for _ in range(c)]
            raw = [(math.cos(t), rng.choice([1, 0.3]) * math.sin(t)) for t in
                   (rng.choice(cen) + rng.gauss(0, 0.15) for _ in range(n))]
        pts = hull([(round(S * x), round(S * y)) for x, y in raw])
        if len(pts) >= max(5, n - 3):
            return pts


def rstar(I, k):
    best = -1
    for c in itertools.combinations(range(I.n), k):
        m = min(I.D[a][b] for a, b in itertools.combinations(c, 2))
        if m > best:
            best = m
    return best


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    trials = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    nmax = int(sys.argv[3]) if len(sys.argv) > 3 else 14
    rng = random.Random(seed)
    tests = fails = 0
    for trial in range(trials):
        kind = rng.choice(["twochain", "valtr", "triangle", "clusters"])
        I = Inst(family(rng, rng.randint(8, nmax), kind))
        n = I.n
        k = rng.randint(3, min(8, n - 1))
        rs = rstar(I, k)
        for r2 in sorted({I.D[a][b] for a, b in itertools.combinations(range(n), 2)}):
            if r2 > rs:
                break
            tests += 1
            if not any(v1(I, q, k, r2) for q in range(n)):
                fails += 1
                print("COUNTEREXAMPLE", kind, "n", n, "k", k, "r2", r2, "r*2", rs)
                print(I.P)
                sys.stdout.flush()
                break
    print(f"seed {seed}: {tests} feasible (P,k,r) tested, {fails} failures")


if __name__ == "__main__":
    main()
