#!/usr/bin/env python3
"""Timing of the k = 4 algorithm (Theorem: O(n^2 log n)) against the capped
DP (O(n^3 log n) for fixed k) on large convex instances.

Writes results/k4_scaling.csv and paper/generated/k4_scaling.dat.
"""

import csv
import math
import os
import random
import statistics
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
EXE = os.path.join(ROOT, "build", "convex_dispersion")


S = 10 ** 8   # coordinates <= 2^29 keep every predicate exact


def circle(rng, n):
    ang = sorted(rng.uniform(0, 2 * math.pi) for _ in range(n))
    return sorted({(round(S * math.cos(a)), round(S * math.sin(a))) for a in ang})


def valtr(rng, n):
    """Uniformly random convex polygon (Valtr's construction), one pass."""
    xs = sorted(rng.random() for _ in range(n))
    ys = sorted(rng.random() for _ in range(n))

    def chain(v):
        a, b = [v[0]], [v[0]]
        for x in v[1:-1]:
            (a if rng.random() < .5 else b).append(x)
        a.append(v[-1]); b.append(v[-1])
        return [a[i + 1] - a[i] for i in range(len(a) - 1)] + \
               [b[i] - b[i + 1] for i in range(len(b) - 1)]
    dx, dy = chain(xs), chain(ys)
    rng.shuffle(dy)
    vec = sorted(zip(dx, dy), key=lambda v: math.atan2(v[1], v[0]))
    x = y = 0.0
    pts = []
    for vx, vy in vec:
        pts.append((round(S * x), round(S * y)))
        x += vx; y += vy
    return pts


def twochain(rng, n):
    """Thin lens: points on two flat parabolic arcs."""
    h = rng.uniform(0.001, 0.08)
    return [(round(S * x), round(S * h * (1 - x * x) * (1 if rng.random() < .5 else -1)))
            for x in (rng.uniform(-1, 1) for _ in range(n))]


GEN = {"valtr": valtr, "twochain": twochain, "circle": circle}


def run(text, algo):
    p = subprocess.run([EXE, f"--algo={algo}", "--stats", "--hull"], input=text,
                       capture_output=True, text=True, check=True)
    out = {}
    for line in (p.stdout + p.stderr).splitlines():
        parts = line.split()
        if len(parts) == 2:
            try:
                out[parts[0]] = float(parts[1])
            except ValueError:
                pass
    return out


def main():
    rng = random.Random(4)
    sizes = [int(x) for x in (sys.argv[1:] or ["100", "200", "400", "800", "1600"])]
    capped_max = 400
    rows = []
    for n in sizes:
        for kind in ("valtr", "twochain", "circle"):
            pts = GEN[kind](rng, n)   # the solver keeps the hull vertices (--hull)
            text = f"{len(pts)} 2 4\n" + "".join(f"{x} {y}\n" for x, y in pts)
            t = run(text, "tight4")
            row = {"n": int(t["hull_size"]), "family": kind, "tight4_ms": t["solve_ms"],
                   "opt": int(t["optimal_squared_distance"])}
            if n <= capped_max:
                c = run(text, "capped")
                assert int(c["optimal_squared_distance"]) == row["opt"], (n, kind)
                row["capped_ms"] = c["solve_ms"]
            else:
                row["capped_ms"] = ""
            rows.append(row)
            print(row)
    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    with open(os.path.join(ROOT, "results", "k4_scaling.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(ROOT, "paper", "generated", "k4_scaling.dat"), "w") as fh:
        fh.write("n tight4 capped\n")
        for n in sizes:
            g = [r for r in rows if abs(r["n"] - n) <= n * 0.05 + 3]
            if not g:
                continue
            cap = [r["capped_ms"] for r in g if r["capped_ms"] != ""]
            fh.write(f"{n} {statistics.median(r['tight4_ms'] for r in g):.3f} "
                     f"{statistics.median(cap) if cap else 'nan'}\n")


if __name__ == "__main__":
    main()
