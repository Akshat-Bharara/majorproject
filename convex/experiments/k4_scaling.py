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
from greedy_stress import family  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
EXE = os.path.join(ROOT, "build", "convex_dispersion")


def circle(rng, n):
    S = 10 ** 7
    ang = sorted(rng.uniform(0, 2 * math.pi) for _ in range(n))
    pts = sorted({(round(S * math.cos(a)), round(S * math.sin(a))) for a in ang})
    return pts


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
    capped_max = 800
    rows = []
    for n in sizes:
        for kind in ("valtr", "twochain", "circle"):
            pts = circle(rng, n) if kind == "circle" else family(rng, n, kind)
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
