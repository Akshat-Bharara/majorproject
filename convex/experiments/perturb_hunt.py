#!/usr/bin/env python3
"""Local search around known counterexamples (to weaker greedy rules).

Streams random perturbations of every instance in results/counterexamples/
and of the raw V5 failures in results/hunts/v5_*.log to stdout, in the
v5check input format. Perturbations: coordinate jitter at several scales,
deleting points, and inserting midpoints of hull edges pushed slightly
outward (the solver takes the convex hull, so every output is valid).

usage: perturb_hunt.py SEED COUNT | build/v5check --v4
"""

import glob
import os
import random
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def seeds():
    out = []
    for path in glob.glob(os.path.join(ROOT, "results", "counterexamples", "*.txt")):
        tok = open(path).read().split()
        n, k = int(tok[0]), int(tok[2])
        out.append(([(int(tok[3 + 2 * i]), int(tok[4 + 2 * i])) for i in range(n)], k))
    for path in glob.glob(os.path.join(ROOT, "results", "hunts", "v5_*.log")):
        lines = open(path).read().splitlines()
        for i, line in enumerate(lines):
            m = re.match(r"FAIL n (\d+) k (\d+)", line)
            if m:
                n, k = int(m.group(1)), int(m.group(2))
                pts = [tuple(map(int, l.split())) for l in lines[i + 1:i + 1 + n]]
                out.append((pts, k))
    return out


def perturb(rng, pts):
    pts = list(pts)
    span = max(max(abs(x), abs(y)) for x, y in pts)
    for _ in range(rng.randint(1, 4)):
        op = rng.random()
        if op < 0.5:
            i = rng.randrange(len(pts))
            s = span * rng.choice([1e-4, 1e-3, 1e-2, 3e-2])
            x, y = pts[i]
            pts[i] = (round(x + rng.gauss(0, s)), round(y + rng.gauss(0, s)))
        elif op < 0.75 and len(pts) > 6:
            pts.pop(rng.randrange(len(pts)))
        else:
            i = rng.randrange(len(pts))
            (x1, y1), (x2, y2) = pts[i], pts[(i + 1) % len(pts)]
            mx, my = (x1 + x2) / 2, (y1 + y2) / 2
            f = 1 + rng.choice([1e-3, 1e-2, 3e-2])
            pts.insert(i + 1, (round(mx * f), round(my * f)))
    return pts


def main():
    seed, count = int(sys.argv[1]), int(sys.argv[2])
    rng = random.Random(seed)
    base = seeds()
    for _ in range(count):
        pts, k = rng.choice(base)
        q = perturb(rng, pts)
        kk = max(3, min(k + rng.choice([-1, 0, 0, 1]), len(q) - 1))
        sys.stdout.write(f"{len(q)} 2 {kk}\n" + "".join(f"{x} {y}\n" for x, y in q))


if __name__ == "__main__":
    main()
