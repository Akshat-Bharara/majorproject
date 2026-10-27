#!/usr/bin/env python3
"""Generate the deterministic cross-verification dataset.

Each instance is written in the solver's input format
    n d k
    x_1 ... x_d     (n lines)
to verify/dataset/instances/<id>.txt, and listed in verify/dataset/manifest.csv.

Families
  uniform    integer points uniform in [0, 1000]^d
  clustered  Gaussian clusters around 3-4 random centres
  grid       sub-sample of an integer lattice (many equal distances / ties)
  duplicate  points drawn from a small pool, so coincident points occur
  collinear  points on a line in the plane (d = 2)
  circle     points rounded onto a circle (nearly equal distances)
  large      uniform coordinates of magnitude up to 1e9 (exercises 128-bit
             squared distances)
"""

import csv
import math
import os
import random

SEED = 20251101
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "dataset")

# Largest n per k, chosen so the exhaustive oracle and the FMM solver with the
# classical matrix-multiplication backend both finish quickly.
N_RANGE = {2: (8, 40), 3: (8, 35), 4: (8, 28), 5: (8, 24),
           6: (8, 20), 7: (9, 18), 8: (10, 16), 9: (10, 15)}


def uniform(rng, n, d, lo=0, hi=1000):
    return [[rng.randint(lo, hi) for _ in range(d)] for _ in range(n)]


def clustered(rng, n, d):
    centres = [[rng.uniform(0, 1000) for _ in range(d)]
               for _ in range(rng.randint(3, 4))]
    pts = []
    for _ in range(n):
        c = rng.choice(centres)
        pts.append([round(rng.gauss(x, 40)) for x in c])
    return pts


def grid(rng, n, d):
    side = max(2, math.ceil(n ** (1.0 / d)) + 1)
    cells = [[]]
    for _ in range(d):
        cells = [c + [v] for c in cells for v in range(side)]
    return [[10 * v for v in c] for c in rng.sample(cells, n)]


def duplicate(rng, n, d):
    pool = uniform(rng, max(3, n // 2), d, 0, 50)
    return [list(rng.choice(pool)) for _ in range(n)]


def collinear(rng, n, d):
    a, b = rng.randint(-5, 5), rng.randint(-50, 50)
    xs = rng.sample(range(-500, 500), n)
    return [[x, a * x + b] for x in xs]


def circle(rng, n, d):
    r = 1000
    off = rng.uniform(0, 2 * math.pi)
    return [[round(r * math.cos(off + 2 * math.pi * i / n + rng.uniform(-0.02, 0.02))),
             round(r * math.sin(off + 2 * math.pi * i / n + rng.uniform(-0.02, 0.02)))]
            for i in range(n)]


def large(rng, n, d):
    return uniform(rng, n, d, -10 ** 9, 10 ** 9)


FAMILIES = [
    ("uniform", uniform, (1, 2, 3)),
    ("clustered", clustered, (2, 3)),
    ("grid", grid, (2, 3)),
    ("duplicate", duplicate, (2,)),
    ("collinear", collinear, (2,)),
    ("circle", circle, (2,)),
    ("large", large, (2, 3)),
]


def main():
    rng = random.Random(SEED)
    inst_dir = os.path.join(ROOT, "instances")
    os.makedirs(inst_dir, exist_ok=True)
    for f in os.listdir(inst_dir):
        os.remove(os.path.join(inst_dir, f))

    rows = []
    for fam, gen, dims in FAMILIES:
        for d in dims:
            for k in range(2, 10):
                lo, hi = N_RANGE[k]
                for n in sorted({lo, (lo + hi) // 2, hi}):
                    pts = gen(rng, n, d)
                    iid = f"{fam}_d{d}_k{k}_n{n}"
                    with open(os.path.join(inst_dir, iid + ".txt"), "w") as fh:
                        fh.write(f"{n} {d} {k}\n")
                        for p in pts:
                            fh.write(" ".join(map(str, p)) + "\n")
                    rows.append({"id": iid, "family": fam, "n": n, "d": d, "k": k})

    with open(os.path.join(ROOT, "manifest.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["id", "family", "n", "d", "k"])
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} instances to {inst_dir}")


if __name__ == "__main__":
    main()
