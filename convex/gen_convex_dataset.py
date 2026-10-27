#!/usr/bin/env python3
"""Generate point sets in convex position: random points -> convex hull.

For every instance we sample N random integer points from a distribution,
compute their (strictly) convex hull, and keep the hull vertices. If the hull
has more than the requested n vertices we keep a uniformly random subset of n
of them, which is still in convex position. Coordinates are integers of
absolute value <= 10^6 so all geometric predicates of the solver are exact.

Distributions
  square    uniform in a square            (hull ~ (8/3) ln N vertices: few)
  disk      uniform in a disk              (hull ~ 3.4 N^{1/3})
  gauss     standard Gaussian              (hull ~ sqrt(ln N))
  annulus   uniform in a thin annulus      (many hull vertices)
  thin      uniform in a 1:20 ellipse      (long thin hulls: far-apart
                                            indices can be geometrically close)
  circle    points rounded onto a circle   (almost all points on the hull)

Output: convex/dataset/<split>/<id>.txt in the solver format "n 2 k" and a
manifest.csv per split.
  verify   small instances (n <= 14) for exhaustive cross-checking
  medium   n in {20, 30, 40}, k in 3..6, checked against the general solvers
  scaling  n up to 200, for timing
"""

import csv
import math
import os
import random
import zlib

SEED = 20260427
S = 10 ** 6
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "dataset")


def cross(o, a, b):
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def hull(pts):
    pts = sorted(set(pts))
    if len(pts) < 3:
        return pts
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


def sample(rng, dist, N):
    out = []
    while len(out) < N:
        if dist == "square":
            x, y = rng.uniform(-1, 1), rng.uniform(-1, 1)
        elif dist == "disk":
            x, y = rng.uniform(-1, 1), rng.uniform(-1, 1)
            if x * x + y * y > 1:
                continue
        elif dist == "gauss":
            x, y = rng.gauss(0, 0.3), rng.gauss(0, 0.3)
        elif dist == "annulus":
            t = rng.uniform(0, 2 * math.pi)
            rad = rng.uniform(0.97, 1.0)
            x, y = rad * math.cos(t), rad * math.sin(t)
        elif dist == "thin":
            t = rng.uniform(0, 2 * math.pi)
            rad = math.sqrt(rng.random())
            x, y = rad * math.cos(t), 0.05 * rad * math.sin(t)
        elif dist == "circle":
            t = rng.uniform(0, 2 * math.pi)
            x, y = math.cos(t), math.sin(t)
        elif dist == "cocircular":
            out.append(rng.choice(LATTICE_CIRCLE))
            continue
        if max(abs(x), abs(y)) <= 1:
            out.append((round(S * x), round(S * y)))
    return out


def convex_instance(rng, dist, n):
    """Return n points in convex position (ccw) drawn as hull vertices."""
    N = max(4 * n, 50)
    while True:
        h = hull(sample(rng, dist, N))
        if len(h) >= n:
            keep = sorted(rng.sample(range(len(h)), n))
            return [h[i] for i in keep], N, len(h)
        N *= 2


def write(split, rows, iid, pts, k, meta):
    d = os.path.join(ROOT, split, "instances")
    os.makedirs(d, exist_ok=True)
    pts = pts[:]
    random.Random(zlib.crc32(iid.encode())).shuffle(pts)  # input order != hull order
    with open(os.path.join(d, iid + ".txt"), "w") as fh:
        fh.write(f"{len(pts)} 2 {k}\n")
        for x, y in pts:
            fh.write(f"{x} {y}\n")
    rows.append({"id": iid, "n": len(pts), "k": k, **meta})


DISTS = ["square", "disk", "gauss", "annulus", "thin", "circle", "cocircular"]

# All 324 integer points on the circle x^2 + y^2 = R^2, R = 5*13*17*29: every
# four of them are cocircular, the most degenerate input for in-circle tests.
R_LATTICE = 5 * 13 * 17 * 29
LATTICE_CIRCLE = [(x, y) for x in range(-R_LATTICE, R_LATTICE + 1)
                  for y in [math.isqrt(R_LATTICE ** 2 - x * x)]
                  if y * y == R_LATTICE ** 2 - x * x
                  for y in ({y, -y})]

# Largest hull size requested per distribution. Hulls of uniform-square and
# Gaussian samples grow only like log N and sqrt(log N), so large convex sets
# cannot be drawn from them in reasonable time.
MAX_N = {"square": 30, "gauss": 14}


def main():
    rng = random.Random(SEED)
    plans = {
        "verify": [(n, k) for n in (8, 10, 12, 14) for k in range(3, 8) if k <= n],
        "medium": [(n, k) for n in (20, 30, 40) for k in (3, 4, 5, 6)],
        "scaling": [(n, k) for n in (25, 50, 75, 100, 150, 200) for k in (4, 8, 16)],
    }
    reps = {"verify": 2, "medium": 1, "scaling": 1}
    for split, plan in plans.items():
        rows = []
        for dist in DISTS:
            for n, k in plan:
                if n > MAX_N.get(dist, n):
                    continue
                for rep in range(reps[split]):
                    pts, N, hsize = convex_instance(rng, dist, n)
                    iid = f"{dist}_n{n}_k{k}_r{rep}"
                    write(split, rows, iid, pts, k,
                          {"dist": dist, "samples": N, "hull_size": hsize})
        with open(os.path.join(ROOT, split, "manifest.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"{split}: {len(rows)} instances")


if __name__ == "__main__":
    main()
