#!/usr/bin/env python3
"""Benchmark instances for k-dispersion in convex position.

Every instance is produced the same way: sample a cloud of random integer
points from a distribution, take its strictly convex hull, and keep n of the
hull vertices (a subset of a set in convex position is in convex position).

    python3 gen_convex_dataset.py [outdir=data]

Splits
    verify   n <= 14            exhaustive search is feasible
    medium   20 <= n <= 40
    scaling  25 <= n <= 400     running-time experiments
    ksweep   n = 200, k = 4..64

File format: "n k" on the first line, then n lines "x y" (counterclockwise).
The generator is deterministic (fixed seeds).
"""
import math
import os
import sys

import numpy as np

R = 1_000_000


def prefilter(pts):
    """Exact speed-up: drop the points strictly inside the convex polygon spanned
    by the extreme points in 64 directions (they cannot be hull vertices)."""
    a = np.asarray(pts, dtype=np.int64)
    if len(a) < 5000:
        return a.tolist()
    ang = 2 * math.pi * np.arange(64) / 64
    ext = {int(np.argmax(a[:, 0] * math.cos(t) + a[:, 1] * math.sin(t))) for t in ang}
    poly = strict_hull(a[sorted(ext)].tolist())
    if len(poly) < 3:
        return a.tolist()
    inside = np.ones(len(a), dtype=bool)
    for i in range(len(poly)):
        (x0, y0), (x1, y1) = poly[i], poly[(i + 1) % len(poly)]
        inside &= (x1 - x0) * (a[:, 1] - y0) - (y1 - y0) * (a[:, 0] - x0) > 0
    return a[~inside].tolist()


def strict_hull(pts):
    """Andrew's monotone chain on exact Python integers; ccw, no collinear points."""
    pts = sorted(set(map(tuple, prefilter(pts))))
    if len(pts) < 3:
        return pts

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def lattice_circle():
    """All 972 lattice points on x^2 + y^2 = 1185665^2 (1185665 = 5*13*17*29*37)."""
    n = 5 * 13 * 17 * 29 * 37
    out = []
    for x in range(0, n + 1):
        y2 = n * n - x * x
        y = math.isqrt(y2)
        if y * y == y2:
            for sx in (1, -1):
                for sy in (1, -1):
                    out.append((sx * x, sy * y))
    return sorted(set(out))


_LATTICE = None


def cloud(dist, n, rng):
    """A point cloud whose hull has (with high probability) at least n vertices."""
    global _LATTICE
    if dist == "square":
        m = 400_000
        return rng.integers(-R, R + 1, size=(m, 2))
    if dist == "gauss":
        m = 400_000
        return np.rint(rng.normal(0, R / 6, size=(m, 2))).astype(np.int64)
    if dist == "disk":
        m = int(max(2000, 2.2 * (n / 3.38) ** 3))
        r = R * np.sqrt(rng.random(m))
        t = 2 * math.pi * rng.random(m)
        return np.rint(np.stack([r * np.cos(t), r * np.sin(t)], 1)).astype(np.int64)
    if dist == "annulus":
        m = int(max(2000, 2.2 * (n / 3.38) ** 3))
        r = R * np.sqrt(0.9025 + 0.0975 * rng.random(m))  # radius in [0.95R, R]
        t = 2 * math.pi * rng.random(m)
        return np.rint(np.stack([r * np.cos(t), r * np.sin(t)], 1)).astype(np.int64)
    if dist == "thin":  # uniform in a 1:20 ellipse
        m = int(max(2000, 2.2 * (n / 3.38) ** 3))
        r = np.sqrt(rng.random(m))
        t = 2 * math.pi * rng.random(m)
        return np.rint(np.stack([R * r * np.cos(t), R / 20 * r * np.sin(t)], 1)).astype(np.int64)
    if dist == "circle":  # rounded points on a circle
        t = 2 * math.pi * rng.random(3 * n + 10)
        return np.rint(np.stack([R * np.cos(t), R * np.sin(t)], 1)).astype(np.int64)
    if dist == "cocircular":  # exactly cocircular lattice points
        if _LATTICE is None:
            _LATTICE = lattice_circle()
        return np.array(_LATTICE, dtype=np.int64)
    if dist == "twochain":  # thin lens bounded by two flat parabolic arcs
        m = 6 * n + 20
        x = (2 * rng.random(m) - 1) * R
        s = 1 - (x / R) ** 2
        h = R / 50 * (0.5 + rng.random())
        y = np.where(rng.random(m) < 0.5, h * s, -h * (0.4 + rng.random()) * s)
        return np.rint(np.stack([x, y], 1)).astype(np.int64)
    raise ValueError(dist)


def instance(dist, n, rng):
    hull = strict_hull(cloud(dist, n, rng).tolist())
    if len(hull) < n:
        return None
    keep = sorted(rng.choice(len(hull), size=n, replace=False).tolist())
    return [hull[i] for i in keep]


def write(path, pts, k):
    with open(path, "w") as f:
        f.write(f"{len(pts)} {k}\n")
        for x, y in pts:
            f.write(f"{x} {y}\n")


ALL = ["square", "gauss", "disk", "annulus", "thin", "circle", "cocircular", "twochain"]
LARGE = ["disk", "annulus", "thin", "circle", "cocircular", "twochain"]  # hulls can be large


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "data"
    count = 0

    def emit(split, dist, n, k, s, rng):
        nonlocal count
        pts = instance(dist, n, rng)
        if pts is None:
            return
        os.makedirs(os.path.join(out, split), exist_ok=True)
        write(os.path.join(out, split, f"{dist}_n{n:03d}_k{k:02d}_s{s}.txt"), pts, k)
        count += 1

    # verify: 40 instances per distribution, 6 <= n <= 14, 3 <= k <= min(n-1, 8)
    for d, dist in enumerate(ALL):
        rng = np.random.default_rng(1000 + d)
        for s in range(40):
            n = int(rng.integers(6, 15))
            k = int(rng.integers(3, min(n - 1, 8) + 1))
            emit("verify", dist, n, k, s, rng)
    # medium
    for d, dist in enumerate(ALL):
        rng = np.random.default_rng(2000 + d)
        for s, (n, k) in enumerate([(20, 4), (20, 8), (25, 5), (30, 4), (30, 12), (35, 5), (40, 4),
                                    (40, 5), (40, 8), (40, 12)]):
            if dist == "gauss" and n > 20:
                continue  # Gaussian hulls are too small
            if dist == "square" and n > 30:
                continue
            emit("medium", dist, n, k, s, rng)
    # scaling
    for d, dist in enumerate(LARGE):
        rng = np.random.default_rng(3000 + d)
        for n in (25, 50, 75, 100, 150, 200, 300, 400):
            for k in (4, 8, 16, 32):
                if k >= n:
                    continue
                for s in range(3 if n <= 200 else 1):
                    emit("scaling", dist, n, k, s, rng)
    # k sweep at n = 200
    for d, dist in enumerate(LARGE):
        rng = np.random.default_rng(4000 + d)
        for k in (4, 8, 16, 32, 64):
            emit("ksweep", dist, 200, k, 0, rng)
    print(f"wrote {count} instances to {out}/")


if __name__ == "__main__":
    main()
