# majorproject

Exact k-dispersion (select k of n points maximizing the minimum pairwise
Euclidean distance) via the FMM-based k-clique reduction used in
K. Chen, A. Dumitrescu, *A Couple of Simple Algorithms for k-Dispersion*
(arXiv:2511.00692), Section 3 / Theorem 2.

## Algorithm

1. Compute all squared pairwise distances (integral coordinates, exact).
2. The optimum is one of the pairwise distances. Binary search over the sorted
   distinct values r: build the threshold graph G_r (edge iff dist >= r);
   a k-set with min distance >= r exists iff G_r contains a K_k.
3. K_k detection (Nešetřil–Poljak / Eisenbrand–Grandoni): split
   k = k1 + k2 + k3, let V_i be the k_i-cliques of G_r, connect two clique
   nodes iff their union is a clique. G_r has a K_k iff this 3-partite graph
   has a triangle, found from P = M12 · M23 checked against M13.

The split follows the rectangular products used in the proof of Theorem 2
(`clique_split` in `dispersion.cpp`):

| k | (k1, k2, k3) | exponent e(k) |
|---|--------------|---------------|
| 4 | (1, 1, 2) | ω(1,1,2) < 3.251 |
| 5 | (2, 1, 2) | ω(2,1,2) < 4.086 |
| 6 | (2, 2, 2) | 2ω < 4.744 |
| 7 | (2, 3, 2) | ω(2,3,2) < 5.590 |
| 8 | (2, 2, 4) | 2ω(1,1,2) < 6.502 |

k = 2, 3 use (0,1,1) and (1,1,1). For k ≥ 9 (not covered by the paper) the
balanced split ⌊k/3⌋, ⌈(k−1)/3⌉, ⌈k/3⌉ is used.

## Matrix multiplication backend

`multiply_blocked` in `matrix.cpp` is an exact, cache-blocked classical
product. The reduction is exact, but the asymptotic bounds above only hold if
this routine is replaced by a fast (rectangular) matrix multiplication
algorithm; the paper does not specify one.

## Build and run

```
g++ -O2 -std=c++17 *.cpp -o dispersion
./dispersion < input.txt
```

Input: `n d k` followed by n lines of d integer coordinates.

Verification against an exhaustive solver: `verify/` (see `report/report.pdf`).

## Points in convex position (`convex/`)

A specialised exact solver for planar points in convex position, with the
research write-up in `convex/paper/paper.pdf`.

- `convex_dispersion.{hpp,cpp}`: the capped Delaunay DP (decision in
  O(n^3 log k) time, O(n^2 k) space, exact integer predicates), the
  Tkachenko–Wang DP (STACS 2025) as a baseline, a seeded greedy, the
  certified algorithm (greedy plus one exact decision; always exact), and an
  O(n^2 log n) algorithm for k = 4 (`--algo=tight4`) based on the oriented
  tight-pair theorem.
- Structural results (proved in the paper): every feasible instance has a
  solution containing a tight pair {q, a(q)} or {b(q), q}. With one fixed
  orientation this holds for k = 3 and fails for k = 4
  (`results/counterexamples/C1_n8_k4.txt`).
- `gen_convex_dataset.py`: random points → convex hull → convex instances.
- `run_convex_experiments.py`: cross-checks all solvers (including the
  exhaustive oracle and the FMM solver) and times them.
- `experiments/`: tests of structural conjectures, counterexample search,
  minimisation and independent re-verification.

```
g++ -O2 -std=c++17 convex/geometry.cpp convex/convex_dispersion.cpp \
    convex/main_convex.cpp -o convex_dispersion
./convex_dispersion [--algo=capped|certified|tw|tight4] [--hull] [--stats] < points.txt
```

The input format is the same (`n 2 k`). The points must be in strictly convex
position; with `--hull` the solver runs on the hull vertices instead.
