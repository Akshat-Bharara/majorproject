# Exact k-dispersion for points in convex position

Code, data and experiments for the paper in `../paper/` (*Exact k-Dispersion for
Points in Convex Position: A Capped Delaunay Dynamic Program and the Limits of
Greedy Structure*).

Given `n` points in strictly convex position and `k`, find `k` points whose
minimum pairwise distance is maximum. Everything is exact: coordinates are
integers with `|x|,|y| <= 2^29`, and the only geometric predicates are
comparisons of 64-bit squared distances and a 128-bit in-circle determinant.

## Quick start

```bash
make            # builds bin/kdisp, bin/bench, bin/hunt, bin/selftest, bin/fmm
make test       # randomised self-test + re-verification of all counterexamples
bin/kdisp < data/tightpair_k4.txt
```

Input format: `n k`, then `n` lines `x y` (counterclockwise order).
`bin/kdisp --algo capped|certified|triple|brute` selects the solver
(default `certified`); `bin/kdisp --mis R2` prints a maximum subset with
pairwise squared distances at least `R2`.

## Files

| File | Content |
|---|---|
| `kdisp.hpp` | header-only library: capped Delaunay DP (`CappedDP`), search window, certified algorithm, maximum r-far subset, baselines (uncapped DP, O(n^4) triple-state DP, exhaustive search), seeded greedies |
| `kdisp_cli.cpp` | command-line solver |
| `bench.cpp` | in-process benchmark, one CSV row per instance |
| `hunt.cpp`, `randgen.hpp` | randomised/exhaustive tests of the structural properties of Section 6 |
| `selftest.cpp` | all solvers against exhaustive search on random (also degenerate) instances |
| `gen_convex_dataset.py` | instance generator (random cloud -> strict convex hull -> n vertices); deterministic |
| `run_convex_experiments.py` | runs every experiment; `make_tables.py` turns `results/` into the tables/plots in `../paper/generated/` |
| `experiments/exactcheck.py` | independent pure-Python exhaustive checks (shares no code with the solvers) |
| `experiments/check_counterexamples.py` | re-verifies every explicit instance of the paper (`data/counterexamples.json`) |
| `experiments/verify_python.py` | re-checks the `verify` split with Python brute force and the FMM solver of the repository root |
| `data/` | all instances (`verify`, `medium`, `scaling`, `ksweep`), counterexamples |
| `results/` | raw outputs behind the numbers in the paper |

## Reproducing the paper

```bash
python3 run_convex_experiments.py          # about an hour on one core; --quick for a short run
cd ../paper && pdflatex paper && pdflatex paper
```

Timings are sequential single-core measurements; do not run other jobs at the same time.

## Main findings

* Decision in `O(n^2 k + n^3/w + T_r log k)` time, optimisation in `O(n^3 log k log n)`.
* For `k = 3` every r-far triple can be turned into one of the form
  `{q, a(q), b(q)}` or `{q, a(q), z}` (`a(q)`, `b(q)`: first points at distance
  `>= r` on either side of `q`).
* For `k >= 4` this structure is gone. `data/tightpair_k4.txt` is an 8-point set
  in which the unique optimal 4-set contains no pair `{q, a(q)}`; similar sets
  for `k = 5, 6, 7` are in `data/counterexamples.json`. This property had passed
  more than 55 million exact random tests (`results/hunt_*.txt`).
