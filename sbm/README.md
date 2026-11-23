# Checking the original paper (Singireddy, Basappa, Mitchell, DAM 386 (2026) 205–216)

Re-implementation and verification of the three algorithms of the paper, and
the head-to-head comparison with our algorithms used in the presentation.

| File | Content |
|---|---|
| `sbm.hpp` | the two-way search tree (Sec. 3, Claim 2 / Lemma 3), Algorithm 1, the DP of Sec. 4 (Thm 6) |
| `sbm_check.py` | independent Python version of the search tree, run on the explicit instances |
| `sbm_test.cpp` | search tree / Algorithm 1 / DP against exhaustive search on random instances |
| `sbm_hunt.cpp` | simulated annealing that looks for inputs on which the search tree fails |
| `sbm_approx.cpp` | the k = 3 approximation (Sec. 5): ratio on random and adversarial inputs |
| `sbm_bench.cpp` | timing: SBM DP and FPT vs our capped DP and certified algorithm |
| `results/` | raw outputs |

## Findings

* **Search tree (Claim 2):** correct on 532,833 random decision tests and 13,951
  optimisations; ~1,400 annealing runs found no failing input. The proof in the
  paper is a picture argument. The *one-sided* variant (always extend in the
  same direction, or the search seeded with {q, a(q)}) is false: the 8-point set
  P* (`../convex/data/tightpair_k4.txt`). So any proof must use both directions.
  Under the arc property the claim holds (`../arc/NOTES.md`), with no branching.
* **Algorithm 1 as printed** (line 23 removes e <= 2r, i.e. also the tested
  median) returns a wrong value on 5,288 of 8,967 random tests. The text's
  version (keep e >= 2r) is correct. One-character fix: `<=` -> `<`.
* **DP (Theorem 6):** correct on 13,951 random tests; O(n^4 k^2). Our capped DP is
  12x-37x faster for n = 20..60, and the gap grows linearly in n.
* **k = 3 approximation:** the proven ratio 1/(2 sqrt 2) ~ 0.354 holds; the worst
  ratio we observed is 0.61 (annealing) and 0.73 (20,000 random inputs).
