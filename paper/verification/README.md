# Independent verification scripts

These scripts are a pure-Python reimplementation that uses exact (arbitrary-precision) integers. They check:

- the capped Delaunay DP (Algorithm 1) and the optimisation algorithm built on it;
- the search window (Lemma 10);
- the counterexamples in Table 3;
- the greedy and tight-pair properties of Section 6.

Each check is made against exhaustive search.

Run them from this directory:

```bash
python3 -I dispersion_check.py
```

`big.py` also cross-checks against the repository's general FMM solver, so build that solver first:

```bash
g++ -O2 -std=c++17 -o fmm ../../main.cpp ../../dispersion.cpp ../../clique.cpp ../../matrix.cpp
```

```bash
python3 -I big.py "$PWD"
```

```bash
python3 -I gk.py "$PWD"
```

`results.txt` holds the output of the runs quoted in Section 7 of the paper.
