#!/usr/bin/env python3
"""Independent brute-force verification (shares no code with the C++ library).

    python3 -I pybrute.py <data-dir> <results/small_cpp.txt>

For every exported small instance:
  * recompute the optimum r*_k for every k by enumerating all k-subsets and
    compare with the values written by the C++ solvers (O2 / adaptive);
  * check every reported witness (k distinct indices, minimum distance = r*);
  * re-check the unimodality flag with an independent implementation;
  * check Corollary A2 directly: at every candidate radius r at which the arc
    property holds (own check), some greedy chain q, a(q), a(a(q)), ... is
    r-far (verified pairwise) if and only if an r-far k-set exists.
Exact Python integers only.
"""
import sys
from itertools import combinations


def read_points(path):
    with open(path) as f:
        n = int(f.readline())
        return [tuple(map(int, f.readline().split())) for _ in range(n)]


def d2(P, i, j):
    return (P[i][0] - P[j][0]) ** 2 + (P[i][1] - P[j][1]) ** 2


def optimum_all_k(P):
    n = len(P)
    D = [[d2(P, i, j) for j in range(n)] for i in range(n)]
    best = {}
    for k in range(2, n + 1):
        b = -1
        for S in combinations(range(n), k):
            m = min(D[a][c] for a, c in combinations(S, 2))
            if m > b:
                b = m
        best[k] = b
    return best


def unimodal(P):
    n = len(P)
    for q in range(n):
        row = [d2(P, q, (q + t) % n) for t in range(1, n)]
        m = max(row)
        i = row.index(m)
        j = len(row) - 1 - row[::-1].index(m)
        if j - i > 1:
            return False
        if any(row[t] >= row[t + 1] for t in range(i)) or any(row[t] <= row[t + 1] for t in range(j, len(row) - 1)):
            return False
    return True


def arc_property(P, r2):
    n = len(P)
    for q in range(n):
        far = [d2(P, q, (q + t) % n) >= r2 for t in range(1, n)]
        runs = sum(1 for t in range(len(far)) if far[t] and (t == 0 or not far[t - 1]))
        if runs > 1:
            return False
    return True


def greedy_chain_exists(P, k, r2):
    n = len(P)
    for q in range(n):
        chain, x, total = [q], q, 0
        for _ in range(k - 1):
            step = next((t for t in range(1, n) if d2(P, x, (x + t) % n) >= r2), None)
            if step is None or total + step >= n:
                chain = None
                break
            total += step
            x = (x + step) % n
            chain.append(x)
        if chain and all(d2(P, a, b) >= r2 for a, b in combinations(chain, 2)):
            return True
    return False


def main():
    data, results = sys.argv[1], sys.argv[2]
    rows, flags = {}, {}
    for line in open(results):
        f = line.split()
        if f[1] == "flags":
            flags[f[0]] = int(f[2])
        else:
            rows.setdefault(f[0], []).append((int(f[1]), f[2], int(f[3]), list(map(int, f[4:]))))
    bad = checked = a2 = 0
    for name in sorted(flags):
        P = read_points(data + "/" + name)
        n = len(P)
        best = optimum_all_k(P)
        if unimodal(P) != bool(flags[name]):
            bad += 1
            print("UNIMODAL FLAG MISMATCH", name)
        for k, method, r2, S in rows.get(name, []):
            checked += 1
            ok = r2 == best[k] and len(set(S)) == k and all(0 <= v < n for v in S)
            ok = ok and min(d2(P, a, b) for a, b in combinations(S, 2)) == best[k]
            if not ok:
                bad += 1
                print("MISMATCH", name, k, method, r2, best[k])
        cand = sorted({d2(P, i, j) for i in range(n) for j in range(i + 1, n)})
        for r2 in cand:
            if not arc_property(P, r2):
                continue
            for k in range(2, n + 1):
                a2 += 1
                if greedy_chain_exists(P, k, r2) != (best[k] >= r2):
                    bad += 1
                    print("A2 VIOLATED", name, k, r2)
    print(f"instances {len(flags)}  solver values checked {checked}  A2 checks {a2}  failures {bad}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
