// Re-implementation of the algorithms of
//   V.R. Singireddy, M. Basappa, J.S.B. Mitchell, "Algorithms for k-dispersion
//   for points in convex position in the plane", Discrete Applied Math. 386
//   (2026) 205-216   [SBM].
// * sbm_decide:   the two-way search tree of Section 3 (Claim 2 / Lemma 3)
// * sbm_optimize: Algorithm 1 (median halving over the candidate set)
// * SbmDP:        the dynamic program of Section 4 (Theorem 6), O(n^4 k^2)
// All in exact integer arithmetic on squared distances (D2 = (2r)^2).
#pragma once
#include "../convex/kdisp.hpp"
#include <algorithm>
#include <vector>

namespace sbm {
using kd::i64;
using kd::Instance;

// ---------------- Section 3: two-way search tree --------------------------
struct Search {
    const Instance& I;
    int k;
    i64 D2;
    long long nodes = 0;
    std::vector<int> S, best;
    Search(const Instance& inst, int k_, i64 D2_) : I(inst), k(k_), D2(D2_) {}
    bool far_all(int v) const {
        for (int u : S) if (I.d(u, v) < D2) return false;
        return true;
    }
    // free vertices: ccw offsets lo+1 .. hi-1 from the root S[0]
    bool rec(int lo, int hi) {
        ++nodes;
        if ((int)S.size() == k) { best = S; return true; }
        const int n = I.n, q = S[0];
        for (int side = 0; side < 2; ++side) {
            const bool ccw = side == 0;
            for (int s = 0; s < hi - lo - 1; ++s) {
                const int t = ccw ? lo + 1 + s : hi - 1 - s;
                const int v = (q + t) % n;
                if (!far_all(v)) continue;
                S.push_back(v);
                const bool ok = ccw ? rec(t, hi) : rec(lo, t);
                S.pop_back();
                if (ok) return true;
                break; // only the first feasible vertex on this side
            }
        }
        return false;
    }
    bool run(std::vector<int>* witness = nullptr) {
        for (int v = 0; v < I.n; ++v) {
            S.assign(1, v);
            if (rec(0, I.n)) { if (witness) { *witness = best; std::sort(witness->begin(), witness->end()); } return true; }
        }
        return false;
    }
};
inline bool sbm_decide(const Instance& I, int k, i64 D2, std::vector<int>* w = nullptr) {
    Search s(I, k, D2);
    return s.run(w);
}
// Algorithm 1 as described in the text (YES keeps elements >= median, NO keeps < median).
inline i64 sbm_optimize(const Instance& I, int k, int* calls = nullptr) {
    std::vector<i64> X = I.candidates();
    int c = 0;
    while (X.size() >= 2) {
        const i64 m = X[X.size() / 2];
        ++c;
        const bool yes = sbm_decide(I, k, m);
        std::vector<i64> Y;
        for (i64 e : X) if (yes ? e >= m : e < m) Y.push_back(e);
        X.swap(Y);
    }
    if (calls) *calls = c;
    return X[0];
}

// ---------------- Section 4: dynamic program (Theorem 6) ------------------
// phi(i, j, l; K): best min-edge of the part of a solution that lies left of
// the directed chord p_i -> p_j (i.e. strictly inside the ccw arc from i to j
// ... see below), with K further points to place there, given that the
// Delaunay triangle on the other side of (i,j) has apex l (l = n: no apex,
// (i,j) is a hull edge of the solution).  Candidates p_m: strictly between,
// not strictly inside the circumcircle of (p_i, p_j, p_l).
//   phi(i,j,l;0) = |p_i p_j|
//   phi(i,j,l;K) = max_{m, 0<=K'<=K-1} min(|p_i p_m|, |p_m p_j|,
//                        phi(i,m,j;K'), phi(m,j,i;K-1-K'))
// Answer: max_{i,j} min(|p_i p_j|, phi(i,j,*;k-2)).
// Points are in ccw order; the region "left of p_i -> p_j" for a chord whose
// other side is already used is the ccw arc from j to i; we orient so that
// the open arc (i, j) in ccw order is the region to fill: phi(i,j,..) fills
// the ccw arc strictly between i and j.  Values are squared distances; -1 = infeasible.
struct SbmDP {
    const Instance& I;
    int n, k;
    std::vector<i64> memo; // (i, j, l, K)
    long long work = 0;
    SbmDP(const Instance& inst, int k_) : I(inst), n(inst.n), k(k_) {}
    size_t id(int i, int j, int l, int K) const { return (((size_t)i * n + j) * (n + 1) + l) * (k - 1) + K; }
    // ccw offset helpers
    int off(int a, int b) const { return ((b - a) % n + n) % n; }
    i64 phi(int i, int j, int l, int K) {
        if (K == 0) return I.d(i, j);
        i64& m = memo[id(i, j, l, K)];
        if (m != -2) return m;
        i64 best = -1;
        const int g = off(i, j);
        if (g - 1 >= K) {
            for (int t = 1; t < g; ++t) {
                const int x = (i + t) % n;
                ++work;
                // x must not be strictly inside circumcircle of (p_i, p_l, p_j)... orientation:
                // (p_j, p_i, p_l)? the apex l lies on the other side of chord (i,j): ccw order i, x, j, l.
                if (l < n && kd::incircle(I.p[i], I.p[x], I.p[j], I.p[l]) > 0) continue; // l inside circle(i,x,j)  <=> x inside circle(i,j,l)
                const i64 e = std::min(I.d(i, x), I.d(x, j));
                if (e <= best) continue;
                const int g1 = t - 1, g2 = g - t - 1;
                for (int K1 = 0; K1 <= K - 1; ++K1) {
                    const int K2 = K - 1 - K1;
                    if (K1 > g1 || K2 > g2) continue;
                    const i64 a = phi(i, x, j, K1);
                    if (std::min(e, a) <= best) continue;
                    const i64 b = phi(x, j, i, K2);
                    const i64 v = std::min({e, a, b});
                    if (a < 0 || b < 0) continue;
                    best = std::max(best, v);
                }
            }
        }
        return m = best;
    }
    i64 solve() {
        if (k == 2) { i64 b = 0; for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) b = std::max(b, I.d(i, j)); return b; }
        memo.assign((size_t)n * n * (n + 1) * (k - 1), -2);
        i64 best = -1;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                if (i == j) continue;
                // hull edge (j -> i) of the solution; the solution fills the ccw arc (i, j)
                const i64 v = phi(i, j, n, k - 2);
                if (v >= 0) best = std::max(best, std::min(v, I.d(i, j)));
            }
        return best;
    }
};
} // namespace sbm
