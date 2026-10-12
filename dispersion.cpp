#include "dispersion.hpp"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace fmm {

static inline int wi(int v) { return v >> 6; }
static inline std::uint64_t wm(int v) { return 1ULL << (v & 63); }

std::string uint128_to_string(unsigned __int128 x) {
    if (x == 0) return "0";
    std::string s;
    while (x) {
        s.push_back(char('0' + (x % 10)));
        x /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

FMMKDispersion::FMMKDispersion(std::vector<Point> points, int k, std::size_t gemm_block)
    : points_(std::move(points)), n_(static_cast<int>(points_.size())), k_(k),
      d_(points_.empty() ? 0 : static_cast<int>(points_[0].x.size())),
      gemm_block_(gemm_block) {
    if (n_ == 0) throw std::invalid_argument("empty point set");
    if (k_ < 2 || k_ > n_) throw std::invalid_argument("require 2 <= k <= n");
    for (const auto& p : points_) {
        if ((int)p.x.size() != d_) throw std::invalid_argument("inconsistent dimension");
    }
}

void FMMKDispersion::build_distances() {
    dist2_.assign(n_, std::vector<unsigned __int128>(n_, 0));
    for (int i = 0; i < n_; ++i) {
        for (int j = i + 1; j < n_; ++j) {
            unsigned __int128 s = 0;
            for (int t = 0; t < d_; ++t) {
                const __int128 diff = (__int128)points_[i].x[t] - (__int128)points_[j].x[t];
                const __int128 sq = diff * diff;
                s += (unsigned __int128)sq;
            }
            dist2_[i][j] = dist2_[j][i] = s;
        }
    }
}

void FMMKDispersion::build_candidates() {
    candidate_thresholds_.clear();
    candidate_thresholds_.reserve((std::size_t)n_ * (n_ - 1) / 2);
    for (int i = 0; i < n_; ++i) {
        for (int j = i + 1; j < n_; ++j) {
            candidate_thresholds_.push_back(dist2_[i][j]);
        }
    }
    std::sort(candidate_thresholds_.begin(), candidate_thresholds_.end());
    candidate_thresholds_.erase(
        std::unique(candidate_thresholds_.begin(), candidate_thresholds_.end()),
        candidate_thresholds_.end());
}

ThresholdGraph FMMKDispersion::build_graph(unsigned __int128 r2) const {
    ThresholdGraph G;
    G.n = n_;
    const int W = (n_ + 63) / 64;
    G.adj.assign(n_, std::vector<std::uint64_t>(W, 0));

    for (int i = 0; i < n_; ++i) {
        for (int j = i + 1; j < n_; ++j) {
            if (dist2_[i][j] >= r2) {
                G.adj[i][wi(j)] |= wm(j);
                G.adj[j][wi(i)] |= wm(i);
            }
        }
    }
    return G;
}

FMMKDispersion::FeasibleResult FMMKDispersion::fmm_clique(const ThresholdGraph& G) const {
    // Exactly the Nešetřil–Poljak / Eisenbrand–Grandoni triangle reduction
    // described in Section 3 of the paper:
    //
    // k1=floor(k/3), k2=ceil((k-1)/3), k3=ceil(k/3).
    // Vi contains all ki-cliques of G.
    // Two clique-nodes are adjacent iff their union is a clique.
    // G contains Kk iff the auxiliary 3-partite graph contains a triangle.
    const int k1 = k_ / 3;
    const int k2 = (k_ + 1) / 3; // ceil((k-1)/3)
    const int k3 = (k_ + 2) / 3; // ceil(k/3)

    const auto V1 = enumerate_cliques(G, k1);
    const auto V2 = enumerate_cliques(G, k2);
    const auto V3 = enumerate_cliques(G, k3);

    if (V1.empty() || V2.empty() || V3.empty()) return {};

    // M12[a,b] = 1 iff V1[a] and V2[b] are compatible.
    // M23[b,c] = 1 iff V2[b] and V3[c] are compatible.
    // M13[a,c] = 1 iff V1[a] and V3[c] are compatible.
    Matrix M12(V1.size(), V2.size(), 0);
    Matrix M23(V2.size(), V3.size(), 0);
    Matrix M13(V1.size(), V3.size(), 0);

    for (std::size_t a = 0; a < V1.size(); ++a) {
        for (std::size_t b = 0; b < V2.size(); ++b) {
            if (compatible(G, V1[a], V2[b])) M12(a, b) = 1;
        }
    }
    for (std::size_t b = 0; b < V2.size(); ++b) {
        for (std::size_t c = 0; c < V3.size(); ++c) {
            if (compatible(G, V2[b], V3[c])) M23(b, c) = 1;
        }
    }
    for (std::size_t a = 0; a < V1.size(); ++a) {
        for (std::size_t c = 0; c < V3.size(); ++c) {
            if (compatible(G, V1[a], V3[c])) M13(a, c) = 1;
        }
    }

    // FMM step: P = M12 M23. P[a,c] counts compatible 2-paths
    // V1[a] -> V2[b] -> V3[c]. A triangle exists exactly when
    // M13[a,c] = 1 and P[a,c] > 0.
    const Matrix P = multiply_blocked(M12, M23, gemm_block_);

    for (std::size_t a = 0; a < V1.size(); ++a) {
        for (std::size_t c = 0; c < V3.size(); ++c) {
            if (M13(a, c) == 0 || P(a, c) == 0) continue;

            // Recover one witness middle clique. This is outside the
            // asymptotic matrix multiplication and is only needed to return
            // the actual selected k points.
            for (std::size_t b = 0; b < V2.size(); ++b) {
                if (M12(a, b) && M23(b, c)) {
                    std::vector<int> witness = V1[a].v;
                    witness.insert(witness.end(), V2[b].v.begin(), V2[b].v.end());
                    witness.insert(witness.end(), V3[c].v.begin(), V3[c].v.end());
                    std::sort(witness.begin(), witness.end());
                    witness.erase(std::unique(witness.begin(), witness.end()), witness.end());
                    if ((int)witness.size() == k_) {
                        return {true, std::move(witness)};
                    }
                }
            }
        }
    }

    return {};
}

DispersionResult FMMKDispersion::solve() {
    build_distances();
    build_candidates();

    if (candidate_thresholds_.empty()) return {false, 0, {}};

    // The optimum minimum pairwise distance is one of the pairwise distances
    // in the input set. Feasibility is monotone as r decreases, so binary
    // search the sorted distinct candidate thresholds.
    int lo = 0;
    int hi = static_cast<int>(candidate_thresholds_.size()) - 1;
    int best = -1;
    std::vector<int> best_witness;

    while (lo <= hi) {
        const int mid = lo + (hi - lo) / 2;
        const auto r2 = candidate_thresholds_[mid];
        const ThresholdGraph G = build_graph(r2);
        const auto res = fmm_clique(G);

        if (res.yes) {
            best = mid;
            best_witness = res.witness;
            lo = mid + 1;
        } else {
            hi = mid - 1;
        }
    }

    if (best < 0) return {false, 0, {}};
    return {true, candidate_thresholds_[best], std::move(best_witness)};
}

} // namespace fmm
