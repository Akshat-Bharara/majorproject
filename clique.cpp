#include "clique.hpp"
#include <algorithm>
#include <sstream>
#include <stdexcept>

namespace fmm {

static inline int word_index(int v) { return v >> 6; }
static inline std::uint64_t word_mask(int v) { return 1ULL << (v & 63); }

bool ThresholdGraph::edge(int u, int v) const noexcept {
    if (u == v) return false;
    return (adj[u][word_index(v)] & word_mask(v)) != 0;
}

static void dfs_cliques(const ThresholdGraph& G,
                        int q,
                        std::vector<int>& chosen,
                        std::vector<int>& candidates,
                        std::vector<Clique>& out) {
    if ((int)chosen.size() == q) {
        out.push_back({chosen});
        return;
    }

    const int need = q - static_cast<int>(chosen.size());
    if ((int)candidates.size() < need) return;

    // The increasing-index candidate discipline guarantees uniqueness.
    for (int pos = 0; pos < (int)candidates.size(); ++pos) {
        if ((int)candidates.size() - pos < need) break;

        const int u = candidates[pos];
        chosen.push_back(u);

        std::vector<int> next;
        next.reserve(candidates.size() - pos - 1);
        for (int t = pos + 1; t < (int)candidates.size(); ++t) {
            const int v = candidates[t];
            if (G.edge(u, v)) next.push_back(v);
        }

        dfs_cliques(G, q, chosen, next, out);
        chosen.pop_back();
    }
}

std::vector<Clique> enumerate_cliques(const ThresholdGraph& G, int q) {
    if (q < 0 || q > G.n) throw std::invalid_argument("invalid clique size");
    if (q == 0) return {Clique{}};

    std::vector<int> candidates(G.n);
    for (int i = 0; i < G.n; ++i) candidates[i] = i;

    std::vector<int> chosen;
    chosen.reserve(q);
    std::vector<Clique> out;
    dfs_cliques(G, q, chosen, candidates, out);
    return out;
}

bool compatible(const ThresholdGraph& G, const Clique& A, const Clique& B) {
    // The auxiliary-graph edge exists iff A union B induces a clique.
    // This also enforces A and B to be vertex-disjoint.
    for (int u : A.v) {
        for (int v : B.v) {
            if (!G.edge(u, v)) return false;
        }
    }
    return true;
}

std::string clique_to_string(const Clique& c) {
    std::ostringstream oss;
    oss << "{";
    for (std::size_t i = 0; i < c.v.size(); ++i) {
        if (i) oss << ",";
        oss << c.v[i];
    }
    oss << "}";
    return oss.str();
}

} // namespace fmm
