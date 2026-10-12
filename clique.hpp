#pragma once

#include <cstdint>
#include <vector>
#include <string>

namespace fmm {

struct ThresholdGraph {
    int n = 0;
    std::vector<std::vector<std::uint64_t>> adj; // symmetric bitset adjacency

    bool edge(int u, int v) const noexcept;
};

struct Clique {
    std::vector<int> v; // sorted vertex IDs
};

// Enumerate every q-clique of G exactly once.
std::vector<Clique> enumerate_cliques(const ThresholdGraph& G, int q);

bool compatible(const ThresholdGraph& G, const Clique& A, const Clique& B);

std::string clique_to_string(const Clique& c);

} // namespace fmm
