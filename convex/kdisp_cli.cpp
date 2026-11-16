// Command-line solver.
//
//   kdisp [--algo capped|certified|triple|brute] [--mis R2] < instance
//
// Input:  "n k" followed by n lines "x y": integer points in strictly convex
//         position in counterclockwise order, |x|,|y| <= 2^29.
// Output: the optimal squared distance and one optimal k-set (0-based indices).
// With --mis R2: a maximum subset with pairwise squared distances >= R2.
#include "kdisp.hpp"
#include <cstring>
#include <iostream>

using namespace kd;

int main(int argc, char** argv) {
    std::string algo = "certified";
    i64 mis = -1;
    for (int a = 1; a < argc; ++a) {
        if (!strcmp(argv[a], "--algo") && a + 1 < argc) algo = argv[++a];
        else if (!strcmp(argv[a], "--mis") && a + 1 < argc) mis = atoll(argv[++a]);
        else { std::cerr << "usage: kdisp [--algo capped|certified|triple|brute] [--mis R2] < instance\n"; return 2; }
    }
    Instance I;
    int k;
    if (!read_instance(std::cin, I, k)) { std::cerr << "error: cannot parse instance\n"; return 1; }
    try {
        I.validate();
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    if (mis >= 0) {
        const std::vector<int> M = max_far_set(I, mis);
        std::cout << "max_far_set_size " << M.size() << "\nselected_indices";
        for (int v : M) std::cout << ' ' << v;
        std::cout << "\n";
        return 0;
    }
    if (k < 2 || k > I.n) { std::cerr << "error: need 2 <= k <= n\n"; return 1; }
    Result R;
    if (algo == "capped") R = solve_capped(I, k);
    else if (algo == "certified") R = solve_certified(I, k);
    else if (algo == "triple") { R = solve_triple(I, k); CappedDP w(I); if (k >= 3) { w.run(R.r2, k); R.S = w.witness(); } }
    else if (algo == "brute") { Brute B(I, k); R.r2 = B.solve(); R.S = B.arg; }
    else { std::cerr << "error: unknown algorithm " << algo << "\n"; return 2; }
    std::cout << "optimal_squared_distance " << R.r2 << "\nselected_indices";
    for (int v : R.S) std::cout << ' ' << v;
    std::cout << "\nexact_decisions " << R.decisions << "\n";
    if (algo == "certified") std::cout << "certified_by_greedy " << (R.certified_by_greedy ? "yes" : "no") << "\n";
    return 0;
}
