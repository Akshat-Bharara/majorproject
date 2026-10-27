#include "dispersion.hpp"
#include <chrono>
#include <iostream>
#include <string>

using namespace fmm;

int main(int argc, char** argv) {
    // Optional flag --stats prints solver instrumentation to stderr.
    const bool want_stats = argc > 1 && std::string(argv[1]) == "--stats";

    // Input format:
    // n d k
    // x_1 ... x_d
    // ...
    // x_n ... x_d
    //
    // Coordinates are integral (signed 64-bit).
    // Output is the optimum squared minimum distance and one optimal set.

    int n, d, k;
    if (!(std::cin >> n >> d >> k)) {
        std::cerr << "usage: n d k followed by n integral d-dimensional points\n";
        return 1;
    }

    if (n <= 0 || d <= 0 || k < 2 || k > n) {
        std::cerr << "invalid n, d, or k\n";
        return 1;
    }

    std::vector<Point> points(n);
    for (int i = 0; i < n; ++i) {
        points[i].x.resize(d);
        for (int j = 0; j < d; ++j) {
            if (!(std::cin >> points[i].x[j])) {
                std::cerr << "invalid/incomplete point data\n";
                return 1;
            }
        }
    }

    try {
        FMMKDispersion solver(points, k);
        const auto t0 = std::chrono::steady_clock::now();
        const auto ans = solver.solve();
        const double solve_ms =
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();

        if (!ans.feasible) {
            std::cout << "NO SOLUTION\n";
            return 0;
        }

        std::cout << "optimal_squared_distance "
                  << uint128_to_string(ans.squared_distance) << '\n';
        std::cout << "selected_indices";
        for (int v : ans.selected) std::cout << ' ' << v;
        std::cout << '\n';

        std::cout << "selected_points\n";
        for (int v : ans.selected) {
            for (int j = 0; j < d; ++j) std::cout << (j ? " " : "") << points[v].x[j];
            std::cout << '\n';
        }

        if (want_stats) {
            const auto& st = ans.stats;
            const auto sp = clique_split(k);
            std::cerr << "split " << sp.k1 << ' ' << sp.k2 << ' ' << sp.k3 << '\n'
                      << "candidate_distances " << st.candidate_distances << '\n'
                      << "feasibility_tests " << st.feasibility_tests << '\n'
                      << "matrix_products " << st.matrix_products << '\n'
                      << "max_V " << st.max_V1 << ' ' << st.max_V2 << ' ' << st.max_V3 << '\n'
                      << "solve_ms " << solve_ms << '\n';
        }
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
