// k-dispersion for points in convex position.
//
// usage: convex_dispersion [--algo=capped|tw|certified|tight4] [--no-bounds] [--hull] [--stats]
// (tight4: the O(n^2 log n) algorithm for k = 4)
//
// Input (same format as the general solver, d must be 2):
//   n 2 k
//   x_1 y_1
//   ...
// The points must be in strictly convex position (every point a vertex of the
// convex hull) unless --hull is given, in which case the solver runs on the
// hull vertices only. Output indices refer to the input order.

#include "convex_dispersion.hpp"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <string>

using namespace cvx;

int main(int argc, char** argv) {
    Algo algo = Algo::CappedDP;
    bool stats = false, bounds = true, take_hull = false;
    for (int a = 1; a < argc; ++a) {
        const std::string s = argv[a];
        if (s == "--algo=capped") algo = Algo::CappedDP;
        else if (s == "--algo=tw") algo = Algo::TWNaive;
        else if (s == "--algo=certified") algo = Algo::Certified;
        else if (s == "--algo=tight4") algo = Algo::TightK4;
        else if (s == "--stats") stats = true;
        else if (s == "--no-bounds") bounds = false;
        else if (s == "--hull") take_hull = true;
        else { std::cerr << "unknown option " << s << '\n'; return 1; }
    }

    int n, d, k;
    if (!(std::cin >> n >> d >> k) || n <= 0 || d != 2 || k < 2) {
        std::cerr << "expected header: n 2 k (k >= 2)\n";
        return 1;
    }
    std::vector<Pt> pts(n);
    for (auto& p : pts)
        if (!(std::cin >> p.x >> p.y)) { std::cerr << "invalid point data\n"; return 1; }

    const std::vector<int> h = convex_hull_ccw(pts);
    if ((int)h.size() != n && !take_hull) {
        std::cerr << "points are not in strictly convex position (" << h.size()
                  << " of " << n << " are hull vertices); use --hull\n";
        return 1;
    }
    if (k > (int)h.size()) { std::cerr << "k exceeds number of hull vertices\n"; return 1; }
    std::vector<Pt> cp;
    for (int v : h) cp.push_back(pts[v]);

    try {
        ConvexDispersion solver(cp, k);
        const auto t0 = std::chrono::steady_clock::now();
        const ConvexResult res = solver.solve(algo, bounds);
        const double ms =
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();

        std::cout << "optimal_squared_distance " << res.squared_distance << '\n';
        std::cout << "selected_indices";
        std::vector<int> orig;
        for (int v : res.selected) orig.push_back(h[v]);
        std::sort(orig.begin(), orig.end());
        for (int v : orig) std::cout << ' ' << v;
        std::cout << '\n';
        if (stats) {
            const auto& s = res.stats;
            std::cerr << "hull_size " << h.size() << '\n'
                      << "candidates_total " << s.candidates_total << '\n'
                      << "candidates_window " << s.candidates_window << '\n'
                      << "decisions " << s.decisions << '\n'
                      << "greedy_decisions " << s.greedy_decisions << '\n'
                      << "certified_by_greedy " << (s.certified_by_greedy ? 1 : 0) << '\n'
                      << "inner_steps " << s.inner_steps << '\n'
                      << "solve_ms " << ms << '\n';
        }
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
