// Checks conjecture V5 on many instances read from stdin.
//
// V5: if P (convex position) has k points pairwise >= r, then for some q the
// alternating greedy seeded with {q, a(q)} finds k such points. a(q) is the
// first point ccw from q at distance >= r; the ccw frontier starts at a(q),
// the cw frontier at q; sides alternate (both starting orders are tried) and
// each step adds the first point beyond the frontier, before the opposite
// frontier, that is far from every chosen point; a step that finds nothing
// fails the run.
//
// With --c1 it instead checks conjecture C1 (tight pair) exactly with the DP.
//
// Input: a sequence of instances "n 2 k" + n points (ccw convex order not
// required). For each instance prints "ok <tests>" or "FAIL r2 <r2>".

#include "../convex_dispersion.hpp"
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

using namespace cvx;

static bool v5_run(const std::vector<Pt>& p, int k, i64 r2, int q, int first) {
    const int n = (int)p.size();
    auto far = [&](int a, int b) { return dist2(p[a], p[b]) >= r2; };
    int a = -1;
    for (int s = 1; s < n && a < 0; ++s) if (far(q, (q + s) % n)) a = (q + s) % n;
    if (a < 0) return false;
    std::vector<int> chosen{q, a};
    int x = a, y = q, side = first;
    while ((int)chosen.size() < k) {
        const int step = side == 0 ? 1 : n - 1;
        const int from = side == 0 ? x : y, to = side == 0 ? y : x;
        int c = (from + step) % n, found = -1;
        for (; c != to; c = (c + step) % n) {
            bool ok = true;
            for (int v : chosen) if (!far(c, v)) { ok = false; break; }
            if (ok) { found = c; break; }
        }
        if (found < 0) return false;
        chosen.push_back(found);
        (side == 0 ? x : y) = found;
        side ^= 1;
    }
    return true;
}

int main(int argc, char** argv) {
    // --c1: test conjecture C1 (some feasible k-set contains a tight pair
    // q, a(q)) with the capped DP restricted to tight base chords.
    const bool c1 = argc > 1 && std::string(argv[1]) == "--c1";
    int n, d, k;
    long long total = 0, fails = 0, inst = 0;
    while (std::cin >> n >> d >> k) {
        std::vector<Pt> raw(n);
        for (auto& q : raw) std::cin >> q.x >> q.y;
        const auto h = convex_hull_ccw(raw);
        std::vector<Pt> p;
        for (int v : h) p.push_back(raw[v]);
        if ((int)p.size() < k || k < 3) continue;
        ++inst;
        ConvexDispersion solver(p, k);
        const i64 opt = solver.solve(Algo::CappedDP).squared_distance;
        std::vector<i64> cand;
        for (int a = 0; a < (int)p.size(); ++a)
            for (int b = a + 1; b < (int)p.size(); ++b) {
                const i64 v = dist2(p[a], p[b]);
                if (v <= opt) cand.push_back(v);
            }
        std::sort(cand.begin(), cand.end());
        cand.erase(std::unique(cand.begin(), cand.end()), cand.end());
        for (i64 r2 : cand) {
            ++total;
            bool ok = false;
            if (c1) ok = solver.decide_capped(r2, nullptr, true);
            for (int q = 0; q < (int)p.size() && !ok && !c1; ++q)
                ok = v5_run(p, k, r2, q, 0) || v5_run(p, k, r2, q, 1);
            if (!ok) {
                ++fails;
                std::cout << "FAIL n " << p.size() << " k " << k << " r2 " << r2
                          << " opt " << opt << "\n";
                for (const auto& q : p) std::cout << q.x << ' ' << q.y << '\n';
                std::cout.flush();
                break;
            }
        }
    }
    std::cout << "instances " << inst << " tests " << total << " failures " << fails << "\n";
}
