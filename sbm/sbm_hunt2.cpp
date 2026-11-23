// Adversarial search for inputs on which the two-way search tree of SBM
// (Claim 2 / Lemma 3) fails.  Simulated annealing on integer point sets:
// minimise (#roots from which the search reaches k points at r = r*, then
// #successful leaves), keeping the points in strictly convex position.
//   sbm_hunt <k> <n> <seed> <iterations> <restarts>
#include "sbm.hpp"
#include "../convex/randgen.hpp"
#include <cmath>
#include <iostream>
using namespace kd;

struct Count {
    const Instance& I; int k; i64 D2; std::vector<int> S; long leaves = 0;
    Count(const Instance& i, int kk, i64 d) : I(i), k(kk), D2(d) {}
    bool far_all(int v) const { for (int u : S) if (I.d(u, v) < D2) return false; return true; }
    void rec(int lo, int hi) {
        if ((int)S.size() == k) { ++leaves; return; }
        const int n = I.n, q = S[0];
        for (int side = 0; side < 2; ++side) {
            const bool ccw = side == 0;
            for (int s = 0; s < hi - lo - 1; ++s) {
                const int t = ccw ? lo + 1 + s : hi - 1 - s;
                const int v = (q + t) % n;
                if (!far_all(v)) continue;
                S.push_back(v);
                if (ccw) rec(t, hi); else rec(lo, t);
                S.pop_back();
                break;
            }
        }
    }
};
// score: (roots that succeed)*1000 + total successful leaves; -1 if invalid
static long score(const std::vector<Pt>& p, int k, i64* optout) {
    Instance I(p);
    try { I.validate(); } catch (...) { return -1; }
    if (I.n < k + 1) return -1;
    const i64 o = I.n <= 14 ? Brute(I, k).solve() : solve_capped(I, k).r2;
    if (optout) *optout = o;
    long roots = 0, leaves = 0;
    for (int v = 0; v < I.n; ++v) {
        Count c(I, k, o);
        c.S = {v};
        c.rec(0, I.n);
        roots += c.leaves > 0;
        leaves += c.leaves;
    }
    return roots * 1000 + leaves;
}
int main(int argc, char** argv) {
    const int k = atoi(argv[1]), n = atoi(argv[2]);
    std::mt19937_64 g(atoll(argv[3]));
    const long iters = atol(argv[4]);
    const int restarts = atoi(argv[5]);
    long best_overall = 1L << 60;
    for (int rs = 0; rs < restarts; ++rs) {
        std::vector<Pt> p;
        const double W = 1000, H = 20 + 400 * urand(g);
        if (argc > 6 && rs % 2 == 0) { // start from a given instance (file: x y per line)
            FILE* f = fopen(argv[6], "r"); long long x, y; p.clear();
            while (fscanf(f, "%lld %lld", &x, &y) == 2) p.push_back({x, y});
            fclose(f);
        }
        while ((int)p.size() != n) {
            p = (rs % 3 == 0) ? two_chain(g, n, W, H) : random_convex(g, rs % 3 == 1 ? F_ELLIPSE : F_CIRCLE, n);
            if (rs % 3 != 0) for (auto& q : p) { q.x /= 1000; q.y /= 1000; }
            p = strict_hull(p);
        }
        long cur = score(p, k, nullptr);
        if (cur < 0) continue;
        double T = 300;
        const double scale = 1 + std::max<double>(1, (double)[&] { kd::i64 m = 0; for (auto& q : p) m = std::max({m, std::llabs(q.x), std::llabs(q.y)}); return m; }() / 30.0);
        for (long it = 0; it < iters; ++it) {
            std::vector<Pt> q = p;
            const int i = g() % n;
            const double st = scale * std::pow(10.0, -2.0 * urand(g));
            q[i].x += (kd::i64)llround((2 * urand(g) - 1) * st);
            q[i].y += (kd::i64)llround((2 * urand(g) - 1) * st);
            q = strict_hull(q);
            if ((int)q.size() != n) continue;
            const long s = score(q, k, nullptr);
            if (s < 0) continue;
            if (s <= cur || urand(g) < std::exp((cur - s) / T)) { p = q; cur = s; }
            T *= 0.9995;
            if (cur < 1000) break; // no root succeeds
        }
        if (cur < best_overall) best_overall = cur;
        if (cur < 1000) {
            i64 o;
            score(p, k, &o);
            std::cout << "FOUND k=" << k << " n=" << n << " opt=" << o << " :";
            for (auto& q : p) std::cout << " " << q.x << " " << q.y;
            std::cout << std::endl;
            return 0;
        }
        if (rs % 10 == 9) std::cerr << "k=" << k << " n=" << n << " restart " << rs << " best score " << best_overall << "\n";
    }
    std::cout << "none; best score " << best_overall << "\n";
}
