// The k = 3 approximation of SBM (Section 5, Theorem 7: ratio 1/(2 sqrt 2)).
// Points: extreme points a (min x), b (max y), c (max x), d (min y).
//  Case 1: for every pair u != v of extremes, e = point farthest from line uv: triple (u, v, e)
//  Case 2: for every pair, f = point nearest to the perpendicular bisector of uv: triple (u, v, f)
//  Case 3: best triple of distinct extreme points
// Returns the best triple (squared min distance).  Farthest/nearest are found by a
// linear scan (same result as the O(log n) search of the paper); exact integers.
//   sbm_approx random <instances> <seed>     ratio statistics vs. exact optimum
//   sbm_approx anneal <n> <seed> <iters> <restarts>   minimise the ratio
#include "../convex/kdisp.hpp"
#include "../convex/randgen.hpp"
#include <cmath>
#include <cstring>
#include <iostream>
using namespace kd;

static i64 tri(const Instance& I, int a, int b, int c) {
    if (a == b || b == c || a == c) return -1;
    return std::min({I.d(a, b), I.d(b, c), I.d(a, c)});
}
i64 sbm_k3(const Instance& I) {
    const int n = I.n;
    int ex[4] = {0, 0, 0, 0}; // a: min x, b: max y, c: max x, d: min y (ties: first)
    for (int i = 0; i < n; ++i) {
        const Pt& p = I.p[i];
        if (p.x < I.p[ex[0]].x) ex[0] = i;
        if (p.y > I.p[ex[1]].y) ex[1] = i;
        if (p.x > I.p[ex[2]].x) ex[2] = i;
        if (p.y < I.p[ex[3]].y) ex[3] = i;
    }
    i64 best = -1;
    for (int s = 0; s < 4; ++s)
        for (int t = s + 1; t < 4; ++t) {
            const int u = ex[s], v = ex[t];
            if (u == v) continue;
            const Pt &U = I.p[u], &V = I.p[v];
            int e = -1, f = -1;
            i128 be = -1, bf = -1;
            for (int w = 0; w < n; ++w) {
                if (w == u || w == v) continue;
                const Pt& W = I.p[w];
                i128 cr = orient(U, V, W); if (cr < 0) cr = -cr;              // distance to line uv (scaled)
                i128 dt = (i128)(2 * W.x - U.x - V.x) * (V.x - U.x) + (i128)(2 * W.y - U.y - V.y) * (V.y - U.y);
                if (dt < 0) dt = -dt;                                          // distance to bisector (scaled)
                if (cr > be) { be = cr; e = w; }
                if (bf < 0 || dt < bf) { bf = dt; f = w; }
            }
            if (e >= 0) best = std::max(best, tri(I, u, v, e));
            if (f >= 0) best = std::max(best, tri(I, u, v, f));
        }
    for (int s = 0; s < 4; ++s) for (int t = s + 1; t < 4; ++t) for (int w = t + 1; w < 4; ++w) best = std::max(best, tri(I, ex[s], ex[t], ex[w]));
    return best;
}
static double ratio(const Instance& I) {
    Brute B(I, 3);
    const i64 o = B.solve();
    return std::sqrt((double)sbm_k3(I) / (double)o);
}
int main(int argc, char** argv) {
    std::mt19937_64 g(atoll(argv[3]));
    const double bound = 1 / (2 * std::sqrt(2.0));
    if (!strcmp(argv[1], "random")) {
        const long cnt = atol(argv[2]);
        double worst = 1; long below = 0, tot = 0;
        for (long it = 0; it < cnt; ++it) {
            Instance I(random_convex(g, it % NUM_FAMILIES, 4 + g() % 30));
            if (I.n < 4) continue;
            const double r = ratio(I);
            ++tot; worst = std::min(worst, r); below += r < bound - 1e-12;
        }
        std::cout << "instances " << tot << " worst ratio " << worst << " below 1/(2sqrt2): " << below << "\n";
        return 0;
    }
    const int n = atoi(argv[2]); const long iters = atol(argv[4]); const int restarts = atoi(argv[5]);
    double best = 1; std::vector<Pt> bp;
    for (int rs = 0; rs < restarts; ++rs) {
        std::vector<Pt> p;
        while ((int)p.size() != n) p = random_convex(g, rs % NUM_FAMILIES, n);
        for (auto& q : p) { q.x /= 1000; q.y /= 1000; }
        p = strict_hull(p);
        if ((int)p.size() != n) continue;
        double cur = ratio(Instance(p)), T = 0.02;
        for (long it = 0; it < iters; ++it) {
            std::vector<Pt> q = p;
            const int i = g() % n;
            const double st = 1000 * std::pow(10.0, -3 * urand(g));
            q[i].x += (i64)llround((2 * urand(g) - 1) * st);
            q[i].y += (i64)llround((2 * urand(g) - 1) * st);
            q = strict_hull(q);
            if ((int)q.size() != n) continue;
            const double r = ratio(Instance(q));
            if (r <= cur || urand(g) < std::exp((cur - r) / T)) { p = q; cur = r; }
            if (cur < best) { best = cur; bp = p; }
            T *= 0.9997;
        }
    }
    std::cout << "n=" << n << " worst ratio " << best << " (bound " << bound << ") :";
    for (auto& q : bp) std::cout << " " << q.x << " " << q.y;
    std::cout << "\n";
}
