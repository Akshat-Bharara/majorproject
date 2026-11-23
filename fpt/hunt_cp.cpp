// Stress test of hypothesis CP: for every r-far k-set T and every endpoint u
// of a closest pair of T, the two-sided search rooted at u reaches k points.
// Annealing minimises the number of successful roots (at r = r* and at random
// feasible radii); every evaluated instance is checked for CP.
//   hunt_cp <k> <n> <seed> <iters> <restarts>
#include "../convex/kdisp.hpp"
#include "../convex/randgen.hpp"
#include <cmath>
#include <iostream>
using namespace kd;
static int n, K; static i64 R2; static const Instance* IP;
static bool far(int a, int b) { return IP->d(a, b) >= R2; }
static bool search(std::vector<int>& C, int lo, int hi) {
    if ((int)C.size() == K) return true;
    const int q = C[0];
    for (int side = 0; side < 2; ++side) {
        const bool ccw = side == 0;
        for (int s = 0; s < hi - lo - 1; ++s) {
            const int t = ccw ? lo + 1 + s : hi - 1 - s, v = (q + t) % n;
            bool ok = true; for (int u : C) ok = ok && far(u, v);
            if (!ok) continue;
            C.push_back(v); bool r = ccw ? search(C, t, hi) : search(C, lo, t); C.pop_back();
            if (r) return true;
            break;
        }
    }
    return false;
}
static std::vector<std::vector<int>> fs;
static void rec(std::vector<int>& S, int next) {
    if ((int)S.size() == K) { fs.push_back(S); return; }
    for (int v = next; v < n; ++v) { bool ok = true; for (int u : S) ok = ok && far(u, v); if (ok) { S.push_back(v); rec(S, v + 1); S.pop_back(); } }
}
long checks = 0, viol = 0;
// returns #successful roots at radius r2; checks CP
static int analyse(const Instance& I, i64 r2) {
    IP = &I; n = I.n; R2 = r2; fs.clear(); std::vector<int> S; rec(S, 0);
    if (fs.empty()) return -1;
    std::vector<char> succ(n);
    int cnt = 0;
    for (int v = 0; v < n; ++v) { std::vector<int> C{v}; succ[v] = search(C, 0, n); cnt += succ[v]; }
    for (auto& T : fs) {
        i64 dmin = INT64_MAX;
        for (size_t a = 0; a < T.size(); ++a) for (size_t b = a + 1; b < T.size(); ++b) dmin = std::min(dmin, I.d(T[a], T[b]));
        for (size_t a = 0; a < T.size(); ++a) for (size_t b = a + 1; b < T.size(); ++b)
            if (I.d(T[a], T[b]) == dmin) {
                ++checks;
                if (!succ[T[a]] || !succ[T[b]]) {
                    ++viol;
                    std::cout << "CP VIOLATION k=" << K << " r2=" << r2 << " T:"; for (int v : T) std::cout << " " << v;
                    std::cout << " pts:"; for (auto& p : I.p) std::cout << " " << p.x << " " << p.y; std::cout << std::endl;
                }
            }
    }
    return cnt;
}
int main(int argc, char** argv) {
    const int k = atoi(argv[1]), nn = atoi(argv[2]); std::mt19937_64 g(atoll(argv[3]));
    const long iters = atol(argv[4]); const int restarts = atoi(argv[5]);
    K = k; int best = 1 << 30;
    for (int rs = 0; rs < restarts && viol == 0; ++rs) {
        std::vector<Pt> p;
        while ((int)p.size() != nn) { p = (rs % 2) ? two_chain(g, nn, 1000, 20 + 400 * urand(g)) : random_convex(g, rs % 4 == 0 ? F_ELLIPSE : F_DISK, nn); if (rs % 2 == 0) for (auto& q : p) { q.x /= 1000; q.y /= 1000; } p = strict_hull(p); }
        auto score = [&](const std::vector<Pt>& q) -> int {
            Instance I(q); try { I.validate(); } catch (...) { return -1; }
            Brute B(I, k); const i64 o = B.solve(); if (o <= 0) return -1;
            int s = analyse(I, o);
            const auto C = I.candidates(); // also a random feasible radius
            const i64 r = C[g() % C.size()]; if (r <= o) analyse(I, r);
            return s;
        };
        int cur = score(p); if (cur < 0) continue;
        double T = 1.0;
        for (long it = 0; it < iters && viol == 0; ++it) {
            std::vector<Pt> q = p; const int i = g() % nn;
            const double st = 40 * std::pow(10.0, -2 * urand(g));
            q[i].x += (i64)llround((2 * urand(g) - 1) * st); q[i].y += (i64)llround((2 * urand(g) - 1) * st);
            q = strict_hull(q); if ((int)q.size() != nn) continue;
            const int s = score(q); if (s < 0) continue;
            if (s <= cur || urand(g) < std::exp((cur - s) / T)) { p = q; cur = s; }
            T *= 0.9995; best = std::min(best, cur);
        }
    }
    std::cout << "k=" << k << " n=" << nn << " CP checks " << checks << " violations " << viol << " min successful roots " << best << std::endl;
}
