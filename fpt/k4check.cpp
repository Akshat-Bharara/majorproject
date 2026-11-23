// Check of the proof of SBM_4: for every r-far 4-set T, some adjacent pair (u,v) with u a nearest
// neighbour of v in T gives a witness-preserving search: root u, first step towards v, then a safe
// second step (either side), then the last step.  In case R (two short opposite edges) the proof
// says one of the roots a (towards b) and c (towards d) works.
//   k4check random <instances> <seed>   |   k4check anneal <n> <seed> <iters> <restarts>
#include "../convex/kdisp.hpp"
#include "../convex/randgen.hpp"
#include <cmath>
#include <cstring>
#include <iostream>
using namespace kd;
static const Instance* IP; static int n; static i64 R2;
static bool far(int a, int b) { return IP->d(a, b) >= R2; }
static int off(int q, int v) { return ((v - q) % n + n) % n; }
// greedy step: C chosen with root q, free offsets (lo,hi); returns point or -1
static int step(const std::vector<int>& C, int q, int lo, int hi, bool ccw, int& t) {
    if (ccw) { for (t = lo + 1; t < hi; ++t) { int v = (q + t) % n; bool ok = true; for (int c : C) ok = ok && far(c, v); if (ok) return v; } }
    else { for (t = hi - 1; t > lo; --t) { int v = (q + t) % n; bool ok = true; for (int c : C) ok = ok && far(c, v); if (ok) return v; } }
    return -1;
}
// T = {T[0..3]} ccw; root T[i], first step towards T[i+1] (dir=+1) or T[i-1] (dir=-1); returns success via safe steps
static bool attempt(const std::vector<int>& T, int i, int dir) {
    const int u = T[i], q = u;
    std::vector<int> rem; for (int j = 1; j < 4; ++j) rem.push_back(T[(i + j) % 4]); // ccw order from u
    std::vector<int> C{u}; int lo = 0, hi = n, t;
    // step 1
    const bool ccw1 = dir == 1;
    int x = step(C, q, lo, hi, ccw1, t); if (x < 0) return false;
    std::vector<int> r1 = rem; if (ccw1) r1.erase(r1.begin()); else r1.pop_back();
    for (int s : r1) if (!far(s, x)) return false;            // first step must be safe
    C.push_back(x); if (ccw1) lo = t; else hi = t;
    // step 2: try both sides
    for (int side = 0; side < 2; ++side) {
        const bool ccw = side == 0; int t2;
        int y = step(C, q, lo, hi, ccw, t2); if (y < 0) continue;
        std::vector<int> r2 = r1; if (ccw) r2.erase(r2.begin()); else r2.pop_back();
        bool safe = true; for (int s : r2) safe = safe && far(s, y);
        if (!safe) continue;
        std::vector<int> C2 = C; C2.push_back(y); int lo2 = ccw ? t2 : lo, hi2 = ccw ? hi : t2, t3;
        if (step(C2, q, lo2, hi2, true, t3) >= 0) return true;
    }
    return false;
}
long tests = 0, fails = 0, rtests = 0, rfails = 0;
static void check(const Instance& I) {
    IP = &I; n = I.n;
    for (i64 r2 : I.candidates()) {
        R2 = r2;
        for (int a = 0; a < n; ++a) for (int b = a + 1; b < n; ++b) if (far(a, b)) for (int c = b + 1; c < n; ++c) if (far(a, c) && far(b, c))
            for (int d = c + 1; d < n; ++d) if (far(a, d) && far(b, d) && far(c, d)) {
                std::vector<int> T{a, b, c, d}; ++tests;
                bool ok = false;
                for (int i = 0; i < 4 && !ok; ++i) for (int dir = -1; dir <= 1 && !ok; dir += 2) {
                    const int v = T[(i + dir + 4) % 4], u = T[i];
                    i64 nn = INT64_MAX; for (int w : T) if (w != v) nn = std::min(nn, I.d(v, w));
                    if (I.d(u, v) == nn) ok = attempt(T, i, dir);
                }
                if (!ok) { ++fails; std::cout << "FAIL r2=" << r2 << " T " << a << b << c << d << " pts:"; for (auto& p : I.p) std::cout << " " << p.x << " " << p.y; std::cout << std::endl; }
                // case R check: rotations where T[i]T[i+1] and T[i+2]T[i+3] are both strictly shorter than the other four pairs
                for (int i = 0; i < 2; ++i) {
                    int A = T[i], B = T[i + 1], Cc = T[(i + 2) % 4], D = T[(i + 3) % 4];
                    i64 s1 = I.d(A, B), s2 = I.d(Cc, D), m = std::min({I.d(A, Cc), I.d(A, D), I.d(B, Cc), I.d(B, D)});
                    if (std::max(s1, s2) < m) { ++rtests; if (!attempt(T, i, 1) && !attempt(T, (i + 2) % 4, 1)) ++rfails; }
                }
            }
    }
}
int main(int argc, char** argv) {
    std::mt19937_64 g(atoll(argv[3]));
    if (!strcmp(argv[1], "random")) {
        for (long it = 0; it < atol(argv[2]); ++it) { Instance I(random_convex(g, it % NUM_FAMILIES, 4 + g() % 11)); if (I.n >= 4) check(I); }
    } else {
        const int nn = atoi(argv[2]); const long iters = atol(argv[4]); const int rs = atoi(argv[5]);
        for (int r = 0; r < rs; ++r) {
            std::vector<Pt> p; while ((int)p.size() != nn) p = strict_hull(two_chain(g, nn, 1000, 20 + 500 * urand(g)));
            for (long it = 0; it < iters; ++it) {
                std::vector<Pt> q = p; int i = g() % nn; double st = 50 * std::pow(10.0, -2 * urand(g));
                q[i].x += llround((2 * urand(g) - 1) * st); q[i].y += llround((2 * urand(g) - 1) * st);
                q = strict_hull(q); if ((int)q.size() != nn) continue; p = q;
                if (it % 10 == 0) check(Instance(p));
            }
        }
    }
    std::cout << "4-set tests " << tests << " failures " << fails << " | case R tests " << rtests << " failures " << rfails << std::endl;
}
