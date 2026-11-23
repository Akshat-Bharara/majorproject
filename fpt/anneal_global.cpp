// Adversarial search: make a root u of an r-far k-set T fail the direct-exchange
// strategy while u's angle in T is as SMALL as possible.  Reports the smallest
// failing-root angle found and any violation of
//   Lemma S: if u fails, some other vertex of T has a smaller angle (=> SA rule).
//   anneal_stuck <k> <n> <seed> <iters> <restarts>
#include "../convex/kdisp.hpp"
#include "../convex/randgen.hpp"
#include <cmath>
#include <iostream>
using namespace kd;
static int n, K; static i64 R2; static const Instance* IP;
static bool far(int a, int b) { return IP->d(a, b) >= R2; }
static std::vector<std::vector<int>> fs;
static void rec(std::vector<int>& S, int next) {
    if ((int)S.size() == K) { fs.push_back(S); return; }
    for (int v = next; v < n; ++v) { bool ok = true; for (int u : S) ok = ok && far(u, v); if (ok) { S.push_back(v); rec(S, v + 1); S.pop_back(); } }
}
static int offq(int q, int v) { return ((v - q) % n + n) % n; }
static bool dfsB(std::vector<int>& C, std::vector<int>& rem, int root, int lo, int hi) {
    if ((int)C.size() == K) return true;
    for (int side = 0; side < 2; ++side) {
        const bool ccw = side == 0; int x = -1, t = -1;
        if (ccw) { for (t = lo + 1; t < hi; ++t) { int v = (root + t) % n; bool ok = true; for (int c : C) ok = ok && far(c, v); if (ok) { x = v; break; } } }
        else { for (t = hi - 1; t > lo; --t) { int v = (root + t) % n; bool ok = true; for (int c : C) ok = ok && far(c, v); if (ok) { x = v; break; } } }
        if (x < 0) continue;
        std::vector<int> r2 = rem; if (ccw) r2.erase(r2.begin()); else r2.pop_back();
        bool valid = true; for (int s2 : r2) valid = valid && far(s2, x);
        if (!valid) continue;
        C.push_back(x); const bool ok = dfsB(C, r2, root, ccw ? t : lo, ccw ? hi : t); C.pop_back();
        if (ok) return true;
    }
    return false;
}
static double angdeg(const Instance& I, int a, int p, int b) {
    const Pt &P = I.p[p], &A = I.p[a], &B = I.p[b];
    const long double ax = A.x - P.x, ay = A.y - P.y, bx = B.x - P.x, by = B.y - P.y;
    return (double)(acosl((ax * bx + ay * by) / sqrtl((ax * ax + ay * ay) * (bx * bx + by * by))) * 180 / 3.14159265358979L);
}
long violS = 0, checksG = 0;
// Global rule: over ALL r-far k-sets T and roots u in T take the pair with the smallest angle theta_u(T).
// Score for annealing: (#roots over all T that succeed) -- smaller = harder; report violations of the global rule.
static double eval(const Instance& I, i64 r2) {
    IP = &I; n = I.n; R2 = r2; fs.clear(); std::vector<int> S; rec(S, 0);
    if (fs.empty()) return 999;
    double gmin = 1e9; int gT = -1, gi = -1; long okRoots = 0, total = 0;
    std::vector<std::vector<char>> succ(fs.size());
    for (size_t ti = 0; ti < fs.size(); ++ti) {
        auto& T = fs[ti]; const int k = T.size(); succ[ti].assign(k, 0);
        for (int i = 0; i < k; ++i) {
            const double th = angdeg(I, T[(i + k - 1) % k], T[i], T[(i + 1) % k]);
            std::vector<int> C{T[i]}, rem; for (int v : T) if (v != T[i]) rem.push_back(v);
            std::sort(rem.begin(), rem.end(), [&](int a, int b) { return offq(T[i], a) < offq(T[i], b); });
            succ[ti][i] = dfsB(C, rem, T[i], 0, n); okRoots += succ[ti][i]; ++total;
            if (th < gmin) { gmin = th; gT = ti; gi = i; }
        }
    }
    ++checksG;
    if (!succ[gT][gi]) {
        ++violS;
        if (violS <= 3) { std::cout << "GLOBAL SMALLEST-ANGLE ROOT FAILS k=" << K << " r2=" << r2 << " pts:"; for (auto& p : I.p) std::cout << " " << p.x << " " << p.y; std::cout << std::endl; }
    }
    return (double)okRoots / total * 100; // percentage of successful (T, root) pairs
}
int main(int argc, char** argv) {
    K = atoi(argv[1]); const int nn = atoi(argv[2]); std::mt19937_64 g(atoll(argv[3]));
    const long iters = atol(argv[4]); const int restarts = atoi(argv[5]);
    double gbest = 999; std::vector<Pt> gp;
    for (int rs = 0; rs < restarts; ++rs) {
        std::vector<Pt> p;
        while ((int)p.size() != nn) { p = two_chain(g, nn, 1000, 30 + 500 * urand(g)); p = strict_hull(p); }
        auto score = [&](const std::vector<Pt>& q) -> double {
            Instance I(q); try { I.validate(); } catch (...) { return -1; }
            const auto Cn = I.candidates(); double b = 999;
            for (int t = 0; t < 4; ++t) b = std::min(b, eval(I, Cn[g() % Cn.size()]));
            return b;
        };
        double cur = score(p); if (cur < 0) continue;
        double T = 20;
        for (long it = 0; it < iters; ++it) {
            std::vector<Pt> q = p; const int i = g() % nn;
            const double st = 60 * std::pow(10.0, -2 * urand(g));
            q[i].x += (i64)llround((2 * urand(g) - 1) * st); q[i].y += (i64)llround((2 * urand(g) - 1) * st);
            q = strict_hull(q); if ((int)q.size() != nn) continue;
            const double s = score(q); if (s < 0) continue;
            if (s <= cur || urand(g) < std::exp((cur - s) / T)) { p = q; cur = s; }
            if (cur < gbest) { gbest = cur; gp = p; }
            T *= 0.999;
        }
    }
    std::cout << "k=" << K << " n=" << nn << " min % successful roots " << gbest << " deg; smallest-angle-root failures " << violS << " :";
    for (auto& q : gp) std::cout << " " << q.x << " " << q.y; std::cout << std::endl;
}
