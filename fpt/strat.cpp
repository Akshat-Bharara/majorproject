// Test witness-guided strategies for the two-sided search (towards a proof of Claim 2).
// State: chosen C (interval around root), witness W = C ∪ Rem (r-far, |W| = k),
// Rem = original points of T not yet replaced, in the free arc.
// Proven: the ccw greedy step x keeps (W - t1 + x) r-far whenever
//   min_{c in C} |c t1| <= min_{t in Rem, t != t1} |t1 t|   ("t1 is attached").
// Strategy A: take an attached side (ccw first); fail if neither side is attached.
#include "../convex/kdisp.hpp"
#include "../convex/randgen.hpp"
#include <iostream>
#include <map>
using namespace kd;
static int n, K; static i64 R2; static const Instance* IP;
static bool far(int a, int b) { return IP->d(a, b) >= R2; }
static std::vector<std::vector<int>> fs;
static void rec(std::vector<int>& S, int next) {
    if ((int)S.size() == K) { fs.push_back(S); return; }
    for (int v = next; v < n; ++v) { bool ok = true; for (int u : S) ok = ok && far(u, v); if (ok) { S.push_back(v); rec(S, v + 1); S.pop_back(); } }
}
static int off(int q, int v) { return ((v - q) % n + n) % n; }
// Strategy B: only direct exchanges with the current witness W.
// mode 0: prefer ccw, 1: prefer cw, 2: DFS over both valid exchanges.
static bool dfs(std::vector<int>& C, std::vector<int>& rem, int root, int lo, int hi, int mode) {
    if ((int)C.size() == K) return true;
    for (int side = 0; side < 2; ++side) {
        const bool ccw = (mode == 1) ? side == 1 : side == 0;
        int x = -1, t = -1;
        if (ccw) { for (t = lo + 1; t < hi; ++t) { int v = (root + t) % n; bool ok = true; for (int c : C) ok = ok && far(c, v); if (ok) { x = v; break; } } }
        else { for (t = hi - 1; t > lo; --t) { int v = (root + t) % n; bool ok = true; for (int c : C) ok = ok && far(c, v); if (ok) { x = v; break; } } }
        if (x < 0) continue;
        std::vector<int> r2 = rem;
        if (ccw) r2.erase(r2.begin()); else r2.pop_back();
        bool valid = true; for (int s2 : r2) valid = valid && far(s2, x);
        if (!valid) continue;
        C.push_back(x);
        const bool ok = dfs(C, r2, root, ccw ? t : lo, ccw ? hi : t, mode);
        C.pop_back();
        if (ok) return true;
        if (mode != 2) return false;
    }
    return false;
}
static int run(const std::vector<int>& T, int root, int mode) {
    std::vector<int> C{root};
    std::vector<int> rem; for (int v : T) if (v != root) rem.push_back(v);
    std::sort(rem.begin(), rem.end(), [&](int a, int b) { return off(root, a) < off(root, b); });
    return dfs(C, rem, root, 0, n, mode) ? 1 : 0;
}
int main(int argc, char** argv) {
    std::mt19937_64 g(atoll(argv[2])); const long cnt = atol(argv[1]);
    std::map<std::string, std::pair<long, long>> H;
    auto rc = [&](const std::string& h, bool ok) { auto& p = H[h]; ++p.first; p.second += !ok; };
    long bugs = 0;
    for (long it = 0; it < cnt; ++it) {
        const int fam = it % 3 == 0 ? F_TWOCHAIN : (int)(g() % NUM_FAMILIES);
        Instance I(random_convex(g, fam, 6 + g() % 7)); if (I.n < 6) continue;
        IP = &I; n = I.n;
        for (i64 r2 : I.candidates()) {
            R2 = r2;
            for (K = 3; K <= std::min(n, 7); ++K) {
                fs.clear(); std::vector<int> S; rec(S, 0); if (fs.empty()) break;
                for (auto& T : fs) {
                    i64 dmin = INT64_MAX; int cu = -1, cw = -1;
                    for (size_t a = 0; a < T.size(); ++a) for (size_t b = a + 1; b < T.size(); ++b) if (I.d(T[a], T[b]) < dmin) { dmin = I.d(T[a], T[b]); cu = T[a]; cw = T[b]; }
                    for (int rule = 0; rule < 3; ++rule) {
                        const int r1 = run(T, cu, rule), r2b = run(T, cw, rule);
                        bugs += (r1 < 0) + (r2b < 0);
                        rc("strategy " + std::to_string(rule) + " from closest-pair endpoint", r1 == 1 && r2b == 1);
                        bool anyRoot = false; for (int v : T) anyRoot = anyRoot || run(T, v, rule) == 1;
                        rc("strategy " + std::to_string(rule) + " from some root of T", anyRoot);
                    }
                }
            }
        }
    }
    for (auto& [h, p] : H) std::cout << h << ": " << p.first << " tests, " << p.second << " failures\n";
    std::cout << "invariant bugs " << bugs << "\n";
}
