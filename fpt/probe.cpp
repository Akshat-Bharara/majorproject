// Exploratory probe of the two-sided search (SBM Claim 2).
#include "../convex/kdisp.hpp"
#include "../convex/randgen.hpp"
#include <iostream>
#include <map>
using namespace kd;
static int n, K; static i64 R2; static const Instance* IP;
static bool far(int a, int b) { return IP->d(a, b) >= R2; }
// search from state (C, lo, hi) offsets relative to root; returns success
static bool search(std::vector<int>& C, int lo, int hi) {
    if ((int)C.size() == K) return true;
    const int q = C[0];
    for (int side = 0; side < 2; ++side) {
        const bool ccw = side == 0;
        for (int s = 0; s < hi - lo - 1; ++s) {
            const int t = ccw ? lo + 1 + s : hi - 1 - s, v = (q + t) % n;
            bool ok = true; for (int u : C) ok = ok && far(u, v);
            if (!ok) continue;
            C.push_back(v);
            bool r = ccw ? search(C, t, hi) : search(C, lo, t);
            C.pop_back();
            if (r) return true;
            break;
        }
    }
    return false;
}
// good state: exists r-far k-set T ⊇ C with T\C in free arc (offsets lo+1..hi-1)
static std::vector<std::vector<int>> farsets;
static bool good(const std::vector<int>& C, int lo, int hi) {
    const int q = C[0];
    for (auto& T : farsets) {
        bool ok = true;
        for (int u : C) ok = ok && std::find(T.begin(), T.end(), u) != T.end();
        if (!ok) continue;
        for (int v : T) {
            if (std::find(C.begin(), C.end(), v) != C.end()) continue;
            const int t = ((v - q) % n + n) % n;
            if (!(t > lo && t < hi)) { ok = false; break; }
        }
        if (ok) return true;
    }
    return false;
}
static long dead = 0, goodstates = 0;
static void explore(std::vector<int>& C, int lo, int hi) { // walk all good states, count dead ones
    if ((int)C.size() == K) return;
    ++goodstates;
    const int q = C[0]; bool any = false;
    for (int side = 0; side < 2; ++side) {
        const bool ccw = side == 0;
        for (int s = 0; s < hi - lo - 1; ++s) {
            const int t = ccw ? lo + 1 + s : hi - 1 - s, v = (q + t) % n;
            bool ok = true; for (int u : C) ok = ok && far(u, v);
            if (!ok) continue;
            C.push_back(v);
            const int nlo = ccw ? t : lo, nhi = ccw ? hi : t;
            if (good(C, nlo, nhi)) { any = true; explore(C, nlo, nhi); }
            C.pop_back();
            break;
        }
    }
    if (!any) ++dead;
}
static void rec(std::vector<int>& S, int next) {
    if ((int)S.size() == K) { farsets.push_back(S); return; }
    for (int v = next; v < n; ++v) { bool ok = true; for (int u : S) ok = ok && far(u, v); if (ok) { S.push_back(v); rec(S, v + 1); S.pop_back(); } }
}
int main(int argc, char** argv) {
    std::mt19937_64 g(atoll(argv[2])); long cnt = atol(argv[1]);
    std::map<std::string, std::pair<long, long>> H;
    auto rec2 = [&](const std::string& h, bool ok) { auto& p = H[h]; ++p.first; p.second += !ok; };
    for (long it = 0; it < cnt; ++it) {
        int fam = it % 3 == 0 ? F_TWOCHAIN : (int)(g() % NUM_FAMILIES);
        Instance I(random_convex(g, fam, 6 + g() % 5)); if (I.n < 6) continue;
        IP = &I; n = I.n;
        for (i64 r2 : I.candidates()) {
            R2 = r2;
            for (K = 3; K <= std::min(n, 7); ++K) {
                farsets.clear(); std::vector<int> S; rec(S, 0);
                if (farsets.empty()) break;
                std::vector<char> succ(n);
                for (int v = 0; v < n; ++v) { std::vector<int> C{v}; succ[v] = search(C, 0, n); }
                bool anyT = false, allT = true, cp = true, cpBoth = true;
                for (auto& T : farsets) {
                    bool inT = false; for (int v : T) inT = inT || succ[v];
                    anyT = anyT || inT; allT = allT && inT;
                    int a = -1, b = -1; for (size_t x = 0; x < T.size(); ++x) for (size_t y = x + 1; y < T.size(); ++y) if (a < 0 || I.d(T[x], T[y]) < I.d(a, b)) { a = T[x]; b = T[y]; }
                    cp = cp && (succ[a] || succ[b]); cpBoth = cpBoth && succ[a] && succ[b];
                }
                rec2("some root succeeds (Claim 2), k=" + std::to_string(K), anyT);
                rec2("every far set T has a successful root in T", allT);
                rec2("every T: an endpoint of its closest pair succeeds", cp);
                rec2("every T: both endpoints of its closest pair succeed", cpBoth);
                int cntS = 0; for (int v = 0; v < n; ++v) cntS += succ[v];
                rec2("at least 2 successful roots", cntS >= 2);
                // exchange-free far sets: no w in a gap can replace an adjacent endpoint
                bool exfreeReach = true; long nfree = 0;
                for (auto& T : farsets) {
                    bool free_ = true;
                    for (size_t j = 0; j < T.size() && free_; ++j) {
                        const int a = T[j], b = T[(j + 1) % T.size()];
                        for (int w = (a + 1) % n; w != b && free_; w = (w + 1) % n) {
                            for (int e : {a, b}) {
                                bool ok = true; for (int u : T) if (u != e) ok = ok && far(u, w);
                                if (ok) free_ = false;
                            }
                        }
                    }
                    if (!free_) continue;
                    ++nfree;
                    bool any = false; for (int v : T) any = any || succ[v];
                    exfreeReach = exfreeReach && any;
                }
                rec2("exchange-free far sets exist", nfree > 0);
                rec2("every exchange-free T has a successful root in T", exfreeReach);
                // invariant: from root v in T, good states never dead?
                long d0 = dead;
                for (auto& T : farsets) for (int v : T) { std::vector<int> C{v}; explore(C, 0, n); }
                rec2("no dead good state", dead == d0);
            }
        }
    }
    for (auto& [h, p] : H) std::cout << h << ": " << p.first << " tests, " << p.second << " failures\n";
    std::cout << "good states " << goodstates << " dead " << dead << "\n";
}
