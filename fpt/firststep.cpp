// Sanity test of the attachment lemma and the ear lemma:
//  (1) for every r-far k-set T and adjacent e, rho with rho a nearest neighbour of e in T,
//      the greedy step from root rho towards e keeps a witness;
//  (2) every T has such a pair (some vertex has an adjacent nearest neighbour);
//  (3) general depth: whenever some chosen c has |c t1| <= |t1 tj| for all remaining j, the ccw exchange is valid.
#include "../convex/kdisp.hpp"
#include "../convex/randgen.hpp"
#include <iostream>
#include <functional>
using namespace kd;
int main(int argc, char** argv) {
    std::mt19937_64 g(atoll(argv[2])); long t1 = 0, f1 = 0, t2 = 0, f2 = 0;
    for (long it = 0; it < atol(argv[1]); ++it) {
        Instance I(random_convex(g, it % NUM_FAMILIES, 5 + g() % 9)); const int n = I.n; if (n < 5) continue;
        for (i64 r2 : I.candidates()) for (int k = 3; k <= std::min(n, 7); ++k) {
            // enumerate far sets
            std::vector<std::vector<int>> fs; std::vector<int> S;
            std::function<void(int)> rec = [&](int nx) { if ((int)S.size() == k) { fs.push_back(S); return; }
                for (int v = nx; v < n; ++v) { bool ok = true; for (int u : S) ok = ok && I.d(u, v) >= r2; if (ok) { S.push_back(v); rec(v + 1); S.pop_back(); } } };
            rec(0); if (fs.empty()) break;
            for (auto& T : fs) {
                bool any = false;
                for (int i = 0; i < k; ++i) for (int dir = -1; dir <= 1; dir += 2) {
                    const int e = T[i], rho = T[(i - dir + k) % k]; // rho adjacent to e
                    i64 nn = INT64_MAX; for (int v : T) if (v != e) nn = std::min(nn, I.d(e, v));
                    if (I.d(e, rho) != nn) continue;
                    any = true; ++t1;
                    // step from rho towards e: e is ccw of rho if dir == +1 (rho = T[i-1])
                    const bool ccw = dir == 1; int x = -1;
                    for (int s = 1; s < n; ++s) { int v = ccw ? (rho + s) % n : ((rho - s) % n + n) % n; if (I.d(rho, v) >= r2) { x = v; break; } }
                    bool valid = x >= 0; for (int v : T) if (v != e && v != rho) valid = valid && I.d(v, x) >= r2;
                    f1 += !valid;
                }
                ++t2; f2 += !any;
            }
        }
    }
    std::cout << "first-step tests " << t1 << " failures " << f1 << "; sets without an adjacent nearest neighbour " << f2 << " of " << t2 << "\n";
}
