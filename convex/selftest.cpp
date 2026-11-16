// Randomised self-test: every solver against exhaustive search.
//   selftest [rounds] [seed]
#include "kdisp.hpp"
#include "randgen.hpp"
#include <iostream>

using namespace kd;

static bool valid(const Instance& I, const std::vector<int>& S, int k, i64 r2) {
    if ((int)S.size() != k) return false;
    for (size_t a = 1; a < S.size(); ++a)
        if (S[a] <= S[a - 1]) return false;
    return min_dist2(I, S) == r2;
}

int main(int argc, char** argv) {
    const int rounds = argc > 1 ? atoi(argv[1]) : 2000;
    std::mt19937_64 rng(argc > 2 ? atoll(argv[2]) : 12345);
    long checks = 0;
    for (int it = 0; it < rounds; ++it) {
        const int fam = it % NUM_FAMILIES;
        const int target = 4 + rng() % 13;
        Instance I(random_convex(rng, fam, target));
        const int n = I.n;
        if (n < 4) continue;
        I.validate();
        for (int k = 2; k <= std::min(n, 9); ++k) {
            Brute B(I, k);
            const i64 opt = B.solve();
            const Result a = solve_capped(I, k, true);
            const Result b = solve_capped(I, k, false);
            const Result c = solve_certified(I, k);
            const Result d = solve_triple(I, k);
            const Result e = solve_capped(I, k, true, n); // uncapped ablation
            if (a.r2 != opt || b.r2 != opt || c.r2 != opt || d.r2 != opt || e.r2 != opt ||
                !valid(I, a.S, k, opt) || !valid(I, b.S, k, opt) || !valid(I, c.S, k, opt) ||
                !valid(I, e.S, k, opt)) {
                std::cout << "MISMATCH family " << fam << " n " << n << " k " << k << " brute " << opt
                          << " capped " << a.r2 << " nowin " << b.r2 << " cert " << c.r2 << " triple "
                          << d.r2 << " uncapped " << e.r2 << "\n";
                for (auto& q : I.p) std::cout << q.x << " " << q.y << "\n";
                return 1;
            }
            ++checks;
        }
        // maximum r-far subset against exhaustive search, a few radii
        const auto C = I.candidates();
        for (int t = 0; t < 4; ++t) {
            const i64 r2 = C[rng() % C.size()];
            const std::vector<int> M = max_far_set(I, r2);
            int alpha = 1;
            for (int k = 2; k <= n; ++k) {
                Brute B(I, k);
                if (B.solve() >= r2) alpha = k;
                else break;
            }
            if ((int)M.size() != alpha || (M.size() > 1 && min_dist2(I, M) < r2)) {
                std::cout << "MIS MISMATCH n " << n << " r2 " << r2 << " got " << M.size() << " want "
                          << alpha << "\n";
                return 1;
            }
            ++checks;
        }
    }
    std::cout << "selftest ok: " << checks << " checks\n";
    return 0;
}
