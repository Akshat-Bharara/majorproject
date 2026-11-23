// sbm_test check <rounds> <seed>: SBM decision/optimum and SBM DP vs brute force on random instances
#include "sbm.hpp"
#include "../convex/randgen.hpp"
#include <iostream>
int main(int argc, char** argv) {
    const int rounds = argc > 1 ? atoi(argv[1]) : 1000;
    std::mt19937_64 g(argc > 2 ? atoll(argv[2]) : 1);
    long dec = 0, decfail = 0, opt = 0, optfail = 0, dp = 0, dpfail = 0;
    for (int it = 0; it < rounds; ++it) {
        kd::Instance I(kd::random_convex(g, it % kd::NUM_FAMILIES, 5 + g() % 9));
        if (I.n < 5) continue;
        for (int k = 3; k <= std::min(I.n, 7); ++k) {
            kd::Brute B(I, k);
            const kd::i64 o = B.solve();
            for (kd::i64 r2 : I.candidates()) {
                const bool truth = r2 <= o;
                ++dec;
                if (sbm::sbm_decide(I, k, r2) != truth) ++decfail;
            }
            ++opt; if (sbm::sbm_optimize(I, k) != o) ++optfail;
            ++dp; sbm::SbmDP D(I, k); if (D.solve() != o) { ++dpfail; if (dpfail < 3) std::cout << "DP mismatch n=" << I.n << " k=" << k << " got " << D.solve() << " want " << o << "\n"; }
        }
    }
    std::cout << "decision tests " << dec << " wrong " << decfail << " | optimum tests " << opt << " wrong " << optfail
              << " | DP tests " << dp << " wrong " << dpfail << "\n";
}
