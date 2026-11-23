// Original (SBM, DAM 2026) vs. ours on the same instances; one CSV row per run.
//   sbm_bench <n> <k> <seed> <family>
// algorithms: SBM-DP (Thm 6, O(n^4 k^2)), SBM-FPT (Alg. 1 + search tree),
// ours: capped DP (O(n^3 log k log n)), certified, and (if unimodal) O2.
#include "sbm.hpp"
#include "../convex/randgen.hpp"
#include "../arc/arcdisp.hpp"
#include <chrono>
#include <iostream>
using namespace kd;
static double ms(std::chrono::steady_clock::time_point a) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
}
int main(int argc, char** argv) {
    const int n = atoi(argv[1]), k = atoi(argv[2]);
    std::mt19937_64 g(atoll(argv[3]));
    const int fam = atoi(argv[4]);
    std::vector<Pt> p;
    for (int t = 0; t < 50 && (int)p.size() != n; ++t) p = random_convex(g, fam, n);
    if ((int)p.size() != n) return 0;
    Instance I(p);
    I.validate();
    auto t = std::chrono::steady_clock::now();
    const i64 ref = solve_capped(I, k).r2;
    const double tc = ms(t);
    t = std::chrono::steady_clock::now();
    const i64 rz = solve_certified(I, k).r2;
    const double tz = ms(t);
    double td = -1, tf = -1;
    i64 rd = -1, rf = -1;
    if (n <= 60) { t = std::chrono::steady_clock::now(); sbm::SbmDP D(I, k); rd = D.solve(); td = ms(t); }
    if (k <= 16 && n <= 400) { t = std::chrono::steady_clock::now(); rf = sbm::sbm_optimize(I, k); tf = ms(t); }
    std::cout << n << ',' << k << ',' << family_name(fam) << ',' << tc << ',' << tz << ',' << td << ',' << tf << ','
              << (rz == ref) << ',' << (rd < 0 || rd == ref) << ',' << (rf < 0 || rf == ref) << std::endl;
}
