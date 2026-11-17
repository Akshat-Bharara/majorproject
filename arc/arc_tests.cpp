// Correctness experiments (must pass 100% before any timing is reported).
//
//   arc_tests <rounds> <seed> <data-dir> <results-dir>
//
// Part A  n <= 16: at every candidate radius r at which the arc property
//         holds, D(r) (three implementations) against exhaustive search, for
//         all k in 2..n; helper routines (unimodality check, a(.), b(.)).
// Part B  n <= 16: O1 / O2 / adaptive against exhaustive search, all k.
// Part C  20 <= n <= 150: O1 / O2 / adaptive against the capped Delaunay DP
//         of ../convex/kdisp.hpp (exact oracle).
// Small instances and the values returned by O1/O2/adaptive are exported for
// the independent Python brute force (pybrute.py).
#include "arcdisp.hpp"
#include "gen.hpp"
#include "../convex/kdisp.hpp"
#include <fstream>
#include <iostream>
#include <map>

using namespace arc;

static std::vector<uint64_t> adj;
static bool has_clique(uint64_t cand, int need) {
    if (need <= 0) return true;
    if (__builtin_popcountll(cand) < need) return false;
    while (cand) {
        const int v = __builtin_ctzll(cand);
        cand &= cand - 1;
        if (has_clique(cand & adj[v], need - 1)) return true;
        if (__builtin_popcountll(cand) < need) return false;
    }
    return false;
}
static void build_adj(const Narrow& A, i64 r2) {
    adj.assign(A.n, 0);
    for (int i = 0; i < A.n; ++i)
        for (int j = 0; j < A.n; ++j)
            if (i != j && A.d(i, j) >= r2) adj[i] |= 1ULL << j;
}
static bool valid(const Narrow& A, std::vector<int> S, int k, i64 r2) {
    std::sort(S.begin(), S.end());
    if ((int)S.size() != k) return false;
    for (size_t i = 1; i < S.size(); ++i) if (S[i] == S[i - 1]) return false;
    return k < 2 || A.min_dist(S) >= r2;
}

struct Fam { long inst = 0, unimodal = 0, radii = 0, arc_radii = 0, dcheck = 0, opt = 0, adaptive_ok = 0, adaptive_na = 0, wrong_unchecked = 0, unchecked = 0; };
static std::map<std::string, Fam> stat;
static long failures = 0;
static void fail(const std::string& what, const Narrow& A, int k, i64 r2) {
    ++failures;
    std::cout << "FAIL " << what << " n=" << A.n << " k=" << k << " r2=" << r2 << " :";
    for (auto& q : A.p) std::cout << " (" << q.x << "," << q.y << ")";
    std::cout << std::endl;
}

static void small_instance(const std::string& fam, std::vector<P64> pts, std::ofstream* exp, const std::string& dir, const std::string& name) {
    Narrow A(pts);
    const int n = A.n;
    if (!A.strictly_convex()) { fail("generator: not strictly convex", A, 0, 0); return; }
    Fam& F = stat[fam];
    ++F.inst;
    const bool uni = A.unimodal_naive();
    if (uni != A.unimodal_fast()) fail("unimodal_fast != unimodal_naive", A, 0, 0);
    F.unimodal += uni;
    if (uni) A.compute_peaks();
    const std::vector<i64> cand = A.candidates();
    std::vector<i64> radii = cand;
    for (i64 c : cand) radii.push_back(c + 1); // also non-candidate radii
    const uint64_t all = (1ULL << n) - 1;
    std::vector<i64> opt(n + 1, 0); // brute-force optimum for every k
    for (i64 r2 : radii) {
        build_adj(A, r2);
        const bool arcp = A.arc_property(r2);
        ++F.radii;
        F.arc_radii += arcp;
        if (uni) {
            if (!arcp) fail("unimodal but arc property fails", A, 0, r2);
            for (int q = 0; q < n; ++q)
                if (A.a_uni(q, r2) != A.a_scan(q, r2) || A.b_uni(q, r2) != A.b_scan(q, r2)) fail("a_uni/b_uni", A, q, r2);
        }
        for (int k = 2; k <= n; ++k) {
            const bool truth = has_clique(all, k);
            if (truth) opt[k] = std::max(opt[k], r2);
            std::vector<int> w1, w2;
            const bool d1 = A.decide_chain(r2, k, false, &w1), d2 = A.decide_lift(r2, k, false, &w2);
            if (arcp) {
                ++F.dcheck;
                if (d1 != truth || d2 != truth) fail("D(r) != brute force under the arc property", A, k, r2);
                if (d1 && !valid(A, w1, k, r2)) fail("decide_chain witness", A, k, r2);
                if (d2 && !valid(A, w2, k, r2)) fail("decide_lift witness", A, k, r2);
                if (uni && (A.decide_chain(r2, k, true) != truth || A.decide_lift(r2, k, true) != truth)) fail("D(r) unimodal variant", A, k, r2);
            } else {
                ++F.unchecked; // what happens if D is used without the hypothesis
                F.wrong_unchecked += (d2 != truth);
            }
        }
    }
    if (exp) {
        std::ofstream f(dir + "/" + name);
        f << n << "\n";
        for (auto& q : pts) f << q.x << " " << q.y << "\n";
    }
    for (int k = 2; k <= n; ++k) {
        ++F.opt;
        if (uni) {
            const auto r1 = A.solve_O1(k, true, false);
            const auto r2 = A.solve_O2(k, 1000 + k);
            if (r1.r2 != opt[k] || !valid(A, r1.S, k, opt[k])) fail("O1 != brute force", A, k, r1.r2);
            if (r2.r2 != opt[k] || !valid(A, r2.S, k, opt[k])) fail("O2 != brute force", A, k, r2.r2);
            if (exp) { *exp << name << " " << k << " O2 " << r2.r2; for (int v : r2.S) *exp << " " << v; *exp << "\n"; }
        }
        const auto ra = A.solve_adaptive(k);
        if (ra.status == Status::OK) {
            ++F.adaptive_ok;
            if (ra.r2 != opt[k] || !valid(A, ra.S, k, opt[k])) fail("adaptive != brute force", A, k, ra.r2);
            if (exp) { *exp << name << " " << k << " adaptive " << ra.r2; for (int v : ra.S) *exp << " " << v; *exp << "\n"; }
        } else ++F.adaptive_na;
        if (uni && ra.status != Status::OK) fail("adaptive not applicable on unimodal input", A, k, 0);
    }
    if (exp) *exp << name << " flags " << (uni ? 1 : 0) << "\n";
}

struct Big { long inst = 0, unimodal = 0, checks = 0, adaptive_ok = 0, adaptive_na = 0; };
static std::map<std::string, Big> bstat;
static void big_instance(const std::string& fam, std::vector<P64> pts, std::mt19937_64& g) {
    Narrow A(pts);
    const int n = A.n;
    if (n < 20 || !A.strictly_convex()) return;
    std::vector<kd::Pt> kp;
    for (auto& q : pts) kp.push_back({q.x, q.y});
    kd::Instance I(kp);
    Big& B = bstat[fam];
    ++B.inst;
    const bool uni = A.unimodal_fast();
    if (n <= 120 && uni != A.unimodal_naive()) fail("unimodal_fast != unimodal_naive (big)", A, 0, 0);
    B.unimodal += uni;
    for (int k : {2, 3, 4, 5, 8, 13, n / 4, n / 2, n - 1, n}) {
        if (k < 2 || k > n) continue;
        const kd::Result ref = kd::solve_capped(I, k);
        if (uni) {
            const auto r1 = A.solve_O1(k, true, false);
            const auto r2 = A.solve_O2(k, g());
            ++B.checks;
            if (r1.r2 != ref.r2 || !valid(A, r1.S, k, ref.r2)) fail("O1 != DP oracle", A, k, r1.r2);
            if (r2.r2 != ref.r2 || !valid(A, r2.S, k, ref.r2)) fail("O2 != DP oracle", A, k, r2.r2);
        }
        const auto ra = A.solve_adaptive(k);
        if (ra.status == Status::OK) {
            ++B.adaptive_ok;
            if (ra.r2 != ref.r2 || !valid(A, ra.S, k, ref.r2)) fail("adaptive != DP oracle", A, k, ra.r2);
        } else ++B.adaptive_na;
    }
}

int main(int argc, char** argv) {
    const int rounds = argc > 1 ? atoi(argv[1]) : 60;
    std::mt19937_64 g(argc > 2 ? atoll(argv[2]) : 1);
    const std::string data = argc > 3 ? argv[3] : "data", res = argc > 4 ? argv[4] : "results";
    std::ofstream exp(res + "/small_cpp.txt");
    // the explicit counterexample of the earlier paper
    {
        Narrow A(pstar());
        long long beta = 0;
        const bool arcp = A.arc_property(10240, &beta);
        std::cout << "P*: strictly convex " << A.strictly_convex() << ", arc property at r2=10240: " << arcp << " (beta = " << beta
                  << "), unimodal: " << A.unimodal_naive() << ", D(10240, k=4) = " << A.decide_lift(10240, 4, false)
                  << ", adaptive(k=4): " << (A.solve_adaptive(4).status == Status::OK ? "OK" : "NOT_APPLICABLE") << "\n";
        if (arcp) { std::cout << "UNEXPECTED: P* has the arc property\n"; ++failures; }
        small_instance("pstar", pstar(), &exp, data, "small/pstar.txt");
    }
    for (int it = 0; it < rounds; ++it)
        for (int fam = 0; fam < NUM_FAMILIES; ++fam) {
            const int n = 4 + g() % 13; // 4..16
            std::vector<P64> pts = generate(g, fam, n);
            if ((int)pts.size() < 4) continue;
            const bool ex = it < 6;
            small_instance(family_name(fam), pts, ex ? &exp : nullptr, data,
                           std::string("small/") + family_name(fam) + "_" + std::to_string(it) + ".txt");
        }
    for (int it = 0; it < std::max(2, rounds / 12); ++it)
        for (int fam = 0; fam < NUM_FAMILIES; ++fam) {
            const int n = 20 + g() % 131;
            if ((fam == SQUARE || fam == GAUSS) && it > 2) continue;
            big_instance(family_name(fam), generate(g, fam, std::min(n, fam == SQUARE ? 40 : fam == GAUSS ? 24 : n)), g);
        }
    std::ofstream out(res + "/correctness.txt");
    auto both = [&](const std::string& s) { std::cout << s; out << s; };
    both("PART A/B (n <= 16)\nfamily instances unimodal radii arc_radii D_checks opt_checks adaptive_ok adaptive_na unchecked_D_tests unchecked_D_wrong\n");
    for (auto& [name, F] : stat)
        both(name + " " + std::to_string(F.inst) + " " + std::to_string(F.unimodal) + " " + std::to_string(F.radii) + " " + std::to_string(F.arc_radii) + " " +
             std::to_string(F.dcheck) + " " + std::to_string(F.opt) + " " + std::to_string(F.adaptive_ok) + " " + std::to_string(F.adaptive_na) + " " +
             std::to_string(F.unchecked) + " " + std::to_string(F.wrong_unchecked) + "\n");
    both("PART C (20 <= n <= 150, DP oracle)\nfamily instances unimodal O1O2_checks adaptive_ok adaptive_na\n");
    for (auto& [name, B] : bstat)
        both(name + " " + std::to_string(B.inst) + " " + std::to_string(B.unimodal) + " " + std::to_string(B.checks) + " " + std::to_string(B.adaptive_ok) + " " +
             std::to_string(B.adaptive_na) + "\n");
    both("FAILURES " + std::to_string(failures) + "\n");
    return failures ? 1 : 0;
}
