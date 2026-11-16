// In-process benchmark.  One CSV row per instance file.
//
//   bench [--reps R] [--triple-max N] [--uncapped-max N] [--brute-max N]
//         [--brute-kmax K] [--no-nowindow] file...
//
// Times are medians over R repetitions, in milliseconds, and include
// everything except reading the input and building the distance matrix
// (which every method shares).
#include "kdisp.hpp"
#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>

using namespace kd;

template <class F> static double time_ms(F&& f, int reps) {
    std::vector<double> t;
    for (int r = 0; r < reps; ++r) {
        const auto a = std::chrono::steady_clock::now();
        f();
        const auto b = std::chrono::steady_clock::now();
        t.push_back(std::chrono::duration<double, std::milli>(b - a).count());
    }
    std::sort(t.begin(), t.end());
    return t[t.size() / 2];
}

static bool valid(const Instance& I, const std::vector<int>& S, int k, i64 r2) {
    if ((int)S.size() != k) return false;
    for (size_t a = 1; a < S.size(); ++a)
        if (S[a] <= S[a - 1]) return false;
    return min_dist2(I, S) == r2;
}

static std::string setstr(const std::vector<int>& S) {
    std::string s;
    for (size_t i = 0; i < S.size(); ++i) s += (i ? " " : "") + std::to_string(S[i]);
    return s;
}

int main(int argc, char** argv) {
    int reps = 3, triple_max = 100, uncapped_max = 200, brute_max = 16, brute_kmax = 64;
    bool nowindow = true;
    std::vector<std::string> files;
    for (int a = 1; a < argc; ++a) {
        if (!strcmp(argv[a], "--reps")) reps = atoi(argv[++a]);
        else if (!strcmp(argv[a], "--triple-max")) triple_max = atoi(argv[++a]);
        else if (!strcmp(argv[a], "--uncapped-max")) uncapped_max = atoi(argv[++a]);
        else if (!strcmp(argv[a], "--brute-max")) brute_max = atoi(argv[++a]);
        else if (!strcmp(argv[a], "--brute-kmax")) brute_kmax = atoi(argv[++a]);
        else if (!strcmp(argv[a], "--no-nowindow")) nowindow = false;
        else files.push_back(argv[a]);
    }
    std::cout << "file,n,k,r2,ok,t_capped,t_nowindow,t_uncapped,t_triple,t_certified,"
                 "dec_window,dec_nowindow,dec_certified,greedy_calls,certified_by_greedy,"
                 "cand_total,cand_window,far_triples_at_opt,incircle_per_decision,"
                 "r2_brute,r2_triple,r2_uncapped,r2_certified,r2_nowindow,greedy_ff_r2,set_capped,set_certified\n";
    for (const std::string& file : files) {
        std::ifstream in(file);
        Instance I;
        int k;
        if (!read_instance(in, I, k)) { std::cerr << "cannot read " << file << "\n"; return 1; }
        I.validate();
        const int n = I.n;
        Result rc, rn, ru, rt, rz;
        const double tc = time_ms([&] { rc = solve_capped(I, k, true); }, reps);
        double tn = -1, tu = -1, tt = -1;
        if (nowindow) tn = time_ms([&] { rn = solve_capped(I, k, false); }, reps);
        if (n <= uncapped_max) tu = time_ms([&] { ru = solve_capped(I, k, true, n); }, reps);
        if (n <= triple_max) tt = time_ms([&] { rt = solve_triple(I, k); }, reps);
        const double tz = time_ms([&] { rz = solve_certified(I, k); }, reps);
        i64 rb = -1;
        if (n <= brute_max || k <= brute_kmax) { Brute B(I, k); rb = B.solve(); }
        bool ok = valid(I, rc.S, k, rc.r2) && valid(I, rz.S, k, rc.r2) && rz.r2 == rc.r2;
        if (nowindow) ok = ok && rn.r2 == rc.r2 && valid(I, rn.S, k, rc.r2);
        if (tu >= 0) ok = ok && ru.r2 == rc.r2 && valid(I, ru.S, k, rc.r2);
        if (tt >= 0) ok = ok && rt.r2 == rc.r2;
        if (rb >= 0) ok = ok && rb == rc.r2;
        // number of r-far triples at the optimum (the T_r of the analysis)
        u64 tr = 0;
        for (int i = 0; i < n; ++i)
            for (int j = i + 2; j < n; ++j)
                if (I.d(i, j) >= rc.r2)
                    for (int l = i + 1; l < j; ++l) tr += I.d(i, l) >= rc.r2 && I.d(l, j) >= rc.r2;
        const i64 gff = min_dist2(I, farthest_first(I, k));
        std::cout << file << ',' << n << ',' << k << ',' << rc.r2 << ',' << ok << ',' << tc << ',' << tn << ','
                  << tu << ',' << tt << ',' << tz << ',' << rc.decisions << ',' << rn.decisions << ','
                  << rz.decisions << ',' << rz.greedy_calls << ',' << rz.certified_by_greedy << ','
                  << rc.cand_total << ',' << rc.cand_window << ',' << tr << ','
                  << (rc.decisions ? rc.incircle / rc.decisions : 0) << ',' << rb << ','
                  << (tt >= 0 ? rt.r2 : -1) << ',' << (tu >= 0 ? ru.r2 : -1) << ',' << rz.r2 << ','
                  << (nowindow ? rn.r2 : -1) << ',' << gff << ',' << setstr(rc.S) << ',' << setstr(rz.S)
                  << std::endl;
    }
    return 0;
}
