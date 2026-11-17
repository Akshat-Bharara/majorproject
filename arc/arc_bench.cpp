// Scaling experiments on unimodal instances.
//
//   arc_bench feas
//       which unimodal instances exist in the 2^29 box (rounded circles,
//       polygons built from all short primitive vectors)
//   arc_bench run <n> <k> <narrow|wide> <reps> <seed>
//       one CSV line per algorithm and repetition:
//       n,k,mode,algo,rep,time_ms,decisions,r2,ok
// narrow: |coord| <= 2^29, 64-bit squared distances; O1 for n <= 10000 and
//         the capped Delaunay DP (../convex/kdisp.hpp) for n <= 400 are run too.
// wide:   rounded circle of radius 2^60, 128-bit squared distances; O2 only.
// Every instance is verified to be strictly convex and unimodal (exact
// O(n log n) test) before timing; every witness is verified by an exact
// closest-pair computation.
#include "arcdisp.hpp"
#include "gen.hpp"
#include "../convex/kdisp.hpp"
#include <chrono>
#include <cstring>
#include <iostream>
#include <numeric>
#include <set>

using namespace arc;

static std::string str128(i128 v) {
    if (v == 0) return "0";
    std::string s;
    bool neg = v < 0;
    if (neg) v = -v;
    while (v) { s.push_back('0' + (int)(v % 10)); v /= 10; }
    if (neg) s.push_back('-');
    std::reverse(s.begin(), s.end());
    return s;
}
// exact minimum squared distance of the points with indices S (sweep line)
template <class A> static i128 closest_pair(const A& I, std::vector<int> S) {
    std::sort(S.begin(), S.end(), [&](int a, int b) { return I.p[a].x != I.p[b].x ? I.p[a].x < I.p[b].x : I.p[a].y < I.p[b].y; });
    auto cmp = [&](int a, int b) { return I.p[a].y != I.p[b].y ? I.p[a].y < I.p[b].y : a < b; };
    std::set<int, decltype(cmp)> act(cmp);
    i128 best = -1;
    size_t tail = 0;
    for (size_t i = 0; i < S.size(); ++i) {
        const int v = S[i];
        while (best >= 0 && tail < i) {
            const i128 dx = (i128)I.p[v].x - I.p[S[tail]].x;
            if (dx * dx > best) act.erase(S[tail++]); else break;
        }
        auto it = act.lower_bound(v);
        for (auto up = it; up != act.end(); ++up) {
            const i128 dy = (i128)I.p[*up].y - I.p[v].y;
            if (best >= 0 && dy * dy > best) break;
            const i128 dd = (i128)I.d(v, *up);
            if (best < 0 || dd < best) best = dd;
        }
        for (auto dn = it; dn != act.begin();) {
            --dn;
            const i128 dy = (i128)I.p[v].y - I.p[*dn].y;
            if (best >= 0 && dy * dy > best) break;
            const i128 dd = (i128)I.d(v, *dn);
            if (best < 0 || dd < best) best = dd;
        }
        act.insert(v);
    }
    return best;
}
static double now_ms() {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
static std::vector<P64> primitive_polygon(int rho) {
    std::vector<P64> v;
    for (int a = -rho; a <= rho; ++a)
        for (int b = -rho; b <= rho; ++b)
            if ((a || b) && a * a + b * b <= rho * rho && std::gcd(std::abs(a), std::abs(b)) == 1) v.push_back({a, b});
    auto half = [](const P64& u) { return u.y > 0 || (u.y == 0 && u.x > 0) ? 0 : 1; };
    std::sort(v.begin(), v.end(), [&](const P64& s, const P64& t) {
        if (half(s) != half(t)) return half(s) < half(t);
        return s.x * t.y - s.y * t.x > 0;
    });
    std::vector<P64> p;
    P64 cur{0, 0};
    for (auto& e : v) { p.push_back(cur); cur.x += e.x; cur.y += e.y; }
    return p;
}

template <class A> static void run(A& I, int k, const char* mode, int reps, uint64_t seed, bool narrow) {
    const int n = I.n;
    if (!I.strictly_convex() || !I.unimodal_fast()) { std::cout << "# instance rejected: not strictly convex / not unimodal\n"; return; }
    std::string ref;
    auto report = [&](const char* algo, int rep, double ms, int dec, const std::string& r2, const std::vector<int>& S, i128 r2v) {
        bool ok = (int)S.size() == k && std::adjacent_find(S.begin(), S.end()) == S.end();
        ok = ok && (k < 2 || closest_pair(I, S) == r2v);
        if (ref.empty()) ref = r2;
        ok = ok && ref == r2;
        std::cout << n << ',' << k << ',' << mode << ',' << algo << ',' << rep << ',' << ms << ',' << dec << ',' << r2 << ',' << ok << std::endl;
    };
    for (int rep = 0; rep < reps; ++rep) {
        I.have_peaks = false; // peaks are part of the measured time
        const double t0 = now_ms();
        const auto R = I.solve_O2(k, seed + 7919 * rep);
        const double t1 = now_ms();
        report("O2", rep, t1 - t0, R.decisions, str128(R.r2), R.S, R.r2);
    }
    if (narrow && n <= 10000) {
        for (int rep = 0; rep < (n <= 1000 ? reps : 1); ++rep) {
            I.have_peaks = false;
            const double t0 = now_ms();
            const auto R = I.solve_O1(k, true, false);
            const double t1 = now_ms();
            report("O1", rep, t1 - t0, R.decisions, str128(R.r2), R.S, R.r2);
        }
    }
}

int main(int argc, char** argv) {
    if (argc >= 2 && !strcmp(argv[1], "feas")) {
        std::mt19937_64 g(5);
        std::cout << "rounded circle, radius 5e8 (inside the 2^29 box): n convex unimodal\n";
        for (int n : {1000, 10000, 20000, 30000, 50000, 100000}) {
            Narrow A(round_circle(g, n, 5.0e8L));
            const bool c = A.strictly_convex();
            std::cout << n << ' ' << c << ' ' << (c ? A.unimodal_fast() : false) << std::endl;
        }
        std::cout << "polygon whose edges are all primitive vectors of norm <= rho: rho n max|coord| convex unimodal\n";
        for (int rho : {8, 23, 72, 229, 724}) {
            Narrow A(primitive_polygon(rho));
            i64 m = 0;
            for (auto& q : A.p) m = std::max({m, std::llabs(q.x), std::llabs(q.y)});
            std::cout << rho << ' ' << A.n << ' ' << m << ' ' << A.strictly_convex() << ' ' << A.unimodal_fast() << std::endl;
        }
        return 0;
    }
    if (argc < 7 || strcmp(argv[1], "run")) { std::cerr << "usage: arc_bench feas | run n k narrow|wide reps seed\n"; return 2; }
    const int n = atoi(argv[2]), k = atoi(argv[3]), reps = atoi(argv[5]);
    const uint64_t seed = atoll(argv[6]);
    std::mt19937_64 g(seed);
    if (!strcmp(argv[4], "narrow")) {
        std::vector<P64> pts = round_circle(g, n, 5.0e8L);
        Narrow A(pts);
        run(A, k, "narrow", reps, seed, true);
        if (n <= 400) { // cubic baseline: capped Delaunay DP
            std::vector<kd::Pt> kp;
            for (auto& q : pts) kp.push_back({q.x, q.y});
            kd::Instance I(kp);
            for (int rep = 0; rep < reps; ++rep) {
                const double t0 = now_ms();
                const kd::Result R = kd::solve_capped(I, k);
                const double t1 = now_ms();
                const auto chk = A.solve_O1(k, true, false);
                std::cout << n << ',' << k << ",narrow,DP," << rep << ',' << (t1 - t0) << ',' << R.decisions << ',' << R.r2 << ',' << (R.r2 == chk.r2) << std::endl;
            }
        }
    } else {
        Wide A(round_circle(g, n, 1152921504606846976.0L)); // 2^60
        run(A, k, "wide", reps, seed, false);
    }
    return 0;
}
