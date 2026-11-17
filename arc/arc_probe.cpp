// Probing the hypotheses of Lemma A1 / Corollary A2.
//
//   arc_probe small <instances> <seed>     exhaustive, n <= 13, k = 3..6
//   arc_probe stats <instances> <seed>     20 <= n <= 60, r* from the DP oracle
//
// greedy(r): some chain q, a(q), a(a(q)), ... of k points does not wrap and
//            is r-far (checked pairwise, no assumption).
// H1  if some r-far k-set S has Far_r(s) arc-shaped for all s in S, then greedy(r).
// H2  if the arc property holds at r = r*, binary search over the candidate
//     radii with greedy(.) as predicate returns r*.   (H2u: with the unchecked D.)
// H3  greedy(r*) holds (no assumption).
#include "arcdisp.hpp"
#include "gen.hpp"
#include "../convex/kdisp.hpp"
#include <iostream>
#include <map>

using namespace arc;

static bool greedy(const Narrow& A, i64 r2, int k) {
    const int n = A.n;
    std::vector<int> ta(n);
    for (int q = 0; q < n; ++q) ta[q] = A.a_scan(q, r2);
    std::vector<int> ch(k);
    for (int q = 0; q < n; ++q) {
        int off = 0, x = q;
        bool ok = true;
        ch[0] = q;
        for (int i = 1; i < k && ok; ++i) {
            if (!ta[x] || off + ta[x] >= n) ok = false;
            else { off += ta[x]; x = A.idx(x, ta[x]); ch[i] = x; }
        }
        for (int i = 0; i < k && ok; ++i)
            for (int j = i + 1; j < k && ok; ++j) ok = A.d(ch[i], ch[j]) >= r2;
        if (ok) return true;
    }
    return false;
}
template <class Pred> static i64 bsearch(const std::vector<i64>& cand, Pred ok) {
    size_t lo = 0, hi = cand.size() - 1;
    while (lo < hi) {
        const size_t mid = lo + (hi - lo + 1) / 2;
        if (ok(cand[mid])) lo = mid; else hi = mid - 1;
    }
    return cand[lo];
}
static std::string dump(const Narrow& A) {
    std::string s;
    for (auto& q : A.p) s += "(" + std::to_string(q.x) + "," + std::to_string(q.y) + ") ";
    return s;
}

struct H { long premise = 0, viol = 0; int best_n = 99; std::string cex; };
static std::map<std::string, H> hyp;
static void rec(const std::string& name, bool violated, const Narrow& A, int k, i64 r2) {
    H& h = hyp[name];
    ++h.premise;
    if (violated) {
        ++h.viol;
        if (A.n < h.best_n) { h.best_n = A.n; h.cex = "n=" + std::to_string(A.n) + " k=" + std::to_string(k) + " r2=" + std::to_string(r2) + " : " + dump(A); }
    }
}

static std::vector<uint64_t> adj;
// does an r-far k-set exist all of whose points are in `good`?  (clique search)
static bool clique(uint64_t cand, int need) {
    if (need <= 0) return true;
    while (cand) {
        if (__builtin_popcountll(cand) < need) return false;
        const int v = __builtin_ctzll(cand);
        cand &= cand - 1;
        if (clique(cand & adj[v], need - 1)) return true;
    }
    return false;
}
static void small(const std::vector<P64>& pts) {
    Narrow A(pts);
    const int n = A.n;
    if (n < 6 || !A.strictly_convex()) return;
    const std::vector<i64> cand = A.candidates();
    const uint64_t all = (1ULL << n) - 1;
    for (int k = 3; k <= std::min(6, n - 1); ++k) {
        i64 opt = 0;
        for (i64 r2 : cand) {
            adj.assign(n, 0);
            for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) if (i != j && A.d(i, j) >= r2) adj[i] |= 1ULL << j;
            if (!clique(all, k)) break;
            opt = r2;
            uint64_t good = 0;
            for (int q = 0; q < n; ++q) if (A.far_arcs(q, r2) <= 1) good |= 1ULL << q;
            if (clique(good, k)) rec("H1 k=" + std::to_string(k), !greedy(A, r2, k), A, k, r2);
        }
        rec("H3 k=" + std::to_string(k), !greedy(A, opt, k), A, k, opt);
        if (A.arc_property(opt)) {
            rec("H2 k=" + std::to_string(k), bsearch(cand, [&](i64 r) { return greedy(A, r, k); }) != opt, A, k, opt);
            rec("H2u k=" + std::to_string(k), bsearch(cand, [&](i64 r) { return A.decide_lift(r, k, false); }) != opt, A, k, opt);
        }
    }
}

struct Row { long inst = 0, arcprop = 0, greedy_ok = 0, adaptive_ok = 0, unchecked_ok = 0; double frac = 0, beta = 0; };
int main(int argc, char** argv) {
    const std::string mode = argc > 1 ? argv[1] : "small";
    const long count = argc > 2 ? atol(argv[2]) : 1000;
    std::mt19937_64 g(argc > 3 ? atoll(argv[3]) : 1);
    if (mode == "small") {
        small(pstar());
        for (long it = 0; it < count; ++it) {
            const int fam = (it % 3 == 0) ? TWOCHAIN : (it % 3 == 1) ? ELLIPSE : (int)(g() % NUM_FAMILIES);
            std::vector<P64> pts = generate(g, fam, 6 + g() % 8);
            if (it % 5 == 4) { // local perturbation of the previous shape family: jitter the points
                for (auto& q : pts) { q.x += (i64)(g() % 2001) - 1000; q.y += (i64)(g() % 201) - 100; }
                pts = strict_hull(pts);
            }
            small(pts);
        }
        std::cout << "hypothesis premise_holds violations smallest_counterexample\n";
        for (auto& [name, h] : hyp) std::cout << name << " " << h.premise << " " << h.viol << " " << (h.viol ? h.cex : "-") << "\n";
        return 0;
    }
    // stats: structure at r = r* on general instances
    std::map<std::string, Row> rows, buckets;
    for (long it = 0; it < count; ++it) {
        const int fam = it % (NUM_FAMILIES + 1);
        std::vector<P64> pts;
        if (fam == NUM_FAMILIES) { // ellipses of moderate aspect ratio 1:1.5 .. 1:10
            const int m = 20 + (int)(g() % 41);
            const double asp = 0.1 + 0.57 * urand(g);
            for (int i = 0; i < m; ++i) {
                const double t = 6.283185307179586 * urand(g);
                pts.push_back({(i64)llround(1e6 * cos(t)), (i64)llround(1e6 * asp * sin(t))});
            }
            pts = strict_hull(pts);
        } else pts = generate(g, fam, fam == GAUSS ? 20 : fam == SQUARE ? 30 : 20 + (int)(g() % 41));
        Narrow A(pts);
        const int n = A.n;
        if (n < 12 || !A.strictly_convex()) continue;
        std::vector<kd::Pt> kp;
        for (auto& q : pts) kp.push_back({q.x, q.y});
        kd::Instance I(kp);
        const std::vector<i64> cand = A.candidates();
        for (int k : {4, 8}) {
            const i64 opt = kd::solve_capped(I, k).r2;
            long long beta = 0;
            int good = 0;
            const bool ap = A.arc_property(opt, &beta, &good);
            const bool gr = greedy(A, opt, k);
            const bool ad = A.solve_adaptive(k).status == Status::OK;
            const bool un = bsearch(cand, [&](i64 r) { return A.decide_lift(r, k, false); }) == opt;
            const std::string b = beta == 0 ? "beta=0" : beta <= 2 ? "beta 1-2" : beta <= 10 ? "beta 3-10" : beta <= 40 ? "beta 11-40" : "beta >40";
            for (Row* r : {&rows[fam == NUM_FAMILIES ? "ellipse-mid" : family_name(fam)], &buckets[b]}) {
                ++r->inst; r->arcprop += ap; r->greedy_ok += gr; r->adaptive_ok += ad; r->unchecked_ok += un;
                r->frac += (double)good / n; r->beta += (double)beta;
            }
        }
    }
    auto print = [](const std::map<std::string, Row>& m) {
        for (auto& [name, r] : m)
            std::cout << name << " | " << r.inst << " | " << r.frac / r.inst << " | " << r.beta / r.inst << " | " << r.arcprop << " | " << r.greedy_ok
                      << " | " << r.adaptive_ok << " | " << r.unchecked_ok << "\n";
    };
    std::cout << "group | tests | mean fraction of arc-shaped points at r* | mean beta | arc property at r* | greedy(r*) | adaptive applicable | unchecked O1 = r*\n";
    print(rows);
    std::cout << "--\n";
    print(buckets);
    return 0;
}
