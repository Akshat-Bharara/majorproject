// Exhaustive / randomised tests of the structural properties of Section 6.
//
//   hunt exh   <instances> <seed> [nmax=12] [kmax=6]
//       small instances; every candidate radius; all properties, decided by
//       exhaustive clique search (independent of the dynamic program).
//   hunt big   <instances> <seed> [nmax=30] [kmax=12]
//       larger instances; radii r*_k and random feasible ones; C1 by the
//       cyclic capped DP, V4/V5 by search.
//   hunt tiny  <instances> <seed> [W=60] [kmax=6]
//       exhaustive tests on thin two-chain polygons with small coordinates.
//   hunt pert  <instances> <seed> <file>
//       like exh+big, on random perturbations of the instances in <file>
//       (blocks "n" + n lines "x y").
//
// Output:  "STAT <mode> <prop> <k> <tests> <failures>" lines and, for the
// smallest failing instance of each (prop,k), a "CEX ..." line.
//
// Properties (r fixed, far = distance >= r, a(q)/b(q) = first far point
// counterclockwise/clockwise of q):
//   G   an r-far k-set exists  <=>  some arc A(i)=[a(i),b(i)] has an r-far (k-1)-set
//   C1  (feasible) some r-far k-set contains q, a(q)
//   C1p (feasible) some r-far k-set contains q, a(q), b(q)
//   V1  alternating greedy from {q,a(q),b(q)} succeeds for some q
//   V3  search tree from {q,a(q),b(q)} succeeds for some q
//   V5  alternating greedy from {q,a(q)} succeeds for some q
//   V4  search tree from {q,a(q)} succeeds for some q
//   X   for EVERY r-far k-set some single exchange s_{i+1} -> a(s_i) keeps it r-far
//   T3  (k=3) for EVERY r-far triple (x,y,z) one of (x,a(x),z), (y,a(y),x),
//       (z,a(z),y) is r-far   [Theorem: always true]
#include "kdisp.hpp"
#include "randgen.hpp"
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

using namespace kd;

struct Stat {
    long tests = 0, fails = 0;
    int best_n = 1 << 30;
    std::string cex;
};
static std::map<std::pair<std::string, int>, Stat> stats;

static std::string dump(const Instance& I, int k, i64 r2) {
    std::ostringstream o;
    o << "n " << I.n << " k " << k << " r2 " << r2 << " :";
    for (auto& q : I.p) o << " (" << q.x << "," << q.y << ")";
    return o.str();
}
static void record(const std::string& prop, int k, bool ok, const Instance& I, i64 r2) {
    Stat& s = stats[{prop, k}];
    ++s.tests;
    if (!ok) {
        ++s.fails;
        if (I.n < s.best_n) { s.best_n = I.n; s.cex = dump(I, k, r2); }
    }
}

// ---- exhaustive clique search on bitmasks (n <= 62) ----
static std::vector<u64> adj;
static bool has_clique(u64 cand, int need) {
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
// X: for EVERY r-far k-set S (ccw s_1..s_k) there is an i such that replacing
// s_{i+1} by a(s_i) gives an r-far k-set.  (True for k = 3; this is T3.)
static bool exchange_ok(const Instance& I, const std::vector<int>& ta, std::vector<int>& S, int next, int k) {
    const int n = I.n;
    if ((int)S.size() == k) {
        for (int i = 0; i < k; ++i) {
            const int a = (S[i] + ta[S[i]]) % n;
            bool ok = true;
            for (int j = 0; j < k && ok; ++j)
                if (j != i && j != (i + 1) % k) ok = a != S[j] && (adj[a] >> S[j] & 1);
            if (ok) return true;
        }
        return false;
    }
    for (int v = next; v < n; ++v) {
        bool ok = true;
        for (int u : S) ok = ok && (adj[u] >> v & 1);
        if (!ok) continue;
        S.push_back(v);
        const bool r = exchange_ok(I, ta, S, v + 1, k);
        S.pop_back();
        if (!r) return false;
    }
    return true;
}

static void build_adj(const Instance& I, i64 r2) {
    adj.assign(I.n, 0);
    for (int i = 0; i < I.n; ++i)
        for (int j = 0; j < I.n; ++j)
            if (i != j && I.d(i, j) >= r2) adj[i] |= 1ULL << j;
}

static void seeded_props(const Instance& I, int k, i64 r2) {
    const int n = I.n;
    Seeded sd(I, r2, k);
    bool v1 = false, v3 = false, v5 = false, v4 = false;
    for (int q = 0; q < n; ++q) {
        if (!v1) v1 = sd.alternating(q, true, true) || sd.alternating(q, true, false);
        if (!v5) v5 = sd.alternating(q, false, true) || sd.alternating(q, false, false);
        if (!v3) v3 = sd.tree(q, true);
        if (!v4) v4 = sd.tree(q, false);
    }
    record("V1", k, v1, I, r2);
    record("V3", k, v3, I, r2);
    record("V5", k, v5, I, r2);
    record("V4", k, v4, I, r2);
}

static void test_exhaustive(const Instance& I, int kmax) {
    const int n = I.n;
    const u64 all = n == 64 ? ~0ULL : (1ULL << n) - 1;
    CappedDP dp(I);
    for (i64 r2 : I.candidates()) {
        build_adj(I, r2);
        std::vector<int> ta(n), tb(n);
        for (int q = 0; q < n; ++q) { ta[q] = first_far(I, r2, q, +1); tb[q] = first_far(I, r2, q, -1); }
        // T3: every far triple can be tightened
        bool t3 = true, any3 = false;
        for (int x = 0; x < n && t3; ++x)
            for (int y = x + 1; y < n && t3; ++y)
                for (int z = y + 1; z < n && t3; ++z) {
                    if (!(adj[x] >> y & 1) || !(adj[x] >> z & 1) || !(adj[y] >> z & 1)) continue;
                    any3 = true;
                    const int ax = (x + ta[x]) % n, ay = (y + ta[y]) % n, az = (z + ta[z]) % n;
                    t3 = (adj[ax] >> z & 1) || (adj[ay] >> x & 1) || (adj[az] >> y & 1);
                }
        if (any3) record("T3", 3, t3, I, r2);
        for (int k = 3; k <= std::min(kmax, n); ++k) {
            const bool feas = has_clique(all, k);
            bool rhs = false;
            for (int i = 0; i < n && !rhs; ++i) {
                if (!ta[i]) continue;
                u64 arc = 0;
                for (int t = ta[i]; t <= n - tb[i]; ++t) arc |= 1ULL << ((i + t) % n);
                rhs = has_clique(arc, k - 1);
            }
            record("G", k, feas == rhs, I, r2);
            if (feas && !rhs) record("G=>", k, false, I, r2); // never happens (Lemma)
            if (!feas) continue;
            bool c1 = false, c1p = false;
            for (int q = 0; q < n; ++q) {
                if (!ta[q]) continue;
                const int a = (q + ta[q]) % n, b = ((q - tb[q]) % n + n) % n;
                if (!c1) c1 = has_clique(adj[q] & adj[a], k - 2);
                if (!c1p && a != b && (adj[a] >> b & 1)) c1p = has_clique(adj[q] & adj[a] & adj[b], k - 3);
            }
            record("C1", k, c1, I, r2);
            record("C1p", k, c1p, I, r2);
            if (k <= 5) { std::vector<int> S; record("X", k, exchange_ok(I, ta, S, 0, k), I, r2); }
            // the DP-based C1 test used in `big` mode must agree
            dp.run(r2, k, true, false);
            bool c1dp = false;
            for (int q = 0; q < n && !c1dp; ++q) {
                if (!ta[q]) continue;
                const int a = (q + ta[q]) % n;
                const int j = q > a ? q : q + n;
                if (j - a >= 2) c1dp = dp.count(a, j) >= k - 2;
            }
            record("C1-dp-agrees", k, c1dp == c1, I, r2);
            seeded_props(I, k, r2);
        }
    }
}

static void test_big(const Instance& I, int kmax, std::mt19937_64& rng) {
    const int n = I.n;
    const auto C = I.candidates();
    CappedDP dp(I);
    for (int k = 3; k <= std::min(kmax, n - 1); ++k) {
        const Result R = solve_capped(I, k);
        const size_t top = std::lower_bound(C.begin(), C.end(), R.r2) - C.begin();
        std::vector<i64> radii{R.r2};
        for (int t = 0; t < 2; ++t) radii.push_back(C[top - std::min<size_t>(top, rng() % 12)]);
        std::sort(radii.begin(), radii.end());
        radii.erase(std::unique(radii.begin(), radii.end()), radii.end());
        for (i64 r2 : radii) {
            dp.run(r2, k, true, false);
            bool c1 = false;
            for (int q = 0; q < n && !c1; ++q) {
                const int t = first_far(I, r2, q, +1);
                if (!t) continue;
                const int a = (q + t) % n;
                const int j = q > a ? q : q + n;
                if (j - a >= 2) c1 = dp.count(a, j) >= k - 2;
            }
            record("C1", k, c1, I, r2);
            Seeded sd(I, r2, k);
            bool v5 = false, v4 = false;
            for (int q = 0; q < n; ++q) {
                if (!v5) v5 = sd.alternating(q, false, true) || sd.alternating(q, false, false);
                if (!v4) v4 = sd.tree(q, false);
            }
            record("V5", k, v5, I, r2);
            record("V4", k, v4, I, r2);
        }
    }
}

static std::vector<std::vector<Pt>> read_bases(const char* file) {
    std::vector<std::vector<Pt>> out;
    std::ifstream in(file);
    int n;
    while (in >> n) {
        std::vector<Pt> p(n);
        for (auto& q : p) in >> q.x >> q.y;
        out.push_back(p);
    }
    return out;
}

static std::vector<Pt> perturb(std::mt19937_64& rng, std::vector<Pt> p) {
    const double scale = std::pow(10.0, 0.5 + 3.5 * urand(rng)); // jitter 3 .. 10^4
    i64 xmin = INT64_MAX, xmax = INT64_MIN, ymin = INT64_MAX, ymax = INT64_MIN;
    for (auto& q : p) {
        xmin = std::min(xmin, q.x); xmax = std::max(xmax, q.x);
        ymin = std::min(ymin, q.y); ymax = std::max(ymax, q.y);
    }
    const int mode = rng() % 4;
    if (mode != 3)
        for (auto& q : p) {
            q.x += (i64)llround((2 * urand(rng) - 1) * scale);
            q.y += (i64)llround((2 * urand(rng) - 1) * scale * (double)(ymax - ymin + 1) / (xmax - xmin + 1) * 4);
        }
    const int dels = rng() % 3, ins = rng() % 5;
    for (int t = 0; t < dels && p.size() > 5; ++t) p.erase(p.begin() + rng() % p.size());
    for (int t = 0; t < ins; ++t) {
        // a point near a random hull edge, pushed outwards
        const size_t e = rng() % p.size();
        const Pt a = p[e], b = p[(e + 1) % p.size()];
        const double u = urand(rng);
        const double nx = (double)(b.y - a.y), ny = -(double)(b.x - a.x);
        const double len = std::sqrt(nx * nx + ny * ny) + 1;
        const double out = urand(rng) * 0.02 * len;
        p.push_back({(i64)llround(a.x + u * (b.x - a.x) + nx / len * out),
                     (i64)llround(a.y + u * (b.y - a.y) + ny / len * out)});
    }
    return strict_hull(p);
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "usage: hunt exh|big|pert <instances> <seed> [...]\n";
        return 2;
    }
    const std::string mode = argv[1];
    const long count = atol(argv[2]);
    std::mt19937_64 rng(atoll(argv[3]));
    if (mode == "exh") {
        const int nmax = argc > 4 ? atoi(argv[4]) : 12, kmax = argc > 5 ? atoi(argv[5]) : 6;
        for (long it = 0; it < count; ++it) {
            Instance I(random_convex(rng, it % NUM_FAMILIES, 5 + rng() % (nmax - 4)));
            if (I.n >= 5) test_exhaustive(I, kmax);
        }
    } else if (mode == "tiny") {
        // thin two-chain polygons with small integer coordinates
        const int W = argc > 4 ? atoi(argv[4]) : 60, kmax = argc > 5 ? atoi(argv[5]) : 6;
        for (long it = 0; it < count; ++it) {
            Instance I(two_chain(rng, 5 + rng() % 7, W * (0.5 + urand(rng)), 1 + (W / 8.0) * urand(rng)));
            if (I.n >= 5) test_exhaustive(I, kmax);
        }
    } else if (mode == "big") {
        const int nmax = argc > 4 ? atoi(argv[4]) : 30, kmax = argc > 5 ? atoi(argv[5]) : 12;
        for (long it = 0; it < count; ++it) {
            // two-chain polygons are where the weaker rules fail: half of the instances
            const int fam = (it & 1) ? F_TWOCHAIN : (int)((it / 2) % NUM_FAMILIES);
            Instance I(random_convex(rng, fam, 13 + rng() % (nmax - 12)));
            if (I.n >= 8) test_big(I, kmax, rng);
        }
    } else if (mode == "pert") {
        const auto bases = read_bases(argv[4]);
        for (long it = 0; it < count; ++it) {
            Instance I(perturb(rng, bases[it % bases.size()]));
            if (I.n < 5) continue;
            if (I.n <= 13) test_exhaustive(I, 8);
            else test_big(I, 10, rng);
        }
    }
    for (auto& [key, s] : stats) {
        std::cout << "STAT " << mode << " " << key.first << " " << key.second << " " << s.tests << " "
                  << s.fails << "\n";
        if (s.fails) std::cout << "CEX " << mode << " " << key.first << " " << s.cex << "\n";
    }
    return 0;
}
