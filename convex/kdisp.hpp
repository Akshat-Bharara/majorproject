// Exact k-dispersion for points in convex position.
//
// Header-only library.  All geometry is done with exact integer predicates:
// input coordinates must satisfy |x|,|y| <= 2^29, so squared distances fit in
// 62 bits and the in-circle determinant fits in a signed 128-bit integer.
//
// Contents
//   Instance            points in strictly convex position, ccw order
//   CappedDP            the capped Delaunay dynamic program (decision)
//   TripleDP            straightforward O(n^4) apex-state baseline (decision)
//   brute_force         exhaustive branch-and-bound (reference)
//   farthest_first      greedy 2-approximation, used for the search window
//   seeded greedies     sound heuristics used by the certified algorithm
//   solve_capped / solve_certified / solve_triple   optimisation drivers
//   max_far_set         maximum r-far subset (maximum independent set)
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <istream>
#include <stdexcept>
#include <string>
#include <vector>

namespace kd {

using i64 = long long;
using i128 = __int128;
using u64 = unsigned long long;

constexpr i64 COORD_LIMIT = 1LL << 29;

struct Pt {
    i64 x = 0, y = 0;
};

inline i64 dist2(const Pt& a, const Pt& b) {
    const i64 dx = a.x - b.x, dy = a.y - b.y;
    return dx * dx + dy * dy;
}

// > 0 iff a,b,c make a left (counterclockwise) turn.
inline i128 orient(const Pt& a, const Pt& b, const Pt& c) {
    return (i128)(b.x - a.x) * (c.y - a.y) - (i128)(b.y - a.y) * (c.x - a.x);
}

// For a,b,c in counterclockwise order: > 0 iff d is strictly inside the
// circle through a,b,c;  = 0 iff cocircular;  < 0 iff strictly outside.
// |entries| <= 2^30 and 2^61, so |det| < 2^124.
inline int incircle(const Pt& a, const Pt& b, const Pt& c, const Pt& d) {
    const i64 adx = a.x - d.x, ady = a.y - d.y;
    const i64 bdx = b.x - d.x, bdy = b.y - d.y;
    const i64 cdx = c.x - d.x, cdy = c.y - d.y;
    const i128 ad = (i128)adx * adx + (i128)ady * ady;
    const i128 bd = (i128)bdx * bdx + (i128)bdy * bdy;
    const i128 cd = (i128)cdx * cdx + (i128)cdy * cdy;
    const i128 det = adx * (bdy * cd - cdy * bd) - ady * (bdx * cd - cdx * bd) +
                     ad * ((i128)bdx * cdy - (i128)cdx * bdy);
    return (det > 0) - (det < 0);
}

// Strict convex hull (Andrew's monotone chain), counterclockwise, no
// collinear points kept.  Duplicates are removed.
inline std::vector<Pt> strict_hull(std::vector<Pt> p) {
    std::sort(p.begin(), p.end(), [](const Pt& a, const Pt& b) {
        return a.x != b.x ? a.x < b.x : a.y < b.y;
    });
    p.erase(std::unique(p.begin(), p.end(),
                        [](const Pt& a, const Pt& b) { return a.x == b.x && a.y == b.y; }),
            p.end());
    const int n = (int)p.size();
    if (n < 3) return p;
    std::vector<Pt> h(2 * n);
    int m = 0;
    for (int i = 0; i < n; ++i) {
        while (m >= 2 && orient(h[m - 2], h[m - 1], p[i]) <= 0) --m;
        h[m++] = p[i];
    }
    for (int i = n - 2, t = m + 1; i >= 0; --i) {
        while (m >= t && orient(h[m - 2], h[m - 1], p[i]) <= 0) --m;
        h[m++] = p[i];
    }
    h.resize(m - 1);
    return h;
}

struct Instance {
    int n = 0;
    std::vector<Pt> p;
    std::vector<i64> D; // squared distances, row-major n x n

    Instance() = default;
    explicit Instance(std::vector<Pt> pts) : n((int)pts.size()), p(std::move(pts)) { build(); }

    i64 d(int i, int j) const { return D[(size_t)i * n + j]; }

    void build() {
        n = (int)p.size();
        D.assign((size_t)n * n, 0);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) D[(size_t)i * n + j] = dist2(p[i], p[j]);
    }

    // Throws unless the points are within range and in strictly convex
    // position, listed counterclockwise.
    void validate() const {
        for (const Pt& q : p)
            if (std::llabs(q.x) > COORD_LIMIT || std::llabs(q.y) > COORD_LIMIT)
                throw std::invalid_argument("coordinate out of range (|x|,|y| <= 2^29)");
        if (n < 3) return;
        // Every turn is a strict left turn and the turning number is one.
        int wraps = 0;
        for (int i = 0; i < n; ++i) {
            const Pt &a = p[i], &b = p[(i + 1) % n], &c = p[(i + 2) % n];
            if (orient(a, b, c) <= 0)
                throw std::invalid_argument("points are not in strictly convex ccw order");
            const Pt e1{b.x - a.x, b.y - a.y}, e2{c.x - b.x, c.y - b.y};
            const bool up1 = e1.y > 0 || (e1.y == 0 && e1.x > 0);
            const bool up2 = e2.y > 0 || (e2.y == 0 && e2.x > 0);
            if (!up1 && up2) ++wraps;
        }
        if (wraps != 1) throw std::invalid_argument("polygon winds more than once");
    }

    long double perimeter() const {
        long double s = 0;
        for (int i = 0; i < n; ++i) s += std::sqrt((long double)d(i, (i + 1) % n));
        return s;
    }

    std::vector<i64> candidates() const {
        std::vector<i64> c;
        c.reserve((size_t)n * (n - 1) / 2);
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) c.push_back(d(i, j));
        std::sort(c.begin(), c.end());
        c.erase(std::unique(c.begin(), c.end()), c.end());
        return c;
    }
};

// Reads "n k" followed by n lines "x y".  Points are used as given.
inline bool read_instance(std::istream& in, Instance& I, int& k) {
    int n;
    if (!(in >> n >> k)) return false;
    std::vector<Pt> p(n);
    for (auto& q : p)
        if (!(in >> q.x >> q.y)) return false;
    I = Instance(std::move(p));
    return true;
}

inline i64 min_dist2(const Instance& I, const std::vector<int>& S) {
    i64 m = INT64_MAX;
    for (size_t a = 0; a < S.size(); ++a)
        for (size_t b = a + 1; b < S.size(); ++b) m = std::min(m, I.d(S[a], S[b]));
    return m;
}

// ---------------------------------------------------------------------------
// Capped Delaunay dynamic program.
//
// Chords are ordered pairs (i,j), 0 <= i < n, i < j < i+n (positions >= n are
// taken modulo n).  In linear mode only chords with j <= n-1 are processed,
// which suffices for the decision problem; cyclic mode processes all n(n-1)
// ordered chords and is used by the tight-pair experiments.
//
// For every admissible chord we store L_1..L_c: L_t is an apex of smallest
// apex angle among the apexes whose (capped) value is at least t.
// ---------------------------------------------------------------------------
struct CappedDP {
    const Instance& I;
    int n;
    int K = 0;       // cap = k - 2
    i64 r2 = 0;
    bool cyclic = false;
    std::vector<size_t> off; // chord -> offset into L
    std::vector<int> cnt;    // chord -> number of valid thresholds
    std::vector<int> L;      // apex positions (in [i+1, j-1])
    std::vector<u64> farbits; // linear mode: bit rows of the far relation
    int W = 0;
    // statistics
    u64 n_incircle = 0, n_triples = 0, n_chords = 0;
    int found_i = -1, found_j = -1;

    explicit CappedDP(const Instance& inst) : I(inst), n(inst.n) {}

    size_t cid(int i, int j) const { return (size_t)i * n + (j - i); }
    bool far(int i, int j) const { return I.d(i % n, j % n) >= r2; }
    const Pt& P(int i) const { return I.p[i >= n ? i - n : i]; }

    int count(int i, int j) const {
        if (i >= n) { i -= n; j -= n; }
        return cnt[cid(i, j)];
    }
    const int* thr(int i, int j) const {
        if (i >= n) { i -= n; j -= n; }
        return L.data() + off[cid(i, j)];
    }

    // number of t with L_t(i,l) not strictly inside the circle (a,b,c)
    inline int query(int ci, int cj, const Pt& a, const Pt& b, const Pt& c) {
        int shift = 0;
        if (ci >= n) { ci -= n; cj -= n; shift = n; }
        const size_t id = cid(ci, cj);
        int hi = cnt[id];
        if (hi == 0) return 0;
        const int* t = L.data() + off[id];
        (void)shift;
        // largest q in [0,hi] with ok(q); ok is monotone (true then false)
        int lo = 0;
        while (lo < hi) {
            const int mid = (lo + hi + 1) >> 1;
            ++n_incircle;
            if (incircle(a, b, c, P(t[mid - 1])) <= 0) lo = mid;
            else hi = mid - 1;
        }
        return lo;
    }

    // Returns true iff an r-far k-set exists (k >= 3).  With early = false
    // the whole table is filled (needed for cyclic queries / max_value()).
    bool run(i64 r2_, int k, bool cyclic_ = false, bool early = true) {
        r2 = r2_;
        K = k - 2;
        cyclic = cyclic_;
        found_i = found_j = -1;
        n_incircle = n_triples = n_chords = 0;
        if (K < 1) throw std::invalid_argument("CappedDP needs k >= 3");
        off.assign((size_t)n * n + 1, 0);
        cnt.assign((size_t)n * n, 0);
        size_t total = 0;
        for (int i = 0; i < n; ++i)
            for (int g = 1; g < n; ++g) {
                off[cid(i, i + g)] = total;
                if (g >= 2 && (cyclic || i + g < n)) total += (size_t)std::min(K, g - 1);
            }
        L.assign(total, -1);

        if (!cyclic) {
            W = (n + 63) / 64;
            farbits.assign((size_t)n * W, 0);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    if (i != j && I.d(i, j) >= r2) farbits[(size_t)i * W + (j >> 6)] |= 1ULL << (j & 63);
        }

        std::vector<int> best(K + 1);
        bool found = false;
        for (int g = 2; g < n; ++g) {
            for (int i = 0; i < n; ++i) {
                const int j = i + g;
                if (!cyclic && j >= n) break;
                if (!far(i, j)) continue;
                ++n_chords;
                const int cap = std::min(K, g - 1);
                std::fill(best.begin(), best.begin() + cap + 1, -1);
                const Pt &pi = P(i), &pj = P(j);
                auto apex = [&](int l) {
                    ++n_triples;
                    const Pt& pl = P(l);
                    // (pi, pl, pj) is counterclockwise
                    int v = 1 + query(i, l, pi, pl, pj);
                    if (v < cap) v += query(l, j, pi, pl, pj);
                    if (v > cap) v = cap;
                    int& b = best[v];
                    if (b < 0) b = l;
                    else {
                        ++n_incircle;
                        // l has the smaller apex angle iff b is strictly inside circle(i,l,j)
                        if (incircle(pi, pl, pj, P(b)) > 0) b = l;
                    }
                };
                if (!cyclic) {
                    const u64* ri = &farbits[(size_t)i * W];
                    const u64* rj = &farbits[(size_t)j * W];
                    const int lo = i + 1, hi = j - 1;
                    for (int w = lo >> 6; w <= (hi >> 6); ++w) {
                        u64 m = ri[w] & rj[w];
                        if (w == (lo >> 6)) m &= ~0ULL << (lo & 63);
                        if (w == (hi >> 6) && (hi & 63) != 63) m &= (1ULL << ((hi & 63) + 1)) - 1;
                        while (m) {
                            apex((w << 6) + __builtin_ctzll(m));
                            m &= m - 1;
                        }
                    }
                } else {
                    for (int l = i + 1; l < j; ++l)
                        if (far(i, l) && far(l, j)) apex(l);
                }
                // suffix minima of the apex angle
                int* out = L.data() + off[cid(i, j)];
                int cur = -1, c = 0;
                for (int v = cap; v >= 1; --v) {
                    const int b = best[v];
                    if (b >= 0) {
                        if (cur < 0) cur = b;
                        else {
                            ++n_incircle;
                            if (incircle(pi, P(b), pj, P(cur)) > 0) cur = b;
                        }
                    }
                    out[v - 1] = cur;
                    if (cur >= 0 && c == 0) c = v;
                }
                cnt[cid(i, j)] = c;
                if (c >= K && !found) {
                    found = true;
                    found_i = i;
                    found_j = j;
                    if (early) return true;
                }
            }
        }
        return found;
    }

    // Largest (capped) value over all processed chords: min(K, F(i,j | pi)).
    int max_value() const {
        int m = 0;
        for (int v : cnt) m = std::max(m, v);
        return m;
    }

    // Reconstruct a witness: chord (i,j) together with `demand` interior points.
    void collect(int i, int j, int demand, std::vector<int>& out) {
        if (demand <= 0) return;
        if (i >= n) { i -= n; j -= n; }
        const int l = thr(i, j)[demand - 1];
        out.push_back(l % n);
        const Pt &pi = P(i), &pl = P(l), &pj = P(j);
        const int q1 = query(i, l, pi, pl, pj);
        const int c1 = std::min(q1, demand - 1);
        collect(i, l, c1, out);
        collect(l, j, demand - 1 - c1, out);
    }

    std::vector<int> witness(int i, int j, int demand) {
        std::vector<int> S{i % n, j % n};
        collect(i, j, demand, S);
        std::sort(S.begin(), S.end());
        return S;
    }
    std::vector<int> witness() { return witness(found_i, found_j, K); }
};

// ---------------------------------------------------------------------------
// Baseline: the direct dynamic program whose state is a chord together with
// its apex.  f(i,j,l) = number of interior points of the best configuration
// of chord (i,j) with apex l.  O(n^4) in-circle tests per decision, O(n^3)
// space.  No thresholds, no capping.
// ---------------------------------------------------------------------------
struct TripleDP {
    const Instance& I;
    int n;
    std::vector<int16_t> f;
    u64 n_incircle = 0;
    explicit TripleDP(const Instance& inst) : I(inst), n(inst.n) {}

    // maximum size of an r-far subset if it is >= 3, otherwise 0
    int max_size(i64 r2) {
        f.assign((size_t)n * n * n, 0);
        n_incircle = 0;
        auto at = [&](int i, int j, int l) -> int16_t& { return f[((size_t)i * n + j) * n + l]; };
        int best = 0;
        for (int g = 2; g < n; ++g)
            for (int i = 0; i + g < n; ++i) {
                const int j = i + g;
                if (I.d(i, j) < r2) continue;
                for (int l = i + 1; l < j; ++l) {
                    if (I.d(i, l) < r2 || I.d(l, j) < r2) continue;
                    int b1 = 0, b2 = 0;
                    for (int x = i + 1; x < l; ++x) {
                        const int v = at(i, l, x);
                        if (v > b1) {
                            ++n_incircle;
                            if (incircle(I.p[i], I.p[l], I.p[j], I.p[x]) <= 0) b1 = v;
                        }
                    }
                    for (int x = l + 1; x < j; ++x) {
                        const int v = at(l, j, x);
                        if (v > b2) {
                            ++n_incircle;
                            if (incircle(I.p[i], I.p[l], I.p[j], I.p[x]) <= 0) b2 = v;
                        }
                    }
                    const int v = 1 + b1 + b2;
                    at(i, j, l) = (int16_t)v;
                    best = std::max(best, v);
                }
            }
        return best ? best + 2 : 0;
    }
    bool decide(i64 r2, int k) { return max_size(r2) >= k; }
};

// ---------------------------------------------------------------------------
// Exhaustive reference: branch and bound over k-subsets.
// ---------------------------------------------------------------------------
struct Brute {
    const Instance& I;
    int k;
    i64 best = -1;
    std::vector<int> cur, arg;
    Brute(const Instance& inst, int k_) : I(inst), k(k_) {}
    void rec(int next, i64 curmin) {
        if ((int)cur.size() == k) {
            if (curmin > best) { best = curmin; arg = cur; }
            return;
        }
        for (int v = next; v <= I.n - (k - (int)cur.size()); ++v) {
            i64 m = curmin;
            for (int u : cur) m = std::min(m, I.d(u, v));
            if (m <= best) continue;
            cur.push_back(v);
            rec(v + 1, m);
            cur.pop_back();
        }
    }
    i64 solve() {
        rec(0, INT64_MAX);
        return best;
    }
};

// ---------------------------------------------------------------------------
// Farthest-first greedy (Gonzalez).  Returns the k chosen points; the
// minimum pairwise distance g of the result satisfies g <= r* <= 2g.
// ---------------------------------------------------------------------------
inline std::vector<int> farthest_first(const Instance& I, int k, int start = 0) {
    const int n = I.n;
    std::vector<i64> near(n, INT64_MAX);
    std::vector<int> S{start};
    for (int t = 1; t < k; ++t) {
        const int last = S.back();
        int arg = -1;
        for (int v = 0; v < n; ++v) {
            near[v] = std::min(near[v], I.d(last, v));
            if (arg < 0 || near[v] > near[arg]) arg = v;
        }
        S.push_back(arg);
    }
    std::sort(S.begin(), S.end());
    return S;
}

// first point counterclockwise (dir=+1) or clockwise (dir=-1) of q at
// distance >= r, as an offset t in 1..n-1 (q + dir*t), or 0 if none exists.
inline int first_far(const Instance& I, i64 r2, int q, int dir) {
    const int n = I.n;
    for (int t = 1; t < n; ++t)
        if (I.d(q, ((q + dir * t) % n + n) % n) >= r2) return t;
    return 0;
}

// ---------------------------------------------------------------------------
// Seeded greedies.  Positions are offsets from the seed point q (offset 0 and
// offset n both denote q).  The free region is the open offset interval
// (lo, hi); a ccw step adds the first free point after lo that is far from
// all chosen points, a cw step the first one before hi.
// ---------------------------------------------------------------------------
struct Seeded {
    const Instance& I;
    i64 r2;
    int k;
    int n;
    Seeded(const Instance& inst, i64 r2_, int k_) : I(inst), r2(r2_), k(k_), n(inst.n) {}

    int pt(int q, int t) const { return (q + t) % n; }
    bool ok(int v, const std::vector<int>& S) const {
        for (int u : S)
            if (I.d(u, v) < r2) return false;
        return true;
    }
    // seed3 = true: {q, a(q), b(q)}, otherwise {q, a(q)}.
    bool seed(int q, bool seed3, std::vector<int>& S, int& lo, int& hi) const {
        const int ta = first_far(I, r2, q, +1);
        if (!ta) return false;
        S = {q, pt(q, ta)};
        lo = ta;
        hi = n;
        if (seed3) {
            const int tb = first_far(I, r2, q, -1);
            if (n - tb <= ta) return false; // a(q) = b(q): fewer than three points
            if (I.d(pt(q, ta), pt(q, n - tb)) < r2) return false;
            S.push_back(pt(q, n - tb));
            hi = n - tb;
        }
        return true;
    }
    bool step(int q, bool ccw, std::vector<int>& S, int& lo, int& hi) const {
        if (ccw) {
            for (int t = lo + 1; t < hi; ++t)
                if (ok(pt(q, t), S)) { S.push_back(pt(q, t)); lo = t; return true; }
        } else {
            for (int t = hi - 1; t > lo; --t)
                if (ok(pt(q, t), S)) { S.push_back(pt(q, t)); hi = t; return true; }
        }
        return false;
    }
    // strictly alternating greedy; first_ccw selects the starting side
    bool alternating(int q, bool seed3, bool first_ccw, std::vector<int>* out = nullptr) const {
        std::vector<int> S;
        int lo, hi;
        if (!seed(q, seed3, S, lo, hi)) return false;
        bool side = first_ccw;
        while ((int)S.size() < k) {
            if (!step(q, side, S, lo, hi)) return false;
            side = !side;
        }
        if (out) *out = S;
        return true;
    }
    // search tree: every ccw/cw pattern
    bool tree_rec(int q, std::vector<int>& S, int lo, int hi) const {
        if ((int)S.size() >= k) return true;
        for (int side = 0; side < 2; ++side) {
            int l2 = lo, h2 = hi;
            const size_t sz = S.size();
            if (step(q, side == 0, S, l2, h2)) {
                if (tree_rec(q, S, l2, h2)) return true;
                S.resize(sz);
            }
        }
        return false;
    }
    bool tree(int q, bool seed3, std::vector<int>* out = nullptr) const {
        std::vector<int> S;
        int lo, hi;
        if (!seed(q, seed3, S, lo, hi)) return false;
        if (!tree_rec(q, S, lo, hi)) return false;
        if (out) *out = S;
        return true;
    }
    // The predicate used by the certified algorithm: alternating greedy from
    // both seed types and both starting sides, over all q.
    bool any_alternating(std::vector<int>* out = nullptr) const {
        for (int q = 0; q < n; ++q)
            for (int s3 = 0; s3 < 2; ++s3)
                for (int fc = 0; fc < 2; ++fc)
                    if (alternating(q, s3, fc, out)) return true;
        return false;
    }
};

// ---------------------------------------------------------------------------
// Optimisation drivers.
// ---------------------------------------------------------------------------
struct Result {
    i64 r2 = -1;             // optimal squared distance
    std::vector<int> S;      // an optimal k-set (sorted indices)
    int decisions = 0;       // exact decision calls
    int greedy_calls = 0;    // seeded-greedy predicate calls (certified only)
    bool certified_by_greedy = false;
    size_t cand_total = 0, cand_window = 0;
    u64 incircle = 0, triples = 0;
};

inline Result solve_k2(const Instance& I) {
    Result R;
    for (int i = 0; i < I.n; ++i)
        for (int j = i + 1; j < I.n; ++j)
            if (I.d(i, j) > R.r2) { R.r2 = I.d(i, j); R.S = {i, j}; }
    return R;
}

// Index range [lo, hi] into the sorted candidate list C that is guaranteed to
// contain r*^2, together with the greedy set attaining C[lo].
inline void search_window(const Instance& I, int k, const std::vector<i64>& C, bool use_window,
                          size_t& lo, size_t& hi, std::vector<int>& G) {
    G = farthest_first(I, k);
    const i64 g2 = min_dist2(I, G);
    lo = std::lower_bound(C.begin(), C.end(), g2) - C.begin();
    hi = C.size() - 1;
    if (!use_window) { lo = 0; return; }
    // r*^2 <= 4 g^2 and r* <= per/k (inflated slightly: any valid upper bound works)
    long double ub = 4.0L * (long double)g2;
    const long double pk = I.perimeter() / k * (1.0L + 1e-12L) + 1e-6L;
    ub = std::min(ub, pk * pk);
    while (hi > lo && (long double)C[hi] > ub) --hi;
}

// Binary search over candidate distances with the capped DP.
inline Result solve_capped(const Instance& I, int k, bool use_window = true, int cap_k = -1) {
    if (k == 2) return solve_k2(I);
    Result R;
    const std::vector<i64> C = I.candidates();
    size_t lo, hi;
    std::vector<int> G;
    search_window(I, k, C, use_window, lo, hi, G);
    R.cand_total = C.size();
    R.cand_window = hi - lo + 1;
    CappedDP dp(I);
    const int kk = cap_k < 0 ? k : cap_k; // cap_k = n disables capping (ablation)
    auto decide = [&](i64 r2) {
        ++R.decisions;
        bool yes;
        if (kk == k) yes = dp.run(r2, k);
        else { dp.run(r2, kk, false, false); yes = dp.max_value() >= k - 2; }
        R.incircle += dp.n_incircle;
        R.triples += dp.n_triples;
        return yes;
    };
    // invariant: C[lo] feasible (trivially when lo = 0: any k-set), answer in [lo, hi]
    R.S = G;
    bool have = use_window;
    if (!use_window) { R.S.clear(); }
    while (lo < hi) {
        const size_t mid = lo + (hi - lo + 1) / 2;
        if (decide(C[mid])) {
            lo = mid;
            if (kk == k) { R.S = dp.witness(); have = true; }
        } else hi = mid - 1;
    }
    R.r2 = C[lo];
    if (!have || min_dist2(I, R.S) != R.r2) {
        // need a witness at the optimum (also covers the ablation mode)
        CappedDP w(I);
        w.run(R.r2, k);
        R.S = w.witness();
    }
    return R;
}

// Certified algorithm: binary search with the (sound) seeded greedy, one
// exact decision to certify, exact binary search above if that fails.
inline Result solve_certified(const Instance& I, int k) {
    if (k == 2) return solve_k2(I);
    Result R;
    const std::vector<i64> C = I.candidates();
    size_t lo, hi;
    std::vector<int> G;
    search_window(I, k, C, true, lo, hi, G);
    R.cand_total = C.size();
    R.cand_window = hi - lo + 1;
    R.S = G; // C[lo] is attained by the farthest-first set
    size_t a = lo, b = hi;
    while (a < b) {
        const size_t mid = a + (b - a + 1) / 2;
        std::vector<int> S;
        ++R.greedy_calls;
        if (Seeded(I, C[mid], k).any_alternating(&S)) { a = mid; R.S = S; }
        else b = mid - 1;
    }
    // a = index of rho; every success so far exhibited an r-far k-set
    CappedDP dp(I);
    auto decide = [&](i64 r2) {
        ++R.decisions;
        const bool yes = dp.run(r2, k);
        R.incircle += dp.n_incircle;
        R.triples += dp.n_triples;
        return yes;
    };
    if (a == hi || !decide(C[a + 1])) {
        R.certified_by_greedy = true;
        R.r2 = C[a];
        std::sort(R.S.begin(), R.S.end());
        return R;
    }
    lo = a + 1;
    R.S = dp.witness();
    while (lo < hi) {
        const size_t mid = lo + (hi - lo + 1) / 2;
        if (decide(C[mid])) { lo = mid; R.S = dp.witness(); }
        else hi = mid - 1;
    }
    R.r2 = C[lo];
    return R;
}

// Binary search (with window) using the O(n^4) baseline as decision.
inline Result solve_triple(const Instance& I, int k) {
    if (k == 2) return solve_k2(I);
    Result R;
    const std::vector<i64> C = I.candidates();
    size_t lo, hi;
    std::vector<int> G;
    search_window(I, k, C, true, lo, hi, G);
    TripleDP dp(I);
    while (lo < hi) {
        const size_t mid = lo + (hi - lo + 1) / 2;
        ++R.decisions;
        const bool yes = dp.decide(C[mid], k);
        R.incircle += dp.n_incircle;
        if (yes) lo = mid;
        else hi = mid - 1;
    }
    R.r2 = C[lo];
    return R;
}

// ---------------------------------------------------------------------------
// Maximum r-far subset (maximum independent set of the "closer than r"
// graph).  A greedy maximal r-far set M gives alpha <= 5|M|, so capping at
// 5|M| loses nothing.  Returns a maximum r-far set.
// ---------------------------------------------------------------------------
inline std::vector<int> max_far_set(const Instance& I, i64 r2, u64* incircle_tests = nullptr) {
    const int n = I.n;
    std::vector<int> M;
    for (int v = 0; v < n; ++v) {
        bool ok = true;
        for (int u : M)
            if (I.d(u, v) < r2) { ok = false; break; }
        if (ok) M.push_back(v);
    }
    if (n < 3) return M;
    const int cap = std::min(n, 5 * (int)M.size());
    if (cap < 3) return M;
    CappedDP dp(I);
    dp.run(r2, cap + 1, false, false); // K = cap - 1 > alpha - 2: never truncates
    if (incircle_tests) *incircle_tests = dp.n_incircle;
    int bi = -1, bj = -1, bv = 0;
    for (int i = 0; i < n; ++i)
        for (int j = i + 2; j < n; ++j)
            if (dp.cnt[dp.cid(i, j)] > bv) { bv = dp.cnt[dp.cid(i, j)]; bi = i; bj = j; }
    if (bv == 0) {
        // no r-far triple: the answer is a far pair if one exists, else a point
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                if (I.d(i, j) >= r2) return {i, j};
        return {0};
    }
    return dp.witness(bi, bj, bv);
}

} // namespace kd
