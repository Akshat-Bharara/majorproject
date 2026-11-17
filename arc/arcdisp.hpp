// Exact k-dispersion for points in convex position under the ARC PROPERTY.
//
// P = (p_0..p_{n-1}) strictly convex, counterclockwise.  For a radius r,
// Far_r(q) = { x != q : |qx| >= r }.  P has the r-arc property if every
// Far_r(q) is a set of cyclically consecutive points.  Under it, an r-far
// k-set exists iff for some q the greedy chain q, a(q), a(a(q)), ... (k
// points, a = first far point counterclockwise) does not wrap around and its
// last point is at distance >= r from q  (see NOTES.md for the proofs).
//
// All decisions use exact integer arithmetic on squared distances; r is
// always passed as r^2 and |qx| >= r is inclusive.
//
// The library is a template over (Coord, Dist):
//   Narrow = ArcDisp<int64_t, int64_t>      |coord| <= 2^29
//   Wide   = ArcDisp<int64_t, __int128>     |coord| <= 2^61
#pragma once

#include <algorithm>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <vector>

namespace arc {

using i64 = long long;
using i128 = __int128;

template <class C> struct Pt { C x, y; };

enum class Status { OK, NOT_APPLICABLE };

template <class Dist> struct Result {
    Status status = Status::OK;
    Dist r2 = -1;            // optimal squared distance (if OK)
    std::vector<int> S;      // witness: an optimal k-set
    Dist bad_r2 = -1;        // adaptive mode: radius at which the arc property failed
    int decisions = 0;       // calls to D
    long long candidates = 0;
};

template <class C, class Dist> class ArcDisp {
public:
    int n;
    std::vector<Pt<C>> p;
    std::vector<int> peak;   // peak[q] = offset t of a maximum of row q (unimodal P only)
    bool have_peaks = false;

    explicit ArcDisp(std::vector<Pt<C>> pts) : n((int)pts.size()), p(std::move(pts)) {}

    // ---- predicates -------------------------------------------------------
    inline Dist d(int i, int j) const {
        const Dist dx = (Dist)p[i].x - (Dist)p[j].x, dy = (Dist)p[i].y - (Dist)p[j].y;
        return dx * dx + dy * dy;
    }
    inline int idx(int q, int t) const { int v = q + t; return v >= n ? v - n : v; }
    inline Dist row(int q, int t) const { return d(q, idx(q, t)); } // |p_q p_{q+t}|^2, 1 <= t <= n-1

    // strictly convex, counterclockwise, winding number one
    bool strictly_convex() const {
        if (n < 3) return false;
        auto cross = [&](int a, int b, int c) {
            return (i128)(p[b].x - p[a].x) * (p[c].y - p[a].y) - (i128)(p[b].y - p[a].y) * (p[c].x - p[a].x);
        };
        int wraps = 0;
        for (int i = 0; i < n; ++i) {
            const int a = i, b = idx(i, 1), c = idx(i, 2);
            if (cross(a, b, c) <= 0) return false;
            const C e1x = p[b].x - p[a].x, e1y = p[b].y - p[a].y, e2x = p[c].x - p[b].x, e2y = p[c].y - p[b].y;
            const bool up1 = e1y > 0 || (e1y == 0 && e1x > 0), up2 = e2y > 0 || (e2y == 0 && e2x > 0);
            if (!up1 && up2) ++wraps;
        }
        return wraps == 1;
    }

    // ---- arc property ------------------------------------------------------
    // number of maximal runs of far points in the linear order q+1, ..., q-1
    int far_arcs(int q, Dist r2) const {
        int runs = 0;
        bool prev = false;
        for (int t = 1; t < n; ++t) {
            const bool f = row(q, t) >= r2;
            if (f && !prev) ++runs;
            prev = f;
        }
        return runs;
    }
    // O(n^2).  beta = sum_q max(0, #arcs(Far(q)) - 1);  arc property <=> beta = 0.
    bool arc_property(Dist r2, long long* beta = nullptr, int* arc_points = nullptr) const {
        long long b = 0;
        int good = 0;
        for (int q = 0; q < n; ++q) {
            const int a = far_arcs(q, r2);
            if (a <= 1) ++good;
            else b += a - 1;
        }
        if (beta) *beta = b;
        if (arc_points) *arc_points = good;
        return b == 0;
    }

    // ---- unimodality -------------------------------------------------------
    // Row q is unimodal if |p_q p_{q+t}|, t = 1..n-1, is strictly increasing
    // and then strictly decreasing, where the maximum may be attained twice
    // (at two consecutive t).  O(n^2) reference implementation.
    bool unimodal_naive() const {
        for (int q = 0; q < n; ++q) {
            int phase = 0; // 0 increasing, 1 after a tie at the top, 2 decreasing
            for (int t = 1; t + 1 < n; ++t) {
                const Dist a = row(q, t), b = row(q, t + 1);
                if (b > a) { if (phase != 0) return false; }
                else if (b == a) { if (phase != 0) return false; phase = 1; }
                else phase = 2;
            }
        }
        return true;
    }
    // O(n log n).  For the hull edge (p_j, p_{j+1}) the points strictly closer
    // to p_j form a clockwise run of c_j points starting at p_j (a line cuts a
    // convex polygon into two arcs); z_j = 1 if the next point is equidistant.
    // All rows are unimodal iff c_j > min(c_{j+1} + z_{j+1} - 2, n - 3) for all j.
    bool unimodal_fast() const {
        if (n < 4) return true;
        std::vector<int> c(n), z(n);
        for (int j = 0; j < n; ++j) {
            const int j1 = idx(j, 1);
            int lo = 1, hi = n - 1; // c_j in [1, n-1]: offsets t < c_j are "plus"
            while (lo < hi) {
                const int mid = (lo + hi + 1) / 2; // is offset mid-1 plus?
                const int q = idx(j, n - (mid - 1));
                if (mid - 1 == 0 || d(q, j) < d(q, j1)) lo = mid;
                else hi = mid - 1;
            }
            c[j] = lo;
            z[j] = 0;
            if (lo <= n - 2) {
                const int q = idx(j, n - lo);
                z[j] = d(q, j) == d(q, j1);
            }
        }
        for (int j = 0; j < n; ++j) {
            const int j1 = idx(j, 1);
            if (c[j] <= std::min(c[j1] + z[j1] - 2, n - 3)) return false;
        }
        return true;
    }
    // peaks of all rows by binary search; requires unimodal P.  O(n log n).
    void compute_peaks() {
        peak.assign(n, 1);
        for (int q = 0; q < n; ++q) {
            int lo = 1, hi = n - 1; // largest m such that row increases strictly on [1, m]
            while (lo < hi) {
                const int mid = (lo + hi + 1) / 2;
                if (mid == 1 || row(q, mid) > row(q, mid - 1)) lo = mid;
                else hi = mid - 1;
            }
            peak[q] = lo;
        }
        have_peaks = true;
    }

    // ---- first far points --------------------------------------------------
    // offsets t in 1..n-1 (ccw offset of a(q); ccw offset of b(q)); 0 if undefined
    int a_scan(int q, Dist r2) const {
        for (int t = 1; t < n; ++t) if (row(q, t) >= r2) return t;
        return 0;
    }
    int b_scan(int q, Dist r2) const {
        for (int t = n - 1; t >= 1; --t) if (row(q, t) >= r2) return t;
        return 0;
    }
    int a_uni(int q, Dist r2) const {
        const int m = peak[q];
        if (row(q, m) < r2) return 0;
        int lo = 1, hi = m;
        while (lo < hi) {
            const int mid = (lo + hi) / 2;
            if (row(q, mid) >= r2) hi = mid; else lo = mid + 1;
        }
        return lo;
    }
    int b_uni(int q, Dist r2) const {
        const int m = peak[q];
        if (row(q, m) < r2) return 0;
        int lo = m, hi = n - 1;
        while (lo < hi) {
            const int mid = (lo + hi + 1) / 2;
            if (row(q, mid) >= r2) lo = mid; else hi = mid - 1;
        }
        return lo;
    }

    // ---- decision ----------------------------------------------------------
    // Chain test exactly as specified: O(k) per seed, O(nk) in total (plus the
    // cost of a(.), b(.): O(n^2) by scanning, O(n log n) for unimodal P).
    //   (i)  the chain does not wrap and every a(.) is defined;
    //   (ii) for every i the cyclic predecessor of x_i lies in Far(x_i) = [a(x_i), b(x_i)].
    // Valid under the r-arc property.
    bool decide_chain(Dist r2, int k, bool unimodal, std::vector<int>* witness = nullptr) const {
        if (k < 2 || k > n) return false;
        std::vector<int> ta(n), tb(n);
        for (int q = 0; q < n; ++q) {
            ta[q] = unimodal ? a_uni(q, r2) : a_scan(q, r2);
            tb[q] = unimodal ? b_uni(q, r2) : b_scan(q, r2);
        }
        std::vector<int> off(k);
        for (int q = 0; q < n; ++q) {
            bool ok = true;
            off[0] = 0;
            for (int i = 1; i < k && ok; ++i) {
                const int x = idx(q, off[i - 1]);
                if (!ta[x] || off[i - 1] + ta[x] >= n) ok = false;
                else off[i] = off[i - 1] + ta[x];
            }
            for (int i = 0; i < k && ok; ++i) {
                const int x = idx(q, off[i]);
                // ccw offset from x_i to its cyclic predecessor
                const int back = i ? n - (off[i] - off[i - 1]) : off[k - 1];
                if (back > tb[x]) ok = false;
            }
            if (ok) {
                if (witness) { witness->clear(); for (int i = 0; i < k; ++i) witness->push_back(idx(q, off[i])); }
                return true;
            }
        }
        return false;
    }
    // Same decision in O(n log k) after a(.): under the arc property, (ii) is
    // equivalent to |x_0 x_{k-1}| >= r.  a^{k-1} by repeated squaring.
    bool decide_lift(Dist r2, int k, bool unimodal, std::vector<int>* witness = nullptr) const {
        if (k < 2 || k > n) return false;
        std::vector<int> ta(n);
        for (int q = 0; q < n; ++q) ta[q] = unimodal ? a_uni(q, r2) : a_scan(q, r2);
        // functions q -> (offset travelled, clamped to n = "wrapped or undefined")
        std::vector<int> cur(n), res(n, 0), tmp(n);
        for (int q = 0; q < n; ++q) cur[q] = ta[q] ? ta[q] : n;
        auto compose = [&](const std::vector<int>& f, const std::vector<int>& g, std::vector<int>& out) {
            for (int q = 0; q < n; ++q) { // first f, then g
                const int o = f[q];
                if (o >= n) { out[q] = n; continue; }
                const int o2 = g[idx(q, o)];
                out[q] = (o2 >= n || o + o2 >= n) ? n : o + o2;
            }
        };
        for (int e = k - 1; e > 0; e >>= 1) {
            if (e & 1) { compose(res, cur, tmp); res.swap(tmp); }
            if (e > 1) { compose(cur, cur, tmp); cur.swap(tmp); }
        }
        for (int q = 0; q < n; ++q) {
            if (res[q] >= n) continue;
            if (d(q, idx(q, res[q])) < r2) continue;
            if (witness) {
                witness->assign(1, q);
                int x = q;
                for (int i = 1; i < k; ++i) { x = idx(x, ta[x]); witness->push_back(x); }
            }
            return true;
        }
        return false;
    }

    Dist min_dist(const std::vector<int>& S) const { // O(|S|^2), for tests
        Dist m = -1;
        for (size_t a = 0; a < S.size(); ++a)
            for (size_t b = a + 1; b < S.size(); ++b) {
                const Dist v = d(S[a], S[b]);
                if (m < 0 || v < m) m = v;
            }
        return m;
    }

    std::vector<Dist> candidates() const {
        std::vector<Dist> c;
        c.reserve((size_t)n * (n - 1) / 2);
        for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) c.push_back(d(i, j));
        std::sort(c.begin(), c.end());
        c.erase(std::unique(c.begin(), c.end()), c.end());
        return c;
    }

    // ---- optimisation ------------------------------------------------------
    // (O1) sort all distinct squared distances, binary search with D.
    //   unimodal = true : P must be unimodal (fast a(.)); always applicable.
    //   check_arc = true: ADAPTIVE mode for general convex P: the arc property
    //                     is verified at every tested radius (O(n^2) each);
    //                     if it fails, NOT_APPLICABLE is returned.
    Result<Dist> solve_O1(int k, bool unimodal, bool check_arc) {
        Result<Dist> R;
        if (unimodal && !have_peaks) compute_peaks();
        const std::vector<Dist> cand = candidates();
        R.candidates = (long long)cand.size();
        size_t lo = 0, hi = cand.size() - 1; // cand[lo] feasible: every k-set is cand[0]-far
        auto D = [&](Dist r2, std::vector<int>* w) {
            ++R.decisions;
            return decide_lift(r2, k, unimodal, w);
        };
        if (check_arc && !arc_property(cand[0])) { R.status = Status::NOT_APPLICABLE; R.bad_r2 = cand[0]; return R; }
        while (lo < hi) {
            const size_t mid = lo + (hi - lo + 1) / 2;
            if (check_arc && !arc_property(cand[mid])) {
                R.status = Status::NOT_APPLICABLE;
                R.bad_r2 = cand[mid];
                return R;
            }
            if (D(cand[mid], nullptr)) lo = mid; else hi = mid - 1;
        }
        R.r2 = cand[lo];
        if (!D(R.r2, &R.S)) throw std::logic_error("O1: no witness at the optimum");
        std::sort(R.S.begin(), R.S.end());
        return R;
    }
    Result<Dist> solve_adaptive(int k) { return solve_O1(k, false, true); }

    // (O2) unimodal P.  Every row splits at its peak into an increasing and a
    // decreasing list; these 2n sorted lists hold every candidate (twice).
    // Randomised binary search: keep (lo, hi) with lo feasible and hi not,
    // pick a uniformly random candidate strictly between and test it.
    // Expected O(log n) rounds of O(n log n) each.
    Result<Dist> solve_O2(int k, uint64_t seed = 1) {
        Result<Dist> R;
        if (!have_peaks) compute_peaks();
        std::mt19937_64 rng(seed);
        // per row q: increasing list t in [1, m], decreasing list t in [m+1, n-1].
        // Active ranges (values strictly between lo and hi):
        //   increasing: t in [iL, iR)     decreasing: t in [dL, dR)
        std::vector<int> iL(n), iR(n), dL(n), dR(n);
        for (int q = 0; q < n; ++q) { iL[q] = 1; iR[q] = peak[q] + 1; dL[q] = peak[q] + 1; dR[q] = n; }
        Dist lo = 0; // every k-set is 0-far
        auto first_inc = [&](int q, int a, int b, Dist v, bool strict) { // first t in [a,b) with row > v (strict) / >= v
            while (a < b) {
                const int mid = (a + b) / 2;
                const Dist x = row(q, mid);
                if (strict ? x > v : x >= v) b = mid; else a = mid + 1;
            }
            return a;
        };
        auto first_dec = [&](int q, int a, int b, Dist v, bool strict) { // first t in [a,b) with row < v (strict) / <= v
            while (a < b) {
                const int mid = (a + b) / 2;
                const Dist x = row(q, mid);
                if (strict ? x < v : x <= v) b = mid; else a = mid + 1;
            }
            return a;
        };
        for (;;) {
            long long total = 0;
            for (int q = 0; q < n; ++q) total += (iR[q] - iL[q]) + (dR[q] - dL[q]);
            if (total == 0) break;
            long long pick = (long long)(rng() % (unsigned long long)total);
            Dist pivot = 0;
            for (int q = 0; q < n; ++q) {
                const int ci = iR[q] - iL[q], cd = dR[q] - dL[q];
                if (pick < ci) { pivot = row(q, iL[q] + (int)pick); break; }
                pick -= ci;
                if (pick < cd) { pivot = row(q, dL[q] + (int)pick); break; }
                pick -= cd;
            }
            ++R.decisions;
            if (decide_lift(pivot, k, true, nullptr)) {
                lo = pivot; // keep values > pivot
                for (int q = 0; q < n; ++q) {
                    iL[q] = first_inc(q, iL[q], iR[q], pivot, true);
                    dR[q] = first_dec(q, dL[q], dR[q], pivot, false);
                }
            } else { // keep values < pivot
                for (int q = 0; q < n; ++q) {
                    iR[q] = first_inc(q, iL[q], iR[q], pivot, false);
                    dL[q] = first_dec(q, dL[q], dR[q], pivot, true);
                }
            }
        }
        R.r2 = lo;
        if (!decide_lift(lo, k, true, &R.S)) throw std::logic_error("O2: no witness at the optimum");
        std::sort(R.S.begin(), R.S.end());
        return R;
    }
};

using Narrow = ArcDisp<i64, i64>;
using Wide = ArcDisp<i64, i128>;

} // namespace arc
