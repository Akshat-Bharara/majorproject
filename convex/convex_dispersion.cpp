#include "convex_dispersion.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

namespace cvx {

ConvexDispersion::ConvexDispersion(std::vector<Pt> pts, int k)
    : p_(std::move(pts)), n_((int)p_.size()), k_(k) {
    if (k_ < 2 || k_ > n_) throw std::invalid_argument("require 2 <= k <= n");
    for (const auto& q : p_)
        if (q.x > kMaxCoord || q.x < -kMaxCoord || q.y > kMaxCoord || q.y < -kMaxCoord)
            throw std::invalid_argument("coordinates must lie in [-2^29, 2^29]");
    for (int i = 0; i < n_ && n_ >= 3; ++i)
        if (orient(p_[i], p_[nxt(i, 1)], p_[nxt(i, 2)]) <= 0)
            throw std::invalid_argument("points are not in strictly convex ccw position");
    d2_.resize((std::size_t)n_ * n_);
    for (int a = 0; a < n_; ++a)
        for (int b = 0; b < n_; ++b) d2_[(std::size_t)a * n_ + b] = dist2(p_[a], p_[b]);
}

// ---------------------------------------------------------------------------
// Capped Delaunay DP (Algorithm 1 of the report).
//
// For an ordered chord (i,j) let P(i,j) be the points strictly between i and j
// in ccw order. For c = 1..K (K = k-2) apex[i,j][c] is a point l of P(i,j)
// minimising the apex angle <(i,l,j) among the points whose value
//   v(l) = 1 + H(i,l | circle(i,l,j)) + H(l,j | circle(i,l,j))
// is at least c, where H(a,b | C) counts the levels of chord (a,b) whose apex
// lies outside the open disk bounded by C. Angles are never computed: for x, y
// on the same side of (i,j), <(i,x,j) > <(i,y,j) iff x is strictly inside the
// circle through i, y, j.
// ---------------------------------------------------------------------------
bool ConvexDispersion::decide_capped(i64 r2, std::vector<int>* witness, bool tight_only) {
    ++stats_.decisions;
    const int n = n_;
    auto far = [&](int a, int b) { return D(a, b) >= r2; };

    if (k_ == 2) {
        for (int a = 0; a < n; ++a)
            for (int b = a + 1; b < n; ++b)
                if (far(a, b)) {
                    if (witness) *witness = {a, b};
                    return true;
                }
        return false;
    }

    const int K = k_ - 2;
    std::vector<std::uint16_t> cnt((std::size_t)n * n, 0);
    std::vector<int> apex((std::size_t)n * n * K, -1);
    auto A = [&](int a, int b, int c) -> int& {  // c in 1..K
        return apex[((std::size_t)a * n + b) * K + (c - 1)];
    };

    // Number of levels of chord (a,b) usable when the neighbouring triangle on
    // the other side of (a,b) has third vertex t (local Delaunay condition).
    auto levels = [&](int a, int b, int t) -> int {
        const int m = cnt[(std::size_t)a * n + b];
        int lo = 0, hi = m;  // predicate "apex c not strictly inside" holds for c <= answer
        while (lo < hi) {
            const int mid = (lo + hi + 1) / 2;
            if (!strictly_inside_circle(p_[a], p_[b], p_[t], p_[A(a, b, mid)])) lo = mid;
            else hi = mid - 1;
        }
        return lo;
    };

    std::vector<int> best(K + 1);
    for (int g = 2; g < n; ++g) {
        for (int i = 0; i < n; ++i) {
            const int j = nxt(i, g);
            if (!far(i, j)) continue;
            std::fill(best.begin(), best.end(), -1);
            for (int s = 1; s < g; ++s) {
                const int l = nxt(i, s);
                if (!far(i, l) || !far(l, j)) continue;
                ++stats_.inner_steps;
                const int q1 = levels(i, l, j);
                const int q2 = levels(l, j, i);
                const int v = std::min(K, 1 + q1 + q2);
                int& b = best[v];
                if (b < 0 || strictly_inside_circle(p_[i], p_[l], p_[j], p_[b])) b = l;
            }
            int cur = -1, top = 0;
            for (int v = K; v >= 1; --v) {
                if (best[v] >= 0 &&
                    (cur < 0 || strictly_inside_circle(p_[i], p_[best[v]], p_[j], p_[cur])))
                    cur = best[v];
                A(i, j, v) = cur;
                if (cur >= 0 && top == 0) top = v;
            }
            cnt[(std::size_t)i * n + j] = (std::uint16_t)top;
        }
    }

    for (int i = 0; i < n; ++i)
        for (int g = 2; g < n; ++g) {
            const int j = nxt(i, g);
            if (!far(i, j) || cnt[(std::size_t)i * n + j] < K) continue;
            if (tight_only) {  // i must be a(j): no point of (j, i) is far from j
                bool tight = true;
                for (int s = 1; s < n - g && tight; ++s) tight = !far(j, nxt(j, s));
                if (!tight) continue;
            }
            if (witness) {
                witness->clear();
                witness->push_back(i);
                witness->push_back(j);
                // Reconstruct exactly c points strictly between a and b.
                std::function<void(int, int, int)> rec = [&](int a, int b, int c) {
                    if (c == 0) return;
                    const int l = A(a, b, c);
                    witness->push_back(l);
                    const int q1 = levels(a, l, b);
                    const int c1 = std::min(q1, c - 1);
                    rec(a, l, c1);
                    rec(l, b, c - 1 - c1);
                };
                rec(i, j, K);
                std::sort(witness->begin(), witness->end());
            }
            return true;
        }
    return false;
}

// ---------------------------------------------------------------------------
// Tkachenko--Wang DP (STACS 2025, Lemma 13), direct O(n^4) evaluation.
// f(i,j,t): max number of points X of P(i,j) outside the open disk through
// i, j, t (t = -1: no disk) such that X + {i,j} is r-far.
// ---------------------------------------------------------------------------
bool ConvexDispersion::decide_tw(i64 r2, std::vector<int>* witness) {
    ++stats_.decisions;
    const int n = n_;
    auto far = [&](int a, int b) { return D(a, b) >= r2; };
    if (k_ == 2) return decide_capped(r2, witness);

    // f stored as f[(i*n + j)*(n+1) + (t+1)], t in -1..n-1.
    std::vector<std::int16_t> f((std::size_t)n * n * (n + 1), -1);
    auto F = [&](int a, int b, int t) -> std::int16_t& {
        return f[((std::size_t)a * n + b) * (n + 1) + (t + 1)];
    };
    auto eval = [&](int i, int j, int t, int* arg) -> int {
        const int g = gap(i, j);
        int bestv = 0;
        if (arg) *arg = -1;
        for (int s = 1; s < g; ++s) {
            const int l = nxt(i, s);
            if (!far(i, l) || !far(l, j)) continue;
            if (t >= 0 && strictly_inside_circle(p_[i], p_[j], p_[t], p_[l])) continue;
            ++stats_.inner_steps;
            const int v = 1 + F(i, l, j) + F(l, j, i);
            if (v > bestv) { bestv = v; if (arg) *arg = l; }
        }
        return bestv;
    };

    for (int g = 1; g < n; ++g)
        for (int i = 0; i < n; ++i) {
            const int j = nxt(i, g);
            if (!far(i, j)) continue;
            F(i, j, -1) = (std::int16_t)eval(i, j, -1, nullptr);
            // canonical triples: t on the other side, far from i and j
            for (int s = 1; s < n - g; ++s) {
                const int t = nxt(j, s);
                if (far(t, i) && far(t, j)) F(i, j, t) = (std::int16_t)eval(i, j, t, nullptr);
            }
        }

    for (int i = 0; i < n; ++i)
        for (int g = 1; g < n; ++g) {
            const int j = nxt(i, g);
            if (!far(i, j) || F(i, j, -1) + 2 < k_) continue;
            if (witness) {
                witness->clear();
                std::function<void(int, int, int)> rec = [&](int a, int b, int t) {
                    int l;
                    if (eval(a, b, t, &l) == 0 || l < 0) return;
                    witness->push_back(l);
                    rec(a, l, b);
                    rec(l, b, a);
                };
                witness->push_back(i);
                witness->push_back(j);
                rec(i, j, -1);
                std::sort(witness->begin(), witness->end());
                witness->resize(k_);  // any k of the >= k points are r-far
            }
            return true;
        }
    return false;
}

bool ConvexDispersion::decide_greedy(i64 r2, std::vector<int>* witness) {
    ++stats_.greedy_decisions;
    const int n = n_;
    auto far = [&](int a, int b) { return D(a, b) >= r2; };
    std::vector<int> chosen;
    chosen.reserve(k_);
    auto far_all = [&](int v) {
        for (int c : chosen) if (!far(v, c)) return false;
        return true;
    };

    for (int q = 0; q < n; ++q) {
        int a = -1, b = -1;
        for (int s = 1; s < n && a < 0; ++s) if (far(q, nxt(q, s))) a = nxt(q, s);
        if (a < 0) continue;
        for (int s = 1; s < n && b < 0; ++s) if (far(q, nxt(q, n - s))) b = nxt(q, n - s);
        // seed 0: {q, a(q), b(q)}; seed 1: {q, a(q)}. order: which side first.
        for (int seed = 0; seed < 2; ++seed)
            for (int first = 0; first < 2; ++first) {
                chosen.clear();
                chosen.push_back(q);
                chosen.push_back(a);
                int x = a, y = q;  // ccw frontier, cw frontier
                if (seed == 0 && b != a) {
                    if (!far(a, b)) continue;
                    chosen.push_back(b);
                    y = b;
                }
                // Scan pointers only move forward: a rejected point stays
                // rejected because the chosen set only grows.
                int px = nxt(x, 1), py = (y == 0 ? n - 1 : y - 1);
                int side = first;
                bool stuck[2] = {false, false};
                while ((int)chosen.size() < k_ && !(stuck[0] && stuck[1])) {
                    if (side == 0 && !stuck[0]) {
                        while (px != y && !far_all(px)) px = nxt(px, 1);
                        if (px == y) stuck[0] = true;
                        else { chosen.push_back(px); x = px; px = nxt(px, 1); }
                    } else if (side == 1 && !stuck[1]) {
                        while (py != x && !far_all(py)) py = (py == 0 ? n - 1 : py - 1);
                        if (py == x) stuck[1] = true;
                        else { chosen.push_back(py); y = py; py = (py == 0 ? n - 1 : py - 1); }
                    }
                    side ^= 1;
                }
                if ((int)chosen.size() >= k_) {
                    if (witness) {
                        *witness = chosen;
                        witness->resize(k_);
                        std::sort(witness->begin(), witness->end());
                    }
                    return true;
                }
            }
    }
    return false;
}

// Farthest pair of points given in strictly convex ccw order (rotating
// calipers, Shamos 1978). Returns the pair or {-1,-1} if fewer than 2 points.
static std::pair<int, int> farthest_pair(const std::vector<Pt>& p, const std::vector<int>& F) {
    const int m = (int)F.size();
    if (m < 2) return {-1, -1};
    if (m == 2) return {F[0], F[1]};
    auto area2 = [&](int a, int b, int c) { return orient(p[F[a]], p[F[b]], p[F[c]]); };
    std::pair<int, int> best{F[0], F[1]};
    i64 bd = dist2(p[F[0]], p[F[1]]);
    auto consider = [&](int a, int b) {
        const i64 v = dist2(p[F[a]], p[F[b]]);
        if (v > bd) { bd = v; best = {F[a], F[b]}; }
    };
    int j = 1;
    for (int i = 0; i < m; ++i) {
        const int ni = (i + 1) % m;
        while (area2(i, ni, (j + 1) % m) > area2(i, ni, j)) j = (j + 1) % m;
        consider(i, j);
        consider(ni, j);
        // with a parallel edge pair both j and j+1 are antipodal to (i, ni)
        consider(i, (j + 1) % m);
        consider(ni, (j + 1) % m);
    }
    return best;
}

bool ConvexDispersion::decide_tight4(i64 r2, std::vector<int>* witness) {
    if (k_ != 4) throw std::invalid_argument("TightK4 requires k = 4");
    ++stats_.decisions;
    const int n = n_;
    auto far = [&](int a, int b) { return D(a, b) >= r2; };
    std::vector<int> F;
    F.reserve(n);
    for (int q = 0; q < n; ++q) {
        for (int orient_ = 0; orient_ < 2; ++orient_) {
            // orient_ 0: pair (q, a(q)), remaining points ccw from a(q) to q.
            // orient_ 1: pair (b(q), q), remaining points ccw from q to b(q).
            const int step = orient_ == 0 ? 1 : n - 1;
            int partner = -1;
            for (int s = 1; s < n && partner < 0; ++s) {
                const int v = (q + step * s) % n;
                if (far(q, v)) partner = v;
            }
            if (partner < 0) break;
            const int from = orient_ == 0 ? partner : q;
            const int to = orient_ == 0 ? q : partner;
            F.clear();
            for (int v = nxt(from, 1); v != to; v = nxt(v, 1))
                if (far(v, q) && far(v, partner)) F.push_back(v);
            stats_.inner_steps += F.size();
            const auto fp = farthest_pair(p_, F);
            if (fp.first >= 0 && far(fp.first, fp.second)) {
                if (witness) {
                    *witness = {q, partner, fp.first, fp.second};
                    std::sort(witness->begin(), witness->end());
                }
                return true;
            }
        }
    }
    return false;
}

std::vector<int> ConvexDispersion::greedy() const {
    int a = 0, b = 1;
    for (int x = 0; x < n_; ++x)
        for (int y = x + 1; y < n_; ++y)
            if (D(x, y) > D(a, b)) { a = x; b = y; }
    std::vector<int> chosen{a};
    std::vector<i64> md(n_);
    for (int x = 0; x < n_; ++x) md[x] = D(a, x);
    while ((int)chosen.size() < k_) {
        int arg = -1;
        for (int x = 0; x < n_; ++x)
            if (md[x] > 0 && (arg < 0 || md[x] > md[arg])) arg = x;
        if (arg < 0) break;
        chosen.push_back(arg);
        for (int x = 0; x < n_; ++x) md[x] = std::min(md[x], D(arg, x));
    }
    return chosen;
}

ConvexResult ConvexDispersion::solve(Algo algo, bool use_bounds) {
    stats_ = {};
    std::vector<i64> cand(d2_.begin(), d2_.end());
    cand.erase(std::remove(cand.begin(), cand.end(), (i64)0), cand.end());
    std::sort(cand.begin(), cand.end());
    cand.erase(std::unique(cand.begin(), cand.end()), cand.end());
    stats_.candidates_total = cand.size();

    // Window [L, U] (Lemma 6 of the report).
    i64 L2 = cand.front(), U2 = cand.back();
    std::vector<int> gsel = greedy();
    if (use_bounds && (int)gsel.size() == k_) {
        i64 g2 = std::numeric_limits<i64>::max();
        for (std::size_t a = 0; a < gsel.size(); ++a)
            for (std::size_t b = a + 1; b < gsel.size(); ++b)
                g2 = std::min(g2, D(gsel[a], gsel[b]));
        long double per = 0;
        for (int i = 0; i < n_; ++i) per += std::sqrt((long double)D(i, nxt(i, 1)));
        const long double pk = per / k_;
        // Conservative rounding up keeps the true optimum inside the window.
        const i64 perim2 = (i64)std::ceil(pk * pk * (1 + 1e-12L)) + 1;
        L2 = g2;
        U2 = std::min({U2, 4 * g2, perim2});
    }
    stats_.lower_bound2 = L2;
    stats_.upper_bound2 = U2;
    std::vector<i64> win;
    for (i64 c : cand) if (c >= L2 && c <= U2) win.push_back(c);
    stats_.candidates_window = win.size();

    auto decide = [&](i64 r2, std::vector<int>* w) {
        switch (algo) {
            case Algo::TWNaive: return decide_tw(r2, w);
            case Algo::TightK4: return decide_tight4(r2, w);
            default: return decide_capped(r2, w);
        }
    };

    // Largest feasible candidate; win.front() (= L2 when bounds are used, or
    // the smallest distance otherwise) is always feasible.
    int lo = 0, hi = (int)win.size() - 1;
    if (algo == Algo::Certified) {
        // Phase 1: binary search with the sound greedy; every success is a
        // valid lower bound, so lo stays feasible.
        int glo = 0, ghi = hi;
        while (glo < ghi) {
            const int mid = (glo + ghi + 1) / 2;
            if (decide_greedy(win[mid], nullptr)) glo = mid; else ghi = mid - 1;
        }
        lo = glo;
        // Phase 2: one exact decision just above the greedy value.
        if (lo == hi || !decide_capped(win[lo + 1], nullptr)) {
            stats_.certified_by_greedy = true;
            hi = lo;
        } else {
            lo = lo + 1;
        }
        algo = Algo::CappedDP;
    }
    while (lo < hi) {
        const int mid = (lo + hi + 1) / 2;
        if (decide(win[mid], nullptr)) lo = mid; else hi = mid - 1;
    }
    ConvexResult res;
    res.squared_distance = win[lo];
    decide(win[lo], &res.selected);
    res.stats = stats_;
    return res;
}

} // namespace cvx
