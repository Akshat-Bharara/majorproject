// Instance generators: integer points in strictly convex position (ccw).
#pragma once
#include "arcdisp.hpp"
#include <cmath>
#include <random>

namespace arc {

using P64 = Pt<i64>;
enum Family { COCIRC, ROUNDCIRC, JITTER, DISK, SQUARE, GAUSS, ANNULUS, ELLIPSE, TWOCHAIN, NUM_FAMILIES };
inline const char* family_name(int f) {
    static const char* nm[] = {"cocircular", "roundcircle", "jitter", "disk", "square", "gauss", "annulus", "ellipse", "twochain"};
    return nm[f];
}
inline double urand(std::mt19937_64& g) { return (g() >> 11) * (1.0 / 9007199254740992.0); }

inline std::vector<P64> strict_hull(std::vector<P64> p) {
    std::sort(p.begin(), p.end(), [](const P64& a, const P64& b) { return a.x != b.x ? a.x < b.x : a.y < b.y; });
    p.erase(std::unique(p.begin(), p.end(), [](const P64& a, const P64& b) { return a.x == b.x && a.y == b.y; }), p.end());
    const int n = (int)p.size();
    if (n < 3) return p;
    auto cr = [](const P64& a, const P64& b, const P64& c) { return (i128)(b.x - a.x) * (c.y - a.y) - (i128)(b.y - a.y) * (c.x - a.x); };
    std::vector<P64> h(2 * n);
    int m = 0;
    for (int i = 0; i < n; ++i) { while (m >= 2 && cr(h[m - 2], h[m - 1], p[i]) <= 0) --m; h[m++] = p[i]; }
    for (int i = n - 2, t = m + 1; i >= 0; --i) { while (m >= t && cr(h[m - 2], h[m - 1], p[i]) <= 0) --m; h[m++] = p[i]; }
    h.resize(m - 1);
    return h;
}
inline std::vector<P64> subsample(std::mt19937_64& g, std::vector<P64> h, int m) {
    while ((int)h.size() > m) h.erase(h.begin() + g() % h.size());
    return h;
}
// all lattice points on x^2 + y^2 = N^2, N = 5*13*17*29*37 = 1185665 (972 points)
inline const std::vector<P64>& lattice_circle() {
    static std::vector<P64> pts;
    if (pts.empty()) {
        const i64 N = 1185665;
        std::vector<P64> q;
        for (i64 x = 0; x <= N; ++x) {
            const i64 y2 = N * N - x * x;
            i64 y = (i64)std::sqrt((long double)y2);
            while (y * y > y2) --y;
            while ((y + 1) * (y + 1) <= y2) ++y;
            if (y * y == y2) for (int sx = -1; sx <= 1; sx += 2) for (int sy = -1; sy <= 1; sy += 2) q.push_back({sx * x, sy * y});
        }
        pts = strict_hull(q);
    }
    return pts;
}
// n points at jittered equally spaced angles on a circle of radius R, rounded
inline std::vector<P64> round_circle(std::mt19937_64& g, int n, long double R) {
    const long double PI = 3.14159265358979323846264338327950288L;
    std::vector<P64> p(n);
    for (int i = 0; i < n; ++i) {
        const long double t = 2 * PI * (i + 0.5L * (long double)urand(g)) / n;
        p[i] = {(i64)llroundl(R * cosl(t)), (i64)llroundl(R * sinl(t))};
    }
    return p;
}
inline std::vector<P64> generate(std::mt19937_64& g, int fam, int n) {
    const double PI = 3.14159265358979323846;
    std::vector<P64> p;
    switch (fam) {
    case COCIRC: return subsample(g, lattice_circle(), n);
    case ROUNDCIRC: return strict_hull(round_circle(g, n, 5.0e8L));
    case JITTER: { // regular polygon, radius 10^4..10^6, integer jitter of random magnitude
        const double R = std::pow(10.0, 4 + 2 * urand(g)), J = R * std::pow(10.0, -4 + 3.3 * urand(g));
        for (int i = 0; i < n; ++i)
            p.push_back({(i64)llround(R * cos(2 * PI * i / n) + J * (2 * urand(g) - 1)), (i64)llround(R * sin(2 * PI * i / n) + J * (2 * urand(g) - 1))});
        return strict_hull(p);
    }
    case DISK:
        while ((int)p.size() < 30 * n) {
            const i64 x = (i64)(g() % 2000001) - 1000000, y = (i64)(g() % 2000001) - 1000000;
            if (x * x + y * y <= 1000000LL * 1000000LL) p.push_back({x, y});
        }
        break;
    case SQUARE:
        for (int i = 0; i < 60 * n; ++i) p.push_back({(i64)(g() % 2000001) - 1000000, (i64)(g() % 2000001) - 1000000});
        break;
    case GAUSS: {
        std::normal_distribution<double> N(0, 2e5);
        for (int i = 0; i < 3000 * n; ++i) p.push_back({(i64)llround(N(g)), (i64)llround(N(g))});
        break;
    }
    case ANNULUS:
        for (int i = 0; i < 4 * n; ++i) {
            const double r = 1e6 * std::sqrt(0.9025 + 0.0975 * urand(g)), t = 2 * PI * urand(g);
            p.push_back({(i64)llround(r * cos(t)), (i64)llround(r * sin(t))});
        }
        break;
    case ELLIPSE: // points on a thin 1:20 ellipse
        for (int i = 0; i < n; ++i) {
            const double t = 2 * PI * urand(g);
            p.push_back({(i64)llround(1e6 * cos(t)), (i64)llround(5e4 * sin(t))});
        }
        break;
    default: { // two-chain lens
        const double h = 2e4 * (0.5 + urand(g)), h2 = h * (0.4 + urand(g));
        for (int i = 0; i < 3 * n; ++i) {
            const double x = (2 * urand(g) - 1) * 1e6, s = 1 - (x / 1e6) * (x / 1e6);
            p.push_back({(i64)llround(x), (i64)llround((g() & 1) ? h * s : -h2 * s)});
        }
        break;
    }
    }
    return subsample(g, strict_hull(p), n);
}
// the 8-point set P* of the earlier paper (no tight pair for k = 4 at r^2 = 10240)
inline std::vector<P64> pstar() {
    return {{16, 48}, {-83, 26}, {-114, 18}, {-30, -39}, {-16, -48}, {83, -26}, {114, -18}, {30, 39}};
}
} // namespace arc
