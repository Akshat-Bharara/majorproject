// Random point sets in strictly convex position (used by selftest and hunt).
#pragma once

#include "kdisp.hpp"
#include <random>

namespace kd {

enum Family { F_SQUARE, F_DISK, F_CIRCLE, F_ELLIPSE, F_TWOCHAIN, F_SMALLGRID, NUM_FAMILIES };

inline const char* family_name(int f) {
    static const char* names[] = {"square", "disk", "circle", "ellipse", "twochain", "smallgrid"};
    return names[f];
}

inline double urand(std::mt19937_64& rng) { return (rng() >> 11) * (1.0 / 9007199254740992.0); }

// keep a uniformly random subset of size m (cyclic order preserved)
inline std::vector<Pt> subsample(std::mt19937_64& rng, std::vector<Pt> h, int m) {
    while ((int)h.size() > m) h.erase(h.begin() + rng() % h.size());
    return h;
}

// Two nearly flat convex chains: a thin lens of width 2W and height ~2h.
inline std::vector<Pt> two_chain(std::mt19937_64& rng, int m, double W, double h) {
    std::vector<Pt> p;
    const double h2 = h * (0.3 + 1.4 * urand(rng));
    for (int i = 0; i < 3 * m; ++i) {
        const double x = (2 * urand(rng) - 1) * W;
        const double s = 1 - (x / W) * (x / W);
        const bool top = rng() & 1;
        p.push_back({(i64)llround(x), (i64)llround(top ? h * s : -h2 * s)});
    }
    return subsample(rng, strict_hull(p), m);
}

// A random convex-position set with roughly `target` points (never more).
inline std::vector<Pt> random_convex(std::mt19937_64& rng, int fam, int target) {
    const double PI = 3.14159265358979323846;
    std::vector<Pt> p;
    switch (fam) {
    case F_SQUARE:
        for (int i = 0; i < 60 * target; ++i)
            p.push_back({(i64)(rng() % 2000001) - 1000000, (i64)(rng() % 2000001) - 1000000});
        break;
    case F_DISK:
        while ((int)p.size() < 30 * target) {
            const i64 x = (i64)(rng() % 2000001) - 1000000, y = (i64)(rng() % 2000001) - 1000000;
            if (x * x + y * y <= 1000000LL * 1000000LL) p.push_back({x, y});
        }
        break;
    case F_CIRCLE:
        for (int i = 0; i < target; ++i) {
            const double t = 2 * PI * urand(rng);
            p.push_back({(i64)llround(1e6 * cos(t)), (i64)llround(1e6 * sin(t))});
        }
        break;
    case F_ELLIPSE: {
        const double b = 1e6 * (0.01 + 0.2 * urand(rng));
        for (int i = 0; i < target; ++i) {
            const double t = 2 * PI * urand(rng);
            p.push_back({(i64)llround(1e6 * cos(t)), (i64)llround(b * sin(t))});
        }
        break;
    }
    case F_TWOCHAIN:
        return two_chain(rng, target, 1e5, 300 + 4000 * urand(rng));
    default: // small integer grid: many equal distances and cocircular quadruples
        for (int i = 0; i < 40 * target; ++i) p.push_back({(i64)(rng() % 41) - 20, (i64)(rng() % 41) - 20});
        break;
    }
    return subsample(rng, strict_hull(p), target);
}

} // namespace kd
