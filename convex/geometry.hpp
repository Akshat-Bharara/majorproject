#pragma once

// Exact planar predicates for integer points.
//
// All coordinates must satisfy |x|, |y| <= kMaxCoord = 2^29. Then every
// coordinate difference fits in 31 bits, squared lengths in 62 bits, and the
// 3x3 in-circle determinant (sum of three products of a 62-bit and a 63-bit
// term) fits in a signed 128-bit integer, so the predicates below are exact.

#include <cstdint>
#include <vector>

namespace cvx {

using i64 = std::int64_t;
using i128 = __int128;

constexpr i64 kMaxCoord = i64(1) << 29;

struct Pt {
    i64 x = 0, y = 0;
};

inline i64 dist2(const Pt& a, const Pt& b) {
    const i64 dx = a.x - b.x, dy = a.y - b.y;
    return dx * dx + dy * dy;
}

// > 0 if a, b, c make a left (counterclockwise) turn, < 0 right, 0 collinear.
inline i128 orient(const Pt& a, const Pt& b, const Pt& c) {
    return (i128)(b.x - a.x) * (c.y - a.y) - (i128)(b.y - a.y) * (c.x - a.x);
}

// Raw in-circle determinant: for a, b, c counterclockwise it is > 0 iff d lies
// strictly inside their circumcircle, = 0 iff d is on it.
inline i128 incircle_det(const Pt& a, const Pt& b, const Pt& c, const Pt& d) {
    const i128 adx = a.x - d.x, ady = a.y - d.y;
    const i128 bdx = b.x - d.x, bdy = b.y - d.y;
    const i128 cdx = c.x - d.x, cdy = c.y - d.y;
    const i128 alift = adx * adx + ady * ady;
    const i128 blift = bdx * bdx + bdy * bdy;
    const i128 clift = cdx * cdx + cdy * cdy;
    return alift * (bdx * cdy - bdy * cdx)
         - blift * (adx * cdy - ady * cdx)
         + clift * (adx * bdy - ady * bdx);
}

// True iff d lies strictly inside the circle through a, b, c (any orientation;
// a, b, c must not be collinear).
inline bool strictly_inside_circle(const Pt& a, const Pt& b, const Pt& c, const Pt& d) {
    const i128 o = orient(a, b, c);
    const i128 det = incircle_det(a, b, c, d);
    return o > 0 ? det > 0 : det < 0;
}

// Indices of the strictly convex hull of pts in counterclockwise order
// (Andrew's monotone chain). Collinear boundary points are not hull vertices.
std::vector<int> convex_hull_ccw(const std::vector<Pt>& pts);

} // namespace cvx
