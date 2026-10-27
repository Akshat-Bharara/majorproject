#include "geometry.hpp"
#include <algorithm>
#include <numeric>

namespace cvx {

std::vector<int> convex_hull_ccw(const std::vector<Pt>& pts) {
    const int n = (int)pts.size();
    std::vector<int> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b) {
        return pts[a].x != pts[b].x ? pts[a].x < pts[b].x : pts[a].y < pts[b].y;
    });
    // Drop exact duplicates (keep the first occurrence).
    idx.erase(std::unique(idx.begin(), idx.end(), [&](int a, int b) {
        return pts[a].x == pts[b].x && pts[a].y == pts[b].y;
    }), idx.end());
    if (idx.size() < 3) return idx;

    std::vector<int> h(2 * idx.size());
    int m = 0;
    for (int v : idx) {  // lower hull
        while (m >= 2 && orient(pts[h[m - 2]], pts[h[m - 1]], pts[v]) <= 0) --m;
        h[m++] = v;
    }
    for (int t = (int)idx.size() - 2, lo = m + 1; t >= 0; --t) {  // upper hull
        const int v = idx[t];
        while (m >= lo && orient(pts[h[m - 2]], pts[h[m - 1]], pts[v]) <= 0) --m;
        h[m++] = v;
    }
    h.resize(m - 1);
    return h;
}

} // namespace cvx
