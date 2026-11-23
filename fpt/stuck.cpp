// Collect statistics on stuck states of the direct-exchange strategy (deterministic: ccw if valid else cw).
#include "../convex/kdisp.hpp"
#include "../convex/randgen.hpp"
#include <cmath>
#include <iostream>
#include <map>
using namespace kd;
static int n, K; static i64 R2; static const Instance* IP;
static bool far(int a, int b) { return IP->d(a, b) >= R2; }
static std::vector<std::vector<int>> fs;
static void rec(std::vector<int>& S, int next) {
    if ((int)S.size() == K) { fs.push_back(S); return; }
    for (int v = next; v < n; ++v) { bool ok = true; for (int u : S) ok = ok && far(u, v); if (ok) { S.push_back(v); rec(S, v + 1); S.pop_back(); } }
}
static int offq(int q, int v) { return ((v - q) % n + n) % n; }
std::map<std::string, long> st; long printed = 0; const std::vector<int>* curT; double minStuck = 1e9, maxMinAngle = 0, minMid = 1e9, minXUY = 1e9, maxA = 0;
// DFS over valid exchanges; record stuck states (both invalid) with their (depth, crossing) type
static bool dfs(std::vector<int>& C, std::vector<int>& rem, int root, int lo, int hi, bool record) {
    if ((int)C.size() == K) return true;
    int xs[2] = {-1, -1}, ts[2]; bool valid[2] = {false, false}; int blk[2] = {-1, -1};
    for (int side = 0; side < 2; ++side) {
        const bool ccw = side == 0; int x = -1, t = -1;
        if (ccw) { for (t = lo + 1; t < hi; ++t) { int v = (root + t) % n; bool ok = true; for (int c : C) ok = ok && far(c, v); if (ok) { x = v; break; } } }
        else { for (t = hi - 1; t > lo; --t) { int v = (root + t) % n; bool ok = true; for (int c : C) ok = ok && far(c, v); if (ok) { x = v; break; } } }
        xs[side] = x; ts[side] = t;
        if (x < 0) continue;
        std::vector<int> r2 = rem; if (ccw) r2.erase(r2.begin()); else r2.pop_back();
        valid[side] = true;
        for (size_t q = 0; q < r2.size(); ++q) if (!far(r2[q], x)) { valid[side] = false; blk[side] = ccw ? (int)q + 1 : (int)q; break; }
    }
    if (!valid[0] && !valid[1]) {
        if (record) {
            const int m = rem.size();
            // angles (cosines; larger cosine = smaller angle) in the ORIGINAL T
            auto cosang = [&](int prev, int p, int next) {
                const Pt &P = IP->p[p], &A = IP->p[prev], &B = IP->p[next];
                const long double ax = A.x - P.x, ay = A.y - P.y, bx = B.x - P.x, by = B.y - P.y;
                return (ax * bx + ay * by) / sqrtl((ax * ax + ay * ay) * (bx * bx + by * by)); };
            const std::vector<int>& T = *curT; const int k = T.size();
            auto idxT = [&](int v) { return (int)(std::find(T.begin(), T.end(), v) - T.begin()); };
            auto thetaT = [&](int v) { int i = idxT(v); return cosang(T[(i + k - 1) % k], v, T[(i + 1) % k]); };
            const long double cu = thetaT(root);
            bool smallerRem = false, smallerBetween = false;
            const int j = blk[0], i2 = blk[1] < 0 ? -1 : blk[1]; // rem indices (0-based): ccw blocked by rem[j], cw by rem[i2]
            for (int l = 0; l < m; ++l) if (thetaT(rem[l]) > cu) { smallerRem = true; if (l >= j && l <= i2) smallerBetween = true; }
            ++st[std::string("stuck: some remaining T-vertex has smaller T-angle than root = ") + (smallerRem ? "yes" : "NO")];
            ++st[std::string("stuck: a T-vertex between the blockers has smaller angle = ") + (smallerBetween ? "yes" : "NO")];
            ++st[std::string("stuck: blockers equal (i = j) = ") + (j == i2 ? "yes" : "no")];
            auto lab = [&](const char* nm, bool b) { ++st[std::string("stuck: ") + nm + (b ? " yes" : " NO")]; };
            lab("theta(t_j) < theta(u):", thetaT(rem[j]) > cu);
            lab("theta(t_i) < theta(u):", thetaT(rem[i2]) > cu);
            lab("theta(t_j) or theta(t_i) < theta(u):", thetaT(rem[j]) > cu || thetaT(rem[i2]) > cu);
            lab("theta(t_1) < theta(u):", thetaT(rem[0]) > cu);
            lab("theta(t_m) < theta(u):", thetaT(rem[m - 1]) > cu);
            lab("some t_l with l in [0,j] or [i,m-1] smaller:", [&]{ for (int l = 0; l < m; ++l) if ((l <= j || l >= i2) && thetaT(rem[l]) > cu) return true; return false; }());
            lab("depth>0 and root T-neighbours both replaced:", C.size() > 1);
            if (C.size() == 1 && m == 4 && printed < 12) { ++printed;
                auto deg = [&](long double c) { return (double)(acosl(c) * 180 / 3.14159265358979L); };
                std::cout << "DEPTH0 m=4 j=" << j + 1 << " i=" << i2 + 1 << " angles u " << deg(cu);
                for (int l = 0; l < m; ++l) std::cout << " t" << l + 1 << " " << deg(thetaT(rem[l]));
                std::cout << "\n";
            }
            { double dg = (double)(acosl(cu) * 180 / 3.14159265358979L); if (dg < minStuck) minStuck = dg;
              long double mx = -2, mn = 2; for (int v : T) { long double c = thetaT(v); mx = std::max(mx, c); mn = std::min(mn, c); }
              double minT = (double)(acosl(mx) * 180 / 3.14159265358979L); if (minT > maxMinAngle) maxMinAngle = minT; }
            if (C.size() == 1) {
                auto angAt = [&](int a, int p, int b) { return (double)(acosl(cosang(a, p, b)) * 180 / 3.14159265358979L); };
                minMid = std::min(minMid, angAt(rem[j], root, rem[i2]));
                minXUY = std::min(minXUY, angAt(xs[0], root, xs[1]));
                maxA = std::max({maxA, angAt(rem[0], root, rem[j]), angAt(rem[i2], root, rem[m-1])});
            }
            if (C.size() == 1) { bool ok1 = false; for (int l = 0; l < m; ++l) ok1 = ok1 || thetaT(rem[l]) > cu;
                lab("depth0: t1 or tm smaller:", thetaT(rem[0]) > cu || thetaT(rem[m-1]) > cu);
                lab("depth0: theta(t1)+theta(tm) < 2 theta(u):", acosl(thetaT(rem[0])) + acosl(thetaT(rem[m-1])) < 2 * acosl(cu)); }
            // blk[0] = index j (in rem, 0-based t_1 = 0) blocking ccw, blk[1] = index i blocking cw
            std::string type = "depth " + std::to_string(C.size() - 1) + (blk[1] < blk[0] ? " crossing" : " non-crossing") + " m=" + std::to_string(m);
            ++st[type];
        }
        return false;
    }
    for (int side = 0; side < 2; ++side) {
        if (!valid[side]) continue;
        const bool ccw = side == 0;
        std::vector<int> r2 = rem; if (ccw) r2.erase(r2.begin()); else r2.pop_back();
        C.push_back(xs[side]);
        const bool ok = dfs(C, r2, root, ccw ? ts[side] : lo, ccw ? hi : ts[side], record);
        C.pop_back();
        if (ok) return true;
    }
    return false;
}
int main(int argc, char** argv) {
    std::mt19937_64 g(atoll(argv[2]));
    for (long it = 0; it < atol(argv[1]); ++it) {
        Instance I(random_convex(g, it % 2 ? F_TWOCHAIN : (int)(g() % NUM_FAMILIES), 6 + g() % 7)); if (I.n < 6) continue;
        IP = &I; n = I.n;
        for (i64 r2 : I.candidates()) { R2 = r2;
            for (K = 3; K <= std::min(n, 7); ++K) {
                fs.clear(); std::vector<int> S; rec(S, 0); if (fs.empty()) break;
                for (auto& T : fs) for (int u : T) { curT = &T;
                    std::vector<int> C{u}, rem; for (int v : T) if (v != u) rem.push_back(v);
                    std::sort(rem.begin(), rem.end(), [&](int a, int b) { return offq(u, a) < offq(u, b); });
                    if (!dfs(C, rem, u, 0, n, true)) ++st["root fails"];
                    else ++st["root ok"];
                }
            }
        }
    }
    for (auto& [k, v] : st) std::cout << k << ": " << v << "\n";
    std::cout << "smallest root angle at a stuck state: " << minStuck << " deg; largest min-angle of T among stuck: " << maxMinAngle << " deg\n";
    std::cout << "depth0: min angle t_j u t_i " << minMid << ", min angle x u y " << minXUY << ", max partial angle " << maxA << "\n";
}
