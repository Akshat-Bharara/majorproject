#pragma once

// k-dispersion for planar points in convex position.
//
// Two exact solvers share the same interface:
//
//  * DelaunayCappedDP (this project): triangulation DP over locally Delaunay
//    triangulations. For a chord (i,j) it stores only the k-2 smallest apex
//    angles needed to reach 1..k-2 points (as apex indices, compared with
//    exact in-circle tests), so one decision costs O(n^3 log k) time and
//    O(k n^2) space.
//
//  * TWNaiveDP: the dynamic program of Tkachenko and Wang (STACS 2025,
//    Lemma 13) with subproblems f(i,j,k) over canonical triples, evaluated
//    directly in O(n^4) time per decision (the O(n^{7/2}) version needs
//    cuttings and is not implemented).
//
// Both optimise with a binary search over the distinct pairwise squared
// distances, restricted to the window [L, U] of Lemma "bounds" in the report:
//   L = min pairwise distance of the greedy (farthest-first) k-set,
//   U = min(2 L, perimeter / k)          (squared values compared exactly).

#include "geometry.hpp"
#include <cstdint>
#include <string>
#include <vector>

namespace cvx {

struct ConvexStats {
    std::size_t candidates_total = 0;   // distinct pairwise distances
    std::size_t candidates_window = 0;  // of those inside [L, U]
    int decisions = 0;                  // exact decision-procedure calls
    int greedy_decisions = 0;           // seeded-greedy calls (Certified)
    bool certified_by_greedy = false;   // Certified: greedy value was optimal
    std::uint64_t inner_steps = 0;      // apex evaluations (work measure)
    i64 lower_bound2 = 0, upper_bound2 = 0;
};

struct ConvexResult {
    i64 squared_distance = 0;
    std::vector<int> selected;          // indices into the ccw-ordered points
    ConvexStats stats;
};

// Certified: binary search with the O(n^2 k) seeded greedy (a lower bound that
// is always valid), then one exact CappedDP decision just above it. If that
// decision fails the greedy value is optimal; otherwise the exact binary search
// continues above it. Always exact.
enum class Algo { CappedDP, TWNaive, Certified };

class ConvexDispersion {
public:
    // pts must be in strictly convex position, counterclockwise order.
    ConvexDispersion(std::vector<Pt> pts, int k);

    ConvexResult solve(Algo algo, bool use_bounds = true);

    // Decision: is there a k-subset with all pairwise squared distances >= r2?
    // On success fills witness (ccw indices).
    // tight_only: accept only solutions containing a tight pair q, a(q)
    // (base chord (a(q), q)); used to test conjecture C1.
    bool decide_capped(i64 r2, std::vector<int>* witness, bool tight_only = false);
    bool decide_tw(i64 r2, std::vector<int>* witness);

    // Seeded two-sided greedy (sound, not complete): seeds {q, a(q), b(q)} and
    // {q, a(q)} for every q, frontiers extended alternately (both orders).
    // True only if it found k points pairwise >= r (returned in witness).
    bool decide_greedy(i64 r2, std::vector<int>* witness);

    // Greedy farthest-first k-set (deterministic start: a diametral point).
    std::vector<int> greedy() const;

private:
    std::vector<Pt> p_;
    int n_, k_;
    std::vector<i64> d2_;     // n*n squared distances
    ConvexStats stats_;

    i64 D(int a, int b) const { return d2_[(std::size_t)a * n_ + b]; }
    int nxt(int i, int s) const { int v = i + s; return v >= n_ ? v - n_ : v; }
    int gap(int i, int j) const { int g = j - i; return g < 0 ? g + n_ : g; }
};

} // namespace cvx
