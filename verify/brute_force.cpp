// Exhaustive k-dispersion oracle used to cross-check the FMM solver.
//
// Enumerates all C(n,k) subsets in lexicographic order and keeps the one
// whose minimum pairwise squared distance is largest (first one on ties).
// It shares no code with the FMM solver, so agreement between the two is
// independent evidence of correctness.
//
// Input/output format is identical to the main driver:
//   n d k
//   x_1 ... x_d   (n lines, integral coordinates)
// prints
//   optimal_squared_distance <value>
//   selected_indices i_1 ... i_k
// and the time spent after reading the input on stderr ("solve_ms <ms>").

#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <algorithm>

using u128 = unsigned __int128;

static std::string to_string_u128(u128 x) {
    if (x == 0) return "0";
    std::string s;
    while (x) {
        s.push_back(char('0' + int(x % 10)));
        x /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

int main() {
    int n, d, k;
    if (!(std::cin >> n >> d >> k) || n <= 0 || d <= 0 || k < 2 || k > n) {
        std::cerr << "invalid header\n";
        return 1;
    }
    std::vector<std::vector<std::int64_t>> p(n, std::vector<std::int64_t>(d));
    for (auto& row : p)
        for (auto& x : row)
            if (!(std::cin >> x)) {
                std::cerr << "invalid/incomplete point data\n";
                return 1;
            }

    const auto t0 = std::chrono::steady_clock::now();
    std::vector<std::vector<u128>> D(n, std::vector<u128>(n, 0));
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            u128 s = 0;
            for (int t = 0; t < d; ++t) {
                const __int128 diff = (__int128)p[i][t] - (__int128)p[j][t];
                s += (u128)(diff * diff);
            }
            D[i][j] = D[j][i] = s;
        }

    // idx = current combination; cur[i] = min distance among idx[0..i].
    std::vector<int> idx(k), best_set;
    std::vector<u128> cur(k);
    const u128 INF = ~(u128)0;
    u128 best = 0;
    bool have = false;
    for (int i = 0; i < k; ++i) idx[i] = i;

    auto prefix_min = [&](int from) {
        for (int i = from; i < k; ++i) {
            u128 m = i ? cur[i - 1] : INF;
            for (int j = 0; j < i; ++j) m = std::min(m, D[idx[i]][idx[j]]);
            cur[i] = m;
        }
    };
    prefix_min(0);

    while (true) {
        if (!have || cur[k - 1] > best) {
            best = cur[k - 1];
            best_set = idx;
            have = true;
        }
        // Advance to the next combination.
        int i = k - 1;
        while (i >= 0 && idx[i] == n - k + i) --i;
        if (i < 0) break;
        ++idx[i];
        for (int j = i + 1; j < k; ++j) idx[j] = idx[j - 1] + 1;
        prefix_min(i);
    }

    const double solve_ms =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    std::cerr << "solve_ms " << solve_ms << '\n';

    std::cout << "optimal_squared_distance " << to_string_u128(best) << '\n';
    std::cout << "selected_indices";
    for (int v : best_set) std::cout << ' ' << v;
    std::cout << '\n';
    return 0;
}
