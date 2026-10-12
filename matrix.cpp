#include "matrix.hpp"
#include <algorithm>

namespace fmm {

Matrix::Matrix(std::size_t rows, std::size_t cols, std::uint64_t value)
    : rows_(rows), cols_(cols), data_(rows * cols, value) {}

Matrix multiply_blocked(const Matrix& A, const Matrix& B, std::size_t block) {
    if (A.cols() != B.rows()) {
        throw std::invalid_argument("matrix dimension mismatch");
    }

    Matrix C(A.rows(), B.cols(), 0);

    const std::size_t m = A.rows();
    const std::size_t k = A.cols();
    const std::size_t n = B.cols();

    for (std::size_t ii = 0; ii < m; ii += block) {
        const std::size_t i_end = std::min(ii + block, m);
        for (std::size_t kk = 0; kk < k; kk += block) {
            const std::size_t k_end = std::min(kk + block, k);
            for (std::size_t jj = 0; jj < n; jj += block) {
                const std::size_t j_end = std::min(jj + block, n);

                for (std::size_t i = ii; i < i_end; ++i) {
                    auto* c = C.row(i);
                    const auto* a = A.row(i);
                    for (std::size_t t = kk; t < k_end; ++t) {
                        const std::uint64_t av = a[t];
                        if (av == 0) continue;
                        const auto* b = B.row(t);
                        for (std::size_t j = jj; j < j_end; ++j) {
                            c[j] += av * b[j];
                        }
                    }
                }
            }
        }
    }
    return C;
}

} // namespace fmm
