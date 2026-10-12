#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>
#include <stdexcept>

namespace fmm {

// Dense row-major matrix. The paper's triangle reduction needs exact
// multiplication over nonnegative integers: C = A * B, followed by C[i,j] > 0.
//
// IMPORTANT: this class is the arithmetic backend of the reduction. The paper
// does not specify a low-level implementation of the rectangular FMM itself;
// see README.md. This default backend is an exact blocked classical GEMM.
class Matrix {
public:
    Matrix() = default;
    Matrix(std::size_t rows, std::size_t cols, std::uint64_t value = 0);

    std::size_t rows() const noexcept { return rows_; }
    std::size_t cols() const noexcept { return cols_; }

    std::uint64_t* row(std::size_t i) { return data_.data() + i * cols_; }
    const std::uint64_t* row(std::size_t i) const { return data_.data() + i * cols_; }

    std::uint64_t& operator()(std::size_t i, std::size_t j) {
        return data_[i * cols_ + j];
    }
    std::uint64_t operator()(std::size_t i, std::size_t j) const {
        return data_[i * cols_ + j];
    }

    const std::vector<std::uint64_t>& data() const noexcept { return data_; }
    std::vector<std::uint64_t>& data() noexcept { return data_; }

private:
    std::size_t rows_ = 0;
    std::size_t cols_ = 0;
    std::vector<std::uint64_t> data_;
};

// Exact dense matrix product. Dimensions:
// A: m x k, B: k x n, C: m x n.
//
// This is deliberately isolated behind one function. To obtain the paper's
// asymptotic FMM bounds, replace this routine by an implementation of the
// corresponding rectangular matrix multiplication algorithm.
Matrix multiply_blocked(const Matrix& A, const Matrix& B, std::size_t block = 64);

} // namespace fmm
