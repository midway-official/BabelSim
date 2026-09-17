#pragma once

#include <Eigen/IterativeLinearSolvers>

#include <stdexcept>
#include <vector>

namespace babelsim::detail {

// Eigen's sparse preconditioner solve expression is correct but creates an
// intermediate vector for each permutation/triangular assignment.  These
// backend-only wrappers keep the same factor and operation order while
// reusing the caller's output and one permutation scratch vector.  Physics and
// the public solver API do not see these types.
class InPlaceIncompleteCholesky final : public Eigen::IncompleteCholesky<double> {
    using Base = Eigen::IncompleteCholesky<double>;

public:
    template <typename MatrixType>
    void compute(const MatrixType& matrix) {
        Base::compute(matrix);
        refreshPermutation();
    }

    template <typename MatrixType>
    void factorize(const MatrixType& matrix) {
        Base::factorize(matrix);
        refreshPermutation();
    }

    void solveInPlace(const Eigen::VectorXd& input, Eigen::VectorXd& output) const {
        if (input.size() != this->rows()) {
            throw std::invalid_argument("incomplete Cholesky vector size is invalid");
        }
        if (&input == &output) {
            throw std::invalid_argument("incomplete Cholesky input aliases output");
        }
        if (output.size() != input.size()) output.resize(input.size());
        const Eigen::Index size = input.size();
        if (permutation_forward.empty()) {
            output = input;
        } else {
            for (Eigen::Index index = 0; index < size; ++index) {
                output[permutation_forward[static_cast<std::size_t>(index)]] = input[index];
            }
        }
        output.array() *= this->scalingS().array();
        this->matrixL().template triangularView<Eigen::Lower>().solveInPlace(output);
        this->matrixL().adjoint().template triangularView<Eigen::Upper>().solveInPlace(output);
        output.array() *= this->scalingS().array();
        if (!permutation_forward.empty()) {
            scratch.resize(input.size());
            for (Eigen::Index index = 0; index < size; ++index) {
                scratch[permutation_inverse[static_cast<std::size_t>(index)]] = output[index];
            }
            output.swap(scratch);
        }
    }

private:
    void refreshPermutation() {
        const auto& permutation = this->permutationP();
        permutation_forward.clear();
        permutation_inverse.clear();
        if (permutation.size() == 0) return;
        permutation_forward.resize(static_cast<std::size_t>(permutation.size()));
        permutation_inverse.resize(static_cast<std::size_t>(permutation.size()));
        for (Eigen::Index index = 0; index < permutation.size(); ++index) {
            permutation_forward[static_cast<std::size_t>(index)] =
                permutation.indices()[index];
            permutation_inverse[static_cast<std::size_t>(
                permutation.indices()[index])] = index;
        }
    }

    std::vector<Eigen::Index> permutation_forward;
    std::vector<Eigen::Index> permutation_inverse;
    mutable Eigen::VectorXd scratch;
};

class InPlaceIncompleteLut final : public Eigen::IncompleteLUT<double> {
    using Base = Eigen::IncompleteLUT<double>;

public:
    template <typename MatrixType>
    void compute(const MatrixType& matrix) {
        Base::compute(matrix);
        refreshPermutation();
    }

    template <typename MatrixType>
    void factorize(const MatrixType& matrix) {
        Base::factorize(matrix);
        refreshPermutation();
    }

    void solveInPlace(const Eigen::VectorXd& input, Eigen::VectorXd& output) const {
        if (input.size() != this->rows()) {
            throw std::invalid_argument("ILUT vector size is invalid");
        }
        if (&input == &output) {
            throw std::invalid_argument("ILUT input aliases output");
        }
        if (output.size() != input.size()) output.resize(input.size());
        const Eigen::Index size = input.size();
        if (permutation_inverse.empty()) {
            output = input;
        } else {
            for (Eigen::Index index = 0; index < size; ++index) {
                output[permutation_inverse[static_cast<std::size_t>(index)]] = input[index];
            }
        }
        this->m_lu.template triangularView<Eigen::UnitLower>().solveInPlace(output);
        this->m_lu.template triangularView<Eigen::Upper>().solveInPlace(output);
        if (!permutation_forward.empty()) {
            scratch.resize(input.size());
            for (Eigen::Index index = 0; index < size; ++index) {
                scratch[permutation_forward[static_cast<std::size_t>(index)]] = output[index];
            }
            output.swap(scratch);
        }
    }

private:
    void refreshPermutation() {
        if (this->m_P.size() != this->m_Pinv.size()) {
            throw std::logic_error("ILUT permutations have different sizes");
        }
        if (this->m_P.size() == 0) return;
        permutation_forward.resize(static_cast<std::size_t>(this->m_P.size()));
        permutation_inverse.resize(static_cast<std::size_t>(this->m_Pinv.size()));
        for (Eigen::Index index = 0; index < this->m_P.size(); ++index) {
            permutation_forward[static_cast<std::size_t>(index)] = this->m_P.indices()[index];
            permutation_inverse[static_cast<std::size_t>(index)] = this->m_Pinv.indices()[index];
        }
    }

    std::vector<Eigen::Index> permutation_forward;
    std::vector<Eigen::Index> permutation_inverse;
    mutable Eigen::VectorXd scratch;
};

}  // namespace babelsim::detail
