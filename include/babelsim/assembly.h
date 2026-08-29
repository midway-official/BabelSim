#pragma once

#include "babelsim/equation.h"

#include <Eigen/Sparse>

#include <array>
#include <vector>

namespace babelsim {

struct LinearSystem {
    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd b;
};

// Precomputes the mesh-dependent compressed sparse structure and the direct
// value-array positions of every LDU coefficient. On a decomposed mesh the
// matrix contains owned rows only; cross-rank face entries are applied by the
// distributed solver's halo matvec. update() therefore performs no topology
// traversal, triplet allocation, sorting, or sparse insertion.
class SparseAssembly {
public:
    explicit SparseAssembly(const Mesh& mesh);

    void update(const ScalarEquation& equation);
    void update(const VectorEquation& equation);
    const Eigen::SparseMatrix<double>& matrix() const { return matrix_; }

private:
    void update(
        const Mesh* equation_mesh,
        const std::vector<double>& diagonal,
        const std::vector<double>& upper,
        const std::vector<double>& lower);

    const Mesh* mesh_;
    Eigen::SparseMatrix<double> matrix_;
    std::vector<Eigen::Index> diagonal_positions_;
    std::vector<Eigen::Index> upper_positions_;
    std::vector<Eigen::Index> lower_positions_;
};

Eigen::SparseMatrix<double> assembleMatrix(const ScalarEquation& equation);
Eigen::SparseMatrix<double> assembleMatrix(const VectorEquation& equation);
LinearSystem assemble(const ScalarEquation& equation);
void assembleSource(const ScalarEquation& equation, Eigen::VectorXd& result);
void assembleSource(
    const VectorEquation& equation,
    std::array<Eigen::VectorXd, 3>& result);
std::array<Eigen::VectorXd, 3> assembleSource(const VectorEquation& equation);

}  // namespace babelsim
