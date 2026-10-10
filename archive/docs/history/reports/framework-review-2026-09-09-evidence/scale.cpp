#include "babelsim/linear_solver.h"
#include <iostream>
using namespace babelsim;
int main() {
 for (auto type : {LinearSolverType::ConjugateGradient, LinearSolverType::BiCGSTAB}) {
  Eigen::SparseMatrix<double> A(1,1); A.insert(0,0)=1e-20;
  Eigen::VectorXd b(1),x; b[0]=1e-20;
  LinearSolverConfig c; c.solver=type; c.preconditioner=PreconditionerType::None;
  c.absolute_tolerance=1e-30; c.relative_tolerance=1e-12;
  auto r=solve(A,b,x,c);
  std::cout << "solver=" << static_cast<int>(type) << " converged=" << r.converged() << " iterations=" << r.iterations << " x=" << x[0] << " expected=1 relativeResidual=" << r.relative_residual << '\n';
 }
}
