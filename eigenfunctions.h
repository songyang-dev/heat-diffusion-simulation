#ifndef __eigenfunctions_h__
#define __eigenfunctions_h__

#include <Eigen/Sparse>

std::pair<Eigen::MatrixXd, Eigen::VectorXd> eigenDecomposition(
    const Eigen::SparseMatrix<double>& matrix, int numEigensToSee);

#endif