#ifndef __examine_h__
#define __examine_h__

#include <Eigen/Sparse>

void examineLaplacian(const Eigen::SparseMatrix<double>& cotangents);
void examineEigenDecomposition(const Eigen::MatrixXd& eigenvectors,
    const Eigen::VectorXd& eigenvalues);
void examineMass(const Eigen::SparseMatrix<double>& mass);
#endif