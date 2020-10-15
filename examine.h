#ifndef __examine_h__
#define __examine_h__

#include <Eigen/Sparse>

#include "utils.h"

void examineLaplacian(const Eigen::SparseMatrix<double>& cotangents);
void examineEigenDecomposition(const Eigen::MatrixXd& eigenvectors,
    const Eigen::VectorXd& eigenvalues);
void examineMass(const DiagonalXd& mass);
void examineTemperatures(const Eigen::MatrixXd& temperatures);
#endif