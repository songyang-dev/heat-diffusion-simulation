#ifndef __heatEquation_h__
#define __heatEquation_h__

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "utils.h"

Eigen::MatrixXd simulateHeat(size_t iterations,
    const Eigen::SparseMatrix<double>& cotangents,
    const DiagonalXd& mass,
    const Eigen::VectorXd& initialTemperature);

#endif