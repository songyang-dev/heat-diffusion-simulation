#if !defined(__utils_h__)
#define __utils_h__

#include "trimesh_types.h"
#include <vector>
#include <Eigen/Dense>

// maps

// From edge index to cotangents
typedef std::vector<std::pair<double, double>> EdgeToCotan;

// Dynamic diagonal matrix
typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic, Eigen::Dynamic> DiagonalXd;

// Covert cosine of an angle to cotangent
double cosToCot(double cosine);

#endif // __utils_h__
