#ifndef __laplacian_h__
#define __laplacian_h__

#include <Eigen/Sparse>
#include "trimesh.h"

typedef Eigen::Triplet<double> Triplet;

Eigen::SparseMatrix<double> computeLaplacian(const trimesh::trimesh_t& mesh);

#endif