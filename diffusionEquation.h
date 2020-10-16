#ifndef __diffusionEquation_h__
#define __diffusionEquation_h__

#include <vector>
#include <Eigen/Dense>

#include "trimesh.h"
#include "utils.h"

// Computes the diffusion equation "smoothing" of a mesh
std::vector<Eigen::MatrixXd> simulateDiffusion(size_t iterations,
    const Eigen::SparseMatrix<double>& cotangents,
    const DiagonalXd& mass,
    const trimesh::trimesh_t& initialMesh);

#endif