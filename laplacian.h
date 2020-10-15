#ifndef __laplacian_h__
#define __laplacian_h__

#include <Eigen/Sparse>
#include "trimesh.h"
#include "utils.h"

typedef Eigen::Triplet<double> Triplet;

Eigen::SparseMatrix<double> computeCotangentLaplacian(const trimesh::trimesh_t& mesh);
DiagonalXd computeMassMatrix(const trimesh::trimesh_t& mesh);

EdgeToCotan computeCotangents(const trimesh::trimesh_t& mesh);
double edgeLength(const trimesh::trimesh_t& mesh, const trimesh::edge_t& edge);
double edgeLength(const trimesh::trimesh_t& mesh, const std::pair<trimesh::index_t, trimesh::index_t>& pair);
double faceArea(const trimesh::trimesh_t& mesh, const size_t& face);
#endif