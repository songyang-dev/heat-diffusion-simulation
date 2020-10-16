#ifndef __iglViewer_h__
#define __iglViewer_h__

#include "trimesh.h"
#include <Eigen/Dense>

void view(const trimesh::trimesh_t& mesh, const Eigen::MatrixXd& coloring);
void viewSmoothing(const std::vector<Eigen::MatrixXd>& vertices, const Eigen::MatrixXi& faces);
#endif