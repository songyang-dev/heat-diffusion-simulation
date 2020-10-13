/**
 * Computes the laplacian of a mesh
 * */

#include "laplacian.h"
#include "utils.h"

#include <Eigen/Sparse>

EdgeToCotan computeCotangents(const trimesh::trimesh_t& mesh);
double edgeLength(const trimesh::trimesh_t& mesh, const trimesh::edge_t& edge);

// Computes the cotangent Laplacian of the mesh
Eigen::SparseMatrix<double> computeLaplacian(const trimesh::trimesh_t& mesh)
{
    // to be filled with triplets
    Eigen::SparseMatrix<double> laplacian(mesh.Vertices.rows(), mesh.Vertices.rows());

    std::vector<Triplet> matrix_elements;
    matrix_elements.reserve((6+1) * mesh.Vertices.rows()); // degree 6 vertices

    // 1. Compute off-diagonal values.
    // 2. Compute diagonal values.

    // 1 Off-diagonal values
    // 1a. Compute the two associated cotangents of every edge
    EdgeToCotan cotangents = computeCotangents(mesh);

    return laplacian;
}

// Computes the cotangents for every edge of the mesh
EdgeToCotan computeCotangents(const trimesh::trimesh_t& mesh)
{
    EdgeToCotan cotans;

    // loop on edges
    for (size_t i = 0; i < mesh.Edges.size(); i++)
    {
        auto edge = mesh.Edges[i];

        auto initialEdgeLength = edgeLength(mesh, edge);
        
    }
    
    return cotans;
}

// Computes the length of the edge
double edgeLength(const trimesh::trimesh_t& mesh, const trimesh::edge_t& edge)
{
    
}