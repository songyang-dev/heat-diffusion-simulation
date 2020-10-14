/**
 * Computes the laplacian of a mesh
 * */

#include "laplacian.h"
#include "utils.h"

#include <Eigen/Sparse>

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
    cotans.reserve(mesh.Edges.size());

    // loop on edges
    for (size_t i = 0; i < mesh.Edges.size(); i++)
    {
        auto edge = mesh.Edges[i];
        auto halfedge = mesh.directed_edge2he_index(edge.start(), edge.end());

        auto initialEdgeLength = edgeLength(mesh, edge);
        
        // alpha
        auto nextEdge = edgeLength(mesh, 
            mesh.he_index2directed_edge(
                mesh.halfedge(halfedge).next_he
            )
        );
        auto nextNextEdge = edgeLength(mesh,
            mesh.he_index2directed_edge(
                mesh.halfedge(mesh.halfedge(halfedge).next_he).next_he
            )  
        );
        double cosineAlpha =  (-std::pow(initialEdgeLength, 2) + std::pow(nextEdge, 2) + std::pow(nextNextEdge, 2))
            / (2 * nextEdge * nextNextEdge);
        auto cotAlpha = cosToCot(cosineAlpha);

        // beta
        auto oppositeNextEdge = edgeLength(mesh,
            mesh.he_index2directed_edge(
                mesh.halfedge(mesh.halfedge(halfedge).opposite_he).next_he
            )
        );
        auto oppositeNextNextEdge = edgeLength(mesh,
            mesh.he_index2directed_edge(
                mesh.halfedge(mesh.halfedge(mesh.halfedge(
                    halfedge).opposite_he).next_he).next_he
            )
        );
        double cosineBeta =  (-std::pow(initialEdgeLength, 2) + std::pow(oppositeNextEdge, 2) + std::pow(oppositeNextNextEdge, 2))
            / (2 * oppositeNextEdge * oppositeNextNextEdge);
        auto cotBeta = cosToCot(cosineBeta);

        cotans.push_back(std::make_pair(cotAlpha, cotBeta));
    }
    
    assert(cotans.size(), mesh.Edges.size());
    return cotans;
}

// Computes the length of the edge
double edgeLength(const trimesh::trimesh_t& mesh, const trimesh::edge_t& edge)
{
    Eigen::Vector3d vector = mesh.Vertices.row(edge.end()) - mesh.Vertices.row(edge.start());
    return vector.norm();
}

// Computes the length of the directed edge, unsigned result
double edgeLength(const trimesh::trimesh_t& mesh, const std::pair<trimesh::index_t, trimesh::index_t>& pair)
{
    Eigen::Vector3d vector = mesh.Vertices.row(pair.second) - mesh.Vertices.row(pair.first);
    return vector.norm();
}