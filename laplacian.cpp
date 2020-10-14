/**
 * Computes the laplacian of a mesh
 * */

#include "laplacian.h"
#include "utils.h"

#include <Eigen/Sparse>

// Computes the cotangent Laplacian of the mesh
Eigen::SparseMatrix<double> computeLaplacian(const trimesh::trimesh_t &mesh)
{
    // to be filled with triplets
    Eigen::SparseMatrix<double> laplacian(mesh.Vertices.rows(), mesh.Vertices.rows());

    std::vector<Triplet> matrix_non_diagonal_elements;
    matrix_non_diagonal_elements.reserve(6 * mesh.Vertices.rows()); // degree 6 vertices

    // 1. Compute off-diagonal values.
    // 2. Compute diagonal values.

    // 1 Off-diagonal values
    // 1a. Compute the two associated cotangents of every edge
    EdgeToCotan cotangents = computeCotangents(mesh);
    // 1b. Compute off-diagonal values
    
    // get neighbors list
    std::vector<std::vector<trimesh::index_t>> neighbors;
    for (size_t i = 0; i < mesh.Vertices.rows(); i++)
    {
        neighbors.push_back(mesh.vertex_vertex_neighbors(i));
    }
    

    for (size_t i = 0; i < mesh.Vertices.rows(); i++)
    {
        for (auto j: neighbors[i])
        {
            if (i == j)
                continue; // diagonal values are computed later

            trimesh::index_t edge = mesh.halfedge(mesh.directed_edge2he_index(i, j)).edge;
            auto cots = cotangents[edge];

            matrix_non_diagonal_elements.push_back(Triplet(i, j,
                                                           -0.5 * (cots.first + cots.second)));

        }
    }
    // use triplets
    Eigen::SparseMatrix<double> offDiagonals(mesh.Vertices.rows(), mesh.Vertices.rows());
    offDiagonals.setFromTriplets(matrix_non_diagonal_elements.begin(), matrix_non_diagonal_elements.end());

    // 2 Diagonals
    std::vector<Triplet> matrix_diagonal_elements;
    matrix_diagonal_elements.reserve(mesh.Vertices.rows());

    for (size_t i = 0; i < mesh.Vertices.rows(); i++)
    {
        double sum = 0;
        for (auto neighbor : neighbors[i])
        {
            sum += offDiagonals.coeff(i, neighbor);
        }
        matrix_diagonal_elements.push_back(Triplet(i,i,0.5*sum));
    }
    // use triplets
    laplacian.setFromTriplets(matrix_diagonal_elements.begin(), matrix_diagonal_elements.end());

    // add non diagonals
    laplacian = laplacian + offDiagonals;

    return laplacian;
}

// Computes the cotangents for every edge of the mesh
EdgeToCotan computeCotangents(const trimesh::trimesh_t &mesh)
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
                                       mesh.halfedge(halfedge).next_he));
        auto nextNextEdge = edgeLength(mesh,
                                       mesh.he_index2directed_edge(
                                           mesh.halfedge(mesh.halfedge(halfedge).next_he).next_he));
        double cosineAlpha = (-std::pow(initialEdgeLength, 2) + std::pow(nextEdge, 2) + std::pow(nextNextEdge, 2)) / (2 * nextEdge * nextNextEdge);
        auto cotAlpha = cosToCot(cosineAlpha);

        // beta
        auto oppositeNextEdge = edgeLength(mesh,
                                           mesh.he_index2directed_edge(
                                               mesh.halfedge(mesh.halfedge(halfedge).opposite_he).next_he));
        auto oppositeNextNextEdge = edgeLength(mesh,
                                               mesh.he_index2directed_edge(
                                                   mesh.halfedge(mesh.halfedge(mesh.halfedge(
                                                                                       halfedge)
                                                                                   .opposite_he)
                                                                     .next_he)
                                                       .next_he));
        double cosineBeta = (-std::pow(initialEdgeLength, 2) + std::pow(oppositeNextEdge, 2) + std::pow(oppositeNextNextEdge, 2)) / (2 * oppositeNextEdge * oppositeNextNextEdge);
        auto cotBeta = cosToCot(cosineBeta);

        cotans.push_back(std::make_pair(cotAlpha, cotBeta));
    }

    assert(cotans.size(), mesh.Edges.size());
    return cotans;
}

// Computes the length of the edge
double edgeLength(const trimesh::trimesh_t &mesh, const trimesh::edge_t &edge)
{
    Eigen::Vector3d vector = mesh.Vertices.row(edge.end()) - mesh.Vertices.row(edge.start());
    return vector.norm();
}

// Computes the length of the directed edge, unsigned result
double edgeLength(const trimesh::trimesh_t &mesh, const std::pair<trimesh::index_t, trimesh::index_t> &pair)
{
    Eigen::Vector3d vector = mesh.Vertices.row(pair.second) - mesh.Vertices.row(pair.first);
    return vector.norm();
}