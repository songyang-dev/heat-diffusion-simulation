/**
 * Computes the laplacian of a mesh
 * */

#include "laplacian.h"
#include "utils.h"

#include <Eigen/Sparse>

// Computes the cotangent Laplacian of the mesh
Eigen::SparseMatrix<double> computeCotangentLaplacian(const trimesh::trimesh_t &mesh)
{
    // to be filled with triplets
    Eigen::SparseMatrix<double> cotangentLaplacian(mesh.Vertices.rows(), mesh.Vertices.rows());

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
        for (auto j : neighbors[i])
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
        matrix_diagonal_elements.push_back(Triplet(i, i, -sum));
    }
    // use triplets
    cotangentLaplacian.setFromTriplets(matrix_diagonal_elements.begin(), matrix_diagonal_elements.end());

    // add non diagonals
    cotangentLaplacian = cotangentLaplacian + offDiagonals;

    return cotangentLaplacian;
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
            mesh.he_index2directed_edge(mesh.halfedge(halfedge).next_he));

        auto nextNextEdge = edgeLength(mesh,
            mesh.he_index2directed_edge(
                mesh.halfedge(mesh.halfedge(halfedge).next_he).next_he));
        double cosineAlpha = (-std::pow(initialEdgeLength, 2) + std::pow(nextEdge, 2) + std::pow(nextNextEdge, 2)) / (2 * nextEdge * nextNextEdge);
        auto cotAlpha = cosToCot(cosineAlpha);

        // beta
        auto oppositeNextEdge = edgeLength(mesh,
            mesh.he_index2directed_edge(
                mesh.halfedge(mesh.halfedge(halfedge).opposite_he).next_he)
        );
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

// Computes the mass matrix of the mesh, using barycentric lumped mass
Eigen::SparseMatrix<double> computeMassMatrix(const trimesh::trimesh_t& mesh)
{
    // 1. Loop on all the faces
    // 2. Find the area of each face
    // 3. Divide the area by 3 and add it to each vertex in the triangle
    // 4. Return a diagonal sparse matrix of each vertex area

    // empty vertex area array
    std::vector<double> areas(mesh.Vertices.size(), 0);

    // 1. Iterating on faces
    for (size_t i = 0; i < mesh.Faces.rows(); i++)
    {
        Eigen::Vector3i face = mesh.Faces.row(i);
        // 2. Area of each face
        auto area = faceArea(mesh, i);

        // 3. Split area
        areas[face.x()] += area/3;
        areas[face.y()] += area/3;
        areas[face.z()] += area/3;
    }
    

    // 4. Diagonal matrix
    // make triplets from the areas matrix
    std::vector<Triplet> triplets;
    triplets.reserve(areas.size());
    for (size_t i = 0; i < areas.size(); i++)
    {
        triplets.push_back(Triplet(i,i, areas[i]));
    }
    // create sparse matrix
    Eigen::SparseMatrix<double> massMatrix(triplets.size(), triplets.size());
    massMatrix.setFromTriplets(triplets.begin(), triplets.end());
    return massMatrix;
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

// Computes the area of the face at the given index
double faceArea(const trimesh::trimesh_t& mesh, const size_t& index){
    
    // use the cross product
    // A = 1/2 * ||P1P2 x P1P3||

    Eigen::Vector3i face = mesh.Faces.row(index);

    Eigen::Vector3d p1p2 = mesh.Vertices.row(face.y()) - mesh.Vertices.row(face.x());
    Eigen::Vector3d p1p3 = mesh.Vertices.row(face.z()) - mesh.Vertices.row(face.x());

    return p1p2.cross(p1p3).norm()/2;
}