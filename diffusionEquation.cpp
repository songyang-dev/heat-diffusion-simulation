#include "diffusionEquation.h"

#define TIME_INTERVAL 0.001
#define SMOOTH_COEFFICIENT 0.1

// Computes the diffusion equation "smoothing" of a mesh
std::vector<Eigen::MatrixXd> simulateDiffusion(size_t iterations,
    const Eigen::SparseMatrix<double>& cotangents,
    const DiagonalXd& mass,
    const trimesh::trimesh_t& initialMesh)
{
    std::vector<Eigen::MatrixXd> smoothed;

    // Form 3 linear systems to solve independently
    // Each linear system solves a coordinate component of the mesh points

    // Laplacian
    Eigen::SparseMatrix<double> laplacian = -1 * mass.inverse() * cotangents;

    // identity matrix
    Eigen::SparseMatrix<double> identity(laplacian.rows(), laplacian.cols());
    identity.setIdentity();

    // BICGSTAB solver
    Eigen::SparseMatrix<double> coeffs(laplacian.rows(), laplacian.cols());
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    solver.compute(identity - TIME_INTERVAL * SMOOTH_COEFFICIENT * laplacian);

    // initial conditions
    smoothed.push_back(initialMesh.Vertices);

    // iterate
    for (size_t i = 1; i < iterations; i++)
    {
        // Smoothed vertices
        Eigen::MatrixXd vertices(initialMesh.Vertices.rows(), 3);

        // solve for each component

        // x components
        vertices.col(0) = solver.solve(smoothed[i-1].col(0));
        // y components
        vertices.col(1) = solver.solve(smoothed[i-1].col(1));
        // z components
        vertices.col(2) = solver.solve(smoothed[i-1].col(2));

        smoothed.push_back(vertices);
    }
    
    return smoothed;
}