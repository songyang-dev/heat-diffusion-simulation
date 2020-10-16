#include "heatEquation.h"
#include "examine.h"

#include <Eigen/IterativeLinearSolvers>

#define TIME_INTERVAL 1
#define HEAT_COEFFICIENT 0.1

// simulates the surface temperature of the mesh after heat dissipates
Eigen::MatrixXd simulateHeat(size_t iterations,
    const Eigen::SparseMatrix<double>& cotangents,
    const DiagonalXd& mass,
    const Eigen::VectorXd& initialTemperature)
{
    assert(cotangents.rows() == mass.rows());
    assert(cotangents.cols() == mass.cols());
    assert(iterations > 0);
    
    Eigen::MatrixXd temperatures(mass.rows(), iterations);

    // Laplacian
    Eigen::SparseMatrix<double> laplacian = -1 * mass.inverse() * cotangents;

    // Initial conditions
    temperatures.col(0) = initialTemperature;

    // identity matrix
    Eigen::SparseMatrix<double> identity(laplacian.rows(), laplacian.cols());
    identity.setIdentity();
    
    // explicit discretization
    // for (size_t i = 1; i < iterations; i++)
    // {
    //     temperatures.col(i) = 
    //         (TIME_INTERVAL * HEAT_COEFFICIENT * laplacian + identity)
    //         * temperatures.col(i-1);
    // }

    // implicit discretization
    Eigen::SparseMatrix<double> coeffs(laplacian.rows(), laplacian.cols());
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    solver.compute(identity - TIME_INTERVAL * HEAT_COEFFICIENT * laplacian);
    for (size_t i = 1; i < iterations; i++)
    {
        temperatures.col(i) = solver.solve(temperatures.col(i-1));
    }
    
    
    // debug, examine
    examineTemperatures(temperatures);

    return temperatures;
}