#include "heatEquation.h"
#include "examine.h"

#define TIME_INTERVAL 0.5
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
    Eigen::SparseMatrix<double> laplacian = mass.inverse() * cotangents;

    // Initial conditions
    temperatures.col(0) = initialTemperature;

    // explicit discretization
    for (size_t i = 1; i < iterations; i++)
    {
        temperatures.col(i) = 
            -TIME_INTERVAL * HEAT_COEFFICIENT * laplacian * temperatures.col(i-1)
            + temperatures.col(i-1);
    }
    
    // debug, examine
    examineTemperatures(temperatures);

    return temperatures;
}