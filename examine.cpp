#include "examine.h"

#include <Eigen/Sparse>
#include <iostream>
#include <fstream>
#include <algorithm>

void examineLaplacian(const Eigen::SparseMatrix<double>& cotangents)
{
    // print a portion of the cotangent laplacian    
    std::ofstream myfile;
    myfile.open("matrix.txt");
    
    int limitRows = std::min(cotangents.rows(), 12);
    int limitCols = std::min(cotangents.cols(), 12);

    myfile << "Cotangents" << std::endl;
    myfile << cotangents.block(0,0, limitRows, limitCols) << std::endl;

    // myfile << cotangents << std::endl;

    // debug row sums
    // for (Eigen::Index i = 0; i < cotangents.rows(); i++)
    // {
    //     std::cout << cotangents.row(i).sum() << std::endl;
    // }

    myfile.close();
}

void examineEigenDecomposition(const Eigen::MatrixXd& eigenvectors,
    const Eigen::VectorXd& eigenvalues)
{
    // print a portion of the eigen decomposition
    std::ofstream myfile;
    myfile.open("matrix.txt", std::ios_base::app);
    
    int limitRows = eigenvectors.rows() <= 12 ? eigenvectors.rows() : 12;
    int limitCols = eigenvectors.cols() <= 12 ? eigenvectors.cols() : 12;

    myfile << std::endl << "Eigenvectors" << std::endl;
    myfile << eigenvectors.block(0,0, limitRows, limitCols) << std::endl;

    myfile << std::endl << "Eigenvalues" << std::endl;
    myfile << eigenvalues << std::endl;

    myfile.close();
}

void examineMass(const DiagonalXd& mass)
{
    // print a portion of the mass matrix
    std::ofstream myfile;
    myfile.open("matrix.txt", std::ios_base::app);
    
    myfile << std::endl << "Barycentric mass matrix" << std::endl;

    if (mass.rows() > 12)
    {
        
        myfile << mass.diagonal().head(12) << std::endl;
    }
    else
    {
        myfile << mass.diagonal() << std::endl;
    }
    
    myfile.close();
}

void examineTemperatures(const Eigen::MatrixXd& temperatures)
{
    // print a portion of the mass matrix
    std::ofstream myfile;
    myfile.open("matrix.txt", std::ios_base::app);
    
    myfile << std::endl << "Temperature cols" << std::endl;

    int limitRows = std::min(temperatures.rows(), 12);
    int limitCols = std::min(temperatures.cols(), 12);

    myfile << temperatures.block(0,0, limitRows, limitCols) << std::endl;
    
    myfile.close();
}