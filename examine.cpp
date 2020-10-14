#include "examine.h"

#include <Eigen/Sparse>
#include <iostream>
#include <fstream>

void examineLaplacian(const Eigen::SparseMatrix<double>& cotangents)
{
    // print a portion of the cotangent laplacian    
    std::ofstream myfile;
    myfile.open("matrix.txt");
    
    int limitRows = cotangents.rows() <= 12 ? cotangents.rows() : 12;
    int limitCols = cotangents.cols() <= 12 ? cotangents.cols() : 12;

    myfile << cotangents.block(0,0, limitRows, limitCols) << std::endl;

    // myfile << cotangents << std::endl;

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