#include "eigenfunctions.h"

#include <igl/eigs.h>

std::pair<Eigen::MatrixXd, Eigen::VectorXd> eigenDecomposition(const Eigen::SparseMatrix<double>& matrix, int numEigensToSee)
{
    Eigen::SparseMatrix<double> identity(matrix.rows(), matrix.cols());
    identity.setIdentity();
    Eigen::MatrixXd eigenvectors;
    Eigen::VectorXd eigenvalues;
    //Eigen::SparseMatrix<double> flipped = matrix * (-1);
    igl::eigs(matrix, identity, numEigensToSee, igl::EIGS_TYPE_SM, eigenvectors, eigenvalues);

    return std::make_pair(eigenvectors, eigenvalues);
}