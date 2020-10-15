#include <igl/eigs.h>
#include <Eigen/Sparse>

#include "laplacian.h"
#include "trimesh.h"
#include "examine.h"
#include "iglViewer.h"

int main(int argc, char * argv[])
{
    // build the mesh from the obj file
    trimesh::trimesh_t mesh(argc == 1 ? "../../input/horse-1.obj" : argv[1]);

    // # of eigenvectors to view
    int numEigensToSee = 1;
    if (argc == 3) numEigensToSee = atoi(argv[2]);

    // get laplacian
    Eigen::SparseMatrix<double> laplacian = computeCotangentLaplacian(mesh);

    // debug, examine
    examineLaplacian(laplacian);

    Eigen::SparseMatrix<double> identity(laplacian.rows(), laplacian.cols());
    identity.setIdentity();
    Eigen::MatrixXd eigenvectors;
    Eigen::VectorXd eigenvalues;
    //Eigen::SparseMatrix<double> flipped = laplacian * (-1);
    igl::eigs(laplacian, identity, numEigensToSee, igl::EIGS_TYPE_SM, eigenvectors, eigenvalues); //mesh.Vertices.block(0,0,mesh.Vertices.rows(), 2);

    // debug, examine
    examineEigenDecomposition(eigenvectors, eigenvalues);

    // get mass matrix
    Eigen::SparseMatrix<double> mass = computeMassMatrix(mesh);

    // debug, examine
    examineMass(mass);

    // igl viewer
    view(mesh, eigenvectors);
}