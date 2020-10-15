#include <Eigen/Sparse>

#include "laplacian.h"
#include "trimesh.h"
#include "examine.h"
#include "iglViewer.h"
#include "eigenfunctions.h"

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

    auto results = eigenDecomposition(laplacian, numEigensToSee);

    // debug, examine
    examineEigenDecomposition(results.first, results.second);

    // get mass matrix
    Eigen::SparseMatrix<double> mass = computeMassMatrix(mesh);

    // debug, examine
    examineMass(mass);

    // igl viewer
    view(mesh, results.second);
}