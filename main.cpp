#include <Eigen/Sparse>

#include "laplacian.h"
#include "trimesh.h"
#include "examine.h"
#include "iglViewer.h"
#include "eigenfunctions.h"
#include "heatEquation.h"
#include "diffusionEquation.h"

int main(int argc, char * argv[])
{
    // build the mesh from the obj file
    trimesh::trimesh_t mesh(argc == 1 ? "../../input/horse-1.obj" : argv[1]);

    // # of eigenvectors to view
    int numEigensToSee = 1;
    if (argc == 3) numEigensToSee = atoi(argv[2]);

    // get cotangents
    Eigen::SparseMatrix<double> cotangents = computeCotangentLaplacian(mesh);

    // debug, examine
    examineLaplacian(cotangents);

    auto results = eigenDecomposition(cotangents, numEigensToSee);

    // debug, examine
    examineEigenDecomposition(results.first, results.second);

    // get mass matrix
    auto mass = computeMassMatrix(mesh);

    // debug, examine
    examineMass(mass);

    // UNCOMMENT HERE AS NEEDED!

    // heat equation simulation
    // Eigen::VectorXd initial(mesh.Vertices.rows());
    // initial.setZero();
    // initial[0] = 1;
    // auto heat = simulateHeat(10, cotangents, mass, results.first.col(0)); // eigenvector conditions
    // auto heat = simulateHeat(10, cotangents, mass, initial); // one hat initial conditions

    // igl viewer
    // view(mesh, results.first); // eigenvectors
    // view(mesh, heat); // heat simulation
    viewSmoothing(simulateDiffusion(10, cotangents, mass, mesh), mesh.Faces); // mesh smoothing
}