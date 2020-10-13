#include <igl/eigs.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <Eigen/Sparse>
#include <iostream>

#include "laplacian.h"
#include "trimesh.h"

int main(int argc, char * argv[])
{
    // build the mesh from the obj file
    trimesh::trimesh_t mesh("../../input/horse-1.obj");
    Eigen::MatrixXd eigenvectors;

    eigenvectors = mesh.Vertices.block(0,0,mesh.Vertices.rows(), 2);

    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(mesh.Vertices,mesh.Faces);
//    viewer.data().compute_normals();

    int selectedcolumn=0;
    viewer.callback_key_down = [&](igl::opengl::glfw::Viewer & viewer,unsigned char key,int)->bool
    {
        switch(key)
        {
            default:
                return false;
            case ' ':
            {
                selectedcolumn = (selectedcolumn + 1) % eigenvectors.cols();
                viewer.data().set_data(eigenvectors.col(selectedcolumn));
                return true;
            }
        }
    };
    viewer.callback_key_down(viewer,' ',0);
    viewer.data().show_lines = false;
    viewer.launch();
}