#include "iglViewer.h"

#include <igl/opengl/glfw/Viewer.h>

void view(const trimesh::trimesh_t& mesh, const Eigen::MatrixXd& eigenvectors)
{
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