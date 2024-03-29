#include "iglViewer.h"

#include <igl/opengl/glfw/Viewer.h>

void view(const trimesh::trimesh_t& mesh, const Eigen::MatrixXd& coloringVectors)
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
                selectedcolumn = (selectedcolumn + 1) % coloringVectors.cols();
                viewer.data().set_data(coloringVectors.col(selectedcolumn));
                return true;
            }
        }
    };
    viewer.callback_key_down(viewer,' ',0);
    viewer.data().show_lines = false;
    viewer.launch();
}

void viewSmoothing(const std::vector<Eigen::MatrixXd>& vertices, 
    const Eigen::MatrixXi& faces)
{
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(vertices[0], faces);
//    viewer.data().compute_normals();

    int selectedSmoothingLevel=0;
    viewer.callback_key_down = [&](igl::opengl::glfw::Viewer & viewer,unsigned char key,int)->bool
    {
        switch(key)
        {
            default:
                return false;
            case ' ':
            {
                viewer.data().set_vertices(vertices[selectedSmoothingLevel]);
                selectedSmoothingLevel = (selectedSmoothingLevel + 1) % vertices.size();
                return true;
            }
        }
    };
    viewer.callback_key_down(viewer,' ',0);
    viewer.data().show_lines = false;
    viewer.launch();
}