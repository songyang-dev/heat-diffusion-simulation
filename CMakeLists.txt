cmake_minimum_required(VERSION 3.1)

# Windows cmake options for vcpkg
if(DEFINED ENV{VCPKG_ROOT} AND NOT DEFINED CMAKE_TOOLCHAIN_FILE)
  set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake"
      CACHE STRING "")
endif()

project(example)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)

find_package(LIBIGL REQUIRED QUIET)

set(
        main_src
        main.cpp
        trimesh.h
        trimesh.cpp
        trimesh_types.h
        laplacian.h
        laplacian.cpp
        utils.h
        utils.cpp
        examine.h
        examine.cpp
        iglViewer.h
        iglViewer.cpp
        eigenfunctions.h
        eigenfunctions.cpp
        heatEquation.h
        heatEquation.cpp
        diffusionEquation.h
        diffusionEquation.cpp
)

# C++ standards
set (CMAKE_CXX_STANDARD 17)

# Add your project files
add_executable(${PROJECT_NAME} ${main_src})
target_link_libraries(${PROJECT_NAME} igl::core igl::opengl_glfw)

