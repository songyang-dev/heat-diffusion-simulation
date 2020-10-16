# IFT 6113 HW 2
Implements Laplacian computations, eigenvector visualization, heat simulation and diffusion equation smoothing of a mesh.

# Usage
The default configuration is to run the diffusion equation smoothing of a mesh.
`./example.exe name-of-obj`

If you want to run the heat simulation or the eigenvector colorings, you need to go into `main.cpp` and uncomment the various function calls. After that, recompile and 
`./example.exe name-of-obj [number-of-eigenvectors-to-get]`
will try to solve the given amount of Laplacian eigenvectors of the mesh. Default is 1. 

Press Space to interact in the viewer, here is the behavior:
* Eigenvector colorings: cycle through colorings, in order of frequency
* Heat simulation: advance through time, winds back after 10 iterations
* Smoothing simulation: increase the smoothing level, winds back after 10 iterations

# Installation
## Prerequisites
For Windows, the following c++ packages must be install with `vcpkg` (https://github.com/Microsoft/vcpkg):
* eigen3
* glfw
* glad
* opengl

Manual installation of `libigl` is required, as elaborated by the [template repository](https://github.com/ivanpuhachov/ift6113_2020/tree/master/hw1_cpp).

*Note*: There is no additional dependency in this project compared to the template repository. So the project should be able to compile on any system that satisfies the template repo.

## Build
### Windows
```batch
mkdir build
cd build
..\..\build.bat release
```
### Other OS
Follow the steps on the template repository.