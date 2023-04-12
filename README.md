# PolyMesh Library

The PolyMesh library is a C++ library for implementing polytopal methods for partial differential equations, with a particular focus on hybrid high-order methods. It allows for the use of very general and possibly curved meshes and for enriched polynomial spaces.

For more information on using the PolyMesh library, see the documentation at https://liamyemm.github.io/PolyMesh/.

## Dependencies

The PolyMesh library depends on the following libraries:

- Eigen library
- Boost library
- CMake

## Installation

To install the PolyMesh library, follow these steps:

1. Clone the PolyMesh library source code from GitHub:

git clone git@github.com:liamyemm/PolyMesh.git


2. Navigate to the PolyMesh directory:

cd PolyMesh

3. Navigate to the directory for the 2D version (or 3D code):

cd 2D


3. Create a build directory:

mkdir build


4. Navigate to the build directory:

cd build


5. Run CMake to generate the build files:

cmake ..

6. Compile the library:

make
