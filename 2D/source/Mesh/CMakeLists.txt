# Library: mesh

FILE(GLOB source_list Polytope.* Vertex.* Edge.* Cell.* Mesh.* MeshBuilder2D.* MeshReaderTyp2.*)
ADD_LIBRARY(mesh ${source_list})
TARGET_LINK_LIBRARIES(mesh function quadrature math)

