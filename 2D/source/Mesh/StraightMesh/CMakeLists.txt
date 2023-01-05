# Library: mesh

FILE(GLOB source_list StraightPolytope.* StraightVertex.* StraightEdge.* StraightCell.* StraightMesh.* MeshBuilder2D.* MeshReaderTyp2.*)
ADD_LIBRARY(straight_mesh ${source_list})
TARGET_LINK_LIBRARIES(straight_mesh function quadrature math)

