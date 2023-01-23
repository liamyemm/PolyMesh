# Library: mesh

FILE(GLOB source_list CurvedPolytope.* CurvedVertex.* CurvedEdge.* CurvedFace.* CurvedCell.* CurvedMesh.*)
ADD_LIBRARY(curved_mesh ${source_list})
TARGET_LINK_LIBRARIES(curved_mesh math)

