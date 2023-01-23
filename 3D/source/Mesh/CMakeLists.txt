# Library: mesh
INCLUDE_DIRECTORIES(CurvedMesh)
INCLUDE_DIRECTORIES(StraightMesh)

ADD_SUBDIRECTORY(CurvedMesh)
ADD_SUBDIRECTORY(StraightMesh)

FILE(GLOB source_list Polytope.* Vertex.* Edge.* Face.* Cell.* Mesh.* mesh_reader.* mesh_builder.*)
ADD_LIBRARY(mesh ${source_list})
TARGET_LINK_LIBRARIES(mesh curved_mesh straight_mesh)

