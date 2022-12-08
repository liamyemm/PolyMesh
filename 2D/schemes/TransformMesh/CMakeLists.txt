# Library: mesh

#FILE(GLOB source_list vertex.* curved_*)
#ADD_LIBRARY(curved_mesh ${source_list})


ADD_EXECUTABLE(transform_mesh TransformMesh.cpp)
TARGET_LINK_LIBRARIES(transform_mesh mesh)
TARGET_INCLUDE_DIRECTORIES(transform_mesh PUBLIC ${EIGEN3_INCLUDE_DIR})

