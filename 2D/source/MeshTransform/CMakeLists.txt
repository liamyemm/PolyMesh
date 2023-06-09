# Library: mesh_transform

FILE(GLOB source_list affine_shift.* mesh_convert.*)
ADD_LIBRARY(mesh_transform ${source_list})
TARGET_LINK_LIBRARIES(mesh_transform math mesh function)

