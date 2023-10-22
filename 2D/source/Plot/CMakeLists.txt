# Library: plot

FILE(GLOB source_list vtk_writer.* vtu_writer.*)
ADD_LIBRARY(plot ${source_list})
TARGET_LINK_LIBRARIES(plot mesh ${VTK_LIBRARIES})
