FILE(GLOB source_list function.* basis.* common.*)
ADD_LIBRARY(function ${source_list})
TARGET_LINK_LIBRARIES(function math)
