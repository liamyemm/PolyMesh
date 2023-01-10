FILE(GLOB source_list HybridCore.* parallel_for.*)
ADD_LIBRARY(hybridcore ${source_list})
TARGET_LINK_LIBRARIES(hybridcore mesh function quadrature)
