FILE(GLOB source_list HybridCore.*)
ADD_LIBRARY(hybridcore ${source_list})
TARGET_LINK_LIBRARIES(hybridcore mesh function quadrature thread)
