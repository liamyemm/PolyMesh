
FILE(GLOB source_list GaussLegendre.* QuadratureRule.*)
ADD_LIBRARY(quadrature ${source_list})
TARGET_LINK_LIBRARIES(quadrature mesh function)
