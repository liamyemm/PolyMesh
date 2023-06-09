
FILE(GLOB source_list QuadratureRule.* generate_quadrature_rule.* QuadHandler.*)
ADD_LIBRARY(quadrature ${source_list})
TARGET_LINK_LIBRARIES(quadrature mesh function math)

ADD_EXECUTABLE(quad_test quadrature_test.cpp)

TARGET_LINK_LIBRARIES(quad_test ${Boost_LIBRARIES})
TARGET_LINK_LIBRARIES(quad_test function mesh quadrature math mesh_transform ${EXTERNAL_LIBRARIES})
