ADD_EXECUTABLE(hho_poisson hho_poisson.cpp)

TARGET_LINK_LIBRARIES(hho_poisson ${Boost_LIBRARIES})
TARGET_LINK_LIBRARIES(hho_poisson mesh quadrature)
TARGET_INCLUDE_DIRECTORIES(hho_poisson PUBLIC ${EIGEN3_INCLUDE_DIR})

