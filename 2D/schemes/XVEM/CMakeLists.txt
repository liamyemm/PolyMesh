ADD_EXECUTABLE(hho_poisson $<TARGET_OBJECTS:TestCase> hho_poisson.cpp)

TARGET_LINK_LIBRARIES(hho_poisson ${Boost_LIBRARIES})
TARGET_LINK_LIBRARIES(hho_poisson function mesh quadrature plot hybridcore)
TARGET_INCLUDE_DIRECTORIES(hho_poisson PUBLIC ${EIGEN3_INCLUDE_DIR})

