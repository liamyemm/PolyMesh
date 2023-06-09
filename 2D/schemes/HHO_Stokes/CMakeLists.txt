ADD_EXECUTABLE(hho_stokes $<TARGET_OBJECTS:StokesCore> $<TARGET_OBJECTS:StokesTests> $<TARGET_OBJECTS:TestCase> $<TARGET_OBJECTS:IntersectMesh> HHO_Stokes.cpp)

TARGET_LINK_LIBRARIES(hho_stokes ${Boost_LIBRARIES})
TARGET_LINK_LIBRARIES(hho_stokes math function mesh quadrature plot mesh_transform thread ${EXTERNAL_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(hho_stokes PUBLIC ${EIGEN3_INCLUDE_DIR})

