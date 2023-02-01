ADD_EXECUTABLE(hho_diffusion $<TARGET_OBJECTS:CutMesh> $<TARGET_OBJECTS:IntersectMesh> $<TARGET_OBJECTS:TestCase> hho_diffusion.cpp)

TARGET_LINK_LIBRARIES(hho_diffusion ${Boost_LIBRARIES})
TARGET_LINK_LIBRARIES(hho_diffusion function mesh quadrature plot hybridcore)
TARGET_INCLUDE_DIRECTORIES(hho_diffusion PUBLIC ${EIGEN3_INCLUDE_DIR})

