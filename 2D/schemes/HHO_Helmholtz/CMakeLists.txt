ADD_EXECUTABLE(hho_helmholtz $<TARGET_OBJECTS:IntersectMesh> $<TARGET_OBJECTS:Enrichment> hho_helmholtz.cpp)

TARGET_LINK_LIBRARIES(hho_helmholtz ${Boost_LIBRARIES})
TARGET_LINK_LIBRARIES(hho_helmholtz function mesh quadrature plot hybridcore ${EXTERNAL_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(hho_helmholtz PUBLIC ${EIGEN3_INCLUDE_DIR})

