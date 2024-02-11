ADD_LIBRARY(xvem_core OBJECT XVEMCore.cpp)
ADD_LIBRARY(poisson_singularity OBJECT PoissonSingularity.cpp)

ADD_EXECUTABLE(xvem xvem_poisson.cpp)
ADD_EXECUTABLE(xvem_core_test Tests/XVEMCoreTest.cpp)
ADD_EXECUTABLE(poisson_singularity_test Tests/PoissonSingularityTest.cpp)

TARGET_LINK_LIBRARIES(poisson_singularity_test ${Boost_LIBRARIES} ${GTEST_BOTH_LIBRARIES} math function poisson_singularity)
TARGET_INCLUDE_DIRECTORIES(poisson_singularity_test PUBLIC ${EIGEN3_INCLUDE_DIR})

TARGET_LINK_LIBRARIES(xvem_core_test ${Boost_LIBRARIES} ${GTEST_BOTH_LIBRARIES} function math mesh mesh_transform thread quadrature xvem_core poisson_singularity)
TARGET_INCLUDE_DIRECTORIES(xvem_core_test PUBLIC ${EIGEN3_INCLUDE_DIR})

TARGET_LINK_LIBRARIES(xvem ${Boost_LIBRARIES} function math mesh mesh_transform thread quadrature xvem_core poisson_singularity plot)
TARGET_INCLUDE_DIRECTORIES(xvem PUBLIC ${EIGEN3_INCLUDE_DIR})

