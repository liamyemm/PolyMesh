ADD_LIBRARY(IntersectMesh OBJECT IntersectMesh/IntersectMesh.cpp)
ADD_LIBRARY(CutMesh OBJECT CutMesh/CutMesh.cpp)
ADD_LIBRARY(TestCase OBJECT HHO_Poisson/TestCase.cpp)

# Executables

INCLUDE_DIRECTORIES(HHO_Poisson)
INCLUDE_DIRECTORIES(HHO_Diffusion)
INCLUDE_DIRECTORIES(TransformMesh)

ADD_SUBDIRECTORY(HHO_Poisson)
ADD_SUBDIRECTORY(HHO_Diffusion)
ADD_SUBDIRECTORY(TransformMesh)
