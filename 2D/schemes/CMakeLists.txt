ADD_LIBRARY(IntersectMesh OBJECT IntersectMesh/IntersectMesh.cpp)
ADD_LIBRARY(IntersectCurvedMesh OBJECT IntersectCurvedMesh/IntersectCurvedMesh.cpp)
ADD_LIBRARY(TestCase OBJECT HHO_Poisson/TestCase.cpp)
ADD_LIBRARY(Enrichment OBJECT HHO_Helmholtz/Enrichment.cpp)

# Executables

INCLUDE_DIRECTORIES(HHO_Poisson)
INCLUDE_DIRECTORIES(HHO_Diffusion)
INCLUDE_DIRECTORIES(HHO_Helmholtz)
#INCLUDE_DIRECTORIES(HHO_Heat)
INCLUDE_DIRECTORIES(TransformMesh)

ADD_SUBDIRECTORY(HHO_Poisson)
ADD_SUBDIRECTORY(HHO_Diffusion)
ADD_SUBDIRECTORY(HHO_Helmholtz)
#ADD_SUBDIRECTORY(HHO_Heat)
ADD_SUBDIRECTORY(TransformMesh)
