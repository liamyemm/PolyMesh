project(PolyMesh)
CMAKE_MINIMUM_REQUIRED(VERSION 3.5)

#------------------------------------------------------------------------------#
#                             Compiler setup
#------------------------------------------------------------------------------#

IF(CMAKE_SYSTEM_NAME=Windows)
  SET(ADDITIONAL_FLAGS "-DWIN32")
ENDIF(CMAKE_SYSTEM_NAME=Windows)

#SET(CMAKE_CXX_COMPILER "c++")
SET(CMAKE_CXX_COMPILER "clang++")
SET(CMAKE_CXX_FLAGS "-std=c++20 -g -Wall")
#SET(CMAKE_CXX_FLAGS "-std=c++20 -g -Wall -Wno-maybe-uninitialized")
SET(CMAKE_CXX_FLAGS_DEBUGALL "${CMAKE_CXX_FLAGS} -fsanitize=address,undefined -pthread" )
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O3 -DOPTIMIZE -pthread" ${ADDITIONAL_FLAGS} )
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -DOPTIMIZE -pthread" ${ADDITIONAL_FLAGS} )
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} -O3 -g -DNDEBUG -DOPTIMIZE -pthread" ${ADDITIONAL_FLAGS} )

SET(CMAKE_BUILD_TYPE Debug)
#SET(CMAKE_BUILD_TYPE Debugall)
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build: Debug Release RelWithDebInfo" FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

#------------------------------------------------------------------------------#
#                             Required libraries
#------------------------------------------------------------------------------#

SET(Boost_DEBUG "ON")
SET(Boost_USE_STATIC_LIBS "ON")
SET(HOME_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
SET(CMAKE_MODULE_PATH ${HOME_DIR}/CMake)

FIND_PACKAGE(Boost COMPONENTS filesystem program_options timer chrono REQUIRED)
SET(EXTERNAL_INCLUDES ${Boost_INCLUDES})
SET(EXTERNAL_LIBRARIES ${Boost_LIBRARIES})

FIND_PACKAGE(Eigen3 3.3.4)
SET(EXTERNAL_INCLUDES ${EXTERNAL_INCLUDES} ${EIGEN3_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIRS})

FIND_PACKAGE(MKL)
IF(MKL_FOUND)
  SET(EXTERNAL_LIBRARIES ${EXTERNAL_LIBRARIES} ${MKL_LIBRARIES})
  SET(EXTERNAL_INCLUDES ${EXTERNAL_INCLUDES} ${MKL_INCLUDES})
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_MKL")  
ENDIF()

INCLUDE_DIRECTORIES(${EXTERNAL_INCLUDES})

#------------------------------------------------------------------------------#
#                             Directories for compiled libraries
#------------------------------------------------------------------------------#

INCLUDE_DIRECTORIES(source/Function)
INCLUDE_DIRECTORIES(source/Quadrature)
INCLUDE_DIRECTORIES(source/Mesh)
INCLUDE_DIRECTORIES(source/Plot)
INCLUDE_DIRECTORIES(source/HybridCore)
INCLUDE_DIRECTORIES(source/Math)
INCLUDE_DIRECTORIES(source/MeshTransform)
INCLUDE_DIRECTORIES(source/Thread)

ADD_SUBDIRECTORY(source/Function)
ADD_SUBDIRECTORY(source/Quadrature)
ADD_SUBDIRECTORY(source/Mesh)
ADD_SUBDIRECTORY(source/Plot)
ADD_SUBDIRECTORY(source/HybridCore)
ADD_SUBDIRECTORY(source/Math)
ADD_SUBDIRECTORY(source/MeshTransform)
ADD_SUBDIRECTORY(source/Thread)

#------------------------------------------------------------------------------#
#                             Directories for schemes
#------------------------------------------------------------------------------#

INCLUDE_DIRECTORIES(schemes)
ADD_SUBDIRECTORY(schemes)

