# =====================================================================
# @file CMakeLists.txt
# @brief for cmake
# @author Pi-Yueh Chuang (pychuang@gwu.edu)
# @version alpha
# @date 2015-09-01
# =====================================================================


CMAKE_MINIMUM_REQUIRED(VERSION 2.8)


MESSAGE("===================================================================")
MESSAGE("Specifying CXX compiler and flags...\n")
MESSAGE("NOTE: clang++ is default C++ compiler. To use other compiler, ")
MESSAGE("please define the environment variable CXX\n")

IF (NOT DEFINED ENV{CXX})
    SET(CMAKE_CXX_COMPILER "clang++")
    MESSAGE("CXX compiler: clang++")
ELSE()
    SET(CMAKE_CXX_COMPILER $ENV{CXX})
    MESSAGE("CXX compiler: " ${CMAKE_CXX_COMPILER})
ENDIF()


SET(CMAKE_CXX_FLAGS "-std=c++11" CACHE STRING "")
SET(CMAKE_CXX_FLAGS_DEBUG "-g" CACHE STRING "")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "")



MESSAGE("\n===================================================================")
MESSAGE("Specifying C compiler and flags...\n")
MESSAGE("NOTE: clang is default C compiler. To use other compiler, ")
MESSAGE("please define the environment variable CXX\n")

IF (NOT DEFINED ENV{CC})
    SET(CMAKE_C_COMPILER "clang")
    MESSAGE("C compiler: clang")
ELSE()
    SET(CMAKE_C_COMPILER $ENV{CC})
    MESSAGE("C compiler: " ${CMAKE_C_COMPILER})
ENDIF()

SET(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "")
SET(CMAKE_C_FLAGS_DEBUG "-g" CACHE STRING "")



MESSAGE("\n===================================================================")
MESSAGE("Specifying other variables...")
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(EXECUTABLE_OUTPUT_PATH ./bin)

IF(NOT DEFINED CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING "")
ENDIF()

MESSAGE("-- Build type: " ${CMAKE_BUILD_TYPE})

SET(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")

MESSAGE("\n===================================================================")
MESSAGE("Starting generator...")
PROJECT(Poisson)


MESSAGE("\n===================================================================")
MESSAGE("Trying to find PETSc...")
MESSAGE("\nNOTE: If it fails to find PETSc, please set the environment variable")
MESSAGE("PETSC_DIR in order to notify CMake the location of PETSc\n")

LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake-modules")

FIND_PACKAGE(PETSc REQUIRED)

IF(${PETSc_FOUND})
    MESSAGE("\tPETSC_DIR: " ${PETSC_DIR})
    MESSAGE("\tPETSC_INCLUDES: " ${PETSC_INCLUDES})
    MESSAGE("\tPETSC_LIBRARIES: " ${PETSC_LIBRARIES})
ENDIF()



MESSAGE("\n===================================================================")
MESSAGE("Trying to find AmgX...")
MESSAGE("\nNOTE: If it fails to find AmgX, please set the environment variable")
MESSAGE("AMGX_DIR in order to notify CMake the location of AmgX\n")


MESSAGE("-- Finding amgx_c.h")
FIND_PATH(AMGX_INCLUDE_DIRS 
          NAMES amgx_c.h  
          PATHS ${AMGX_DIR} $ENV{AMGX_DIR} 
          PATH_SUFFIXES "include")

IF(${AMGX_INCLUDE_DIRS} STREQUAL "AMGX_INCLUDE_DIRS-NOTFOUND")
    MESSAGE(FATAL_ERROR "Could not find amgx_c.h!!!")
ELSE()
    MESSAGE("-- Finding amgx_c.h - found in " ${AMGX_INCLUDE_DIRS})
ENDIF()


MESSAGE("-- Finding amgx_capi.h")
FIND_FILE(AMGX_CAPI_FOUND 
          NAMES amgx_capi.h 
          PATHS ${AMGX_INCLUDE_DIRS} 
          NO_DEFAULT_PATH)

IF(${AMGX_CAPI_FOUND} STREQUAL "AMGX_CAPI_FOUND-NOTFOUND")
    MESSAGE(FATAL_ERROR "amgx_capi.h was not found in: " ${AMGX_INCLUDE_DIRS})
ELSE()
    MESSAGE("-- Finding amgx_capi.h - found")
ENDIF()


MESSAGE("-- Finding amgx_config.h")
FIND_FILE(AMGX_CONFIG_FOUND 
          NAMES amgx_config.h 
          PATHS ${AMGX_INCLUDE_DIRS} 
          NO_DEFAULT_PATH)

IF(${AMGX_CONFIG_FOUND} STREQUAL "AMGX_CONFIG_FOUND-NOTFOUND")
    MESSAGE(FATAL_ERROR "amgx_config.h was not found in: " ${AMGX_INCLUDE_DIRS})
ELSE()
    MESSAGE("-- Finding amgx_config.h - found")
ENDIF()


MESSAGE("-- Finding libamgxsh.so")
FIND_LIBRARY(AMGX_LIBRARIES
             NAMES amgxsh
             PATHS ${AMGX_DIR} $ENV{AMGX_DIR}
             PATH_SUFFIXES "lib")

IF(${AMGX_LIBRARIES} STREQUAL "AMGX_LIBRARIES-NOTFOUND")
    MESSAGE(FATAL_ERROR "libamgxsh.so was not found")
ELSE()
    MESSAGE("-- Finding libamgxsh.so - found: " ${AMGX_LIBRARIES})
ENDIF()



MESSAGE("\n===================================================================")
MESSAGE("Trying to find CUDA 6.5...")
MESSAGE("\nNOTE: If you have several different versions of CUDA, and if CMake ")
MESSAGE("can not find CUDA 6.5, please define CUDA_TOOLKIT_ROOT_DIR manually ")
MESSAGE("when calling CMake. For example:") 
MESSAGE("\n\t$ cmake . -DCUDA_TOOLKIT_ROOT_DIR=/opt/cuda65\n")

FIND_PACKAGE(CUDA 6.5 EXACT REQUIRED)

IF(${CUDA_FOUND})
    SET(CUDA_LIBRARY_DIRS ${CUDA_TOOLKIT_ROOT_DIR}/lib64)
    MESSAGE("\tCUDA_TOOLKIT_ROOT_DIR: " ${CUDA_TOOLKIT_ROOT_DIR})
    MESSAGE("\tCUDA_VERSION: " ${CUDA_VERSION_STRING})
    MESSAGE("\tCUDA_INCLUDE_DIRS: " ${CUDA_INCLUDE_DIRS})
    MESSAGE("\tCUDA_LIBRARY_DIRS: " ${CUDA_LIBRARY_DIRS})
    MESSAGE("\tCUDA_LIBRARIES: " ${CUDA_LIBRARIES})
ENDIF()


MESSAGE("\n===================================================================")
MESSAGE("Setting the information for target and source files...")


SET(SRC ${CMAKE_SOURCE_DIR}/src)
SET(AmgXWrapperSRC ${CMAKE_SOURCE_DIR}/../../src)

SET(SOURCE
    ${SRC}/Poisson.cpp 
    ${SRC}/generateA.cpp
    ${SRC}/generateGrid.cpp
    ${SRC}/generateRHS.cpp
    ${SRC}/generateExt.cpp
    ${SRC}/createKSP.cpp
    ${SRC}/solve.cpp
    ${SRC}/applyNeumannBC.cpp
    ${AmgXWrapperSRC}/AmgXSolver.cpp)

INCLUDE_DIRECTORIES(
    ${PETSC_INCLUDES} 
    ${AMGX_INCLUDE_DIRS} 
    ${CUDA_INCLUDE_DIRS} 
    ${AmgXWrapperSRC})

LINK_DIRECTORIES(
    ${CUDA_LIBRARY_DIRS} 
    ${AMGX_LIBRARY_DIRS}
    ${PETSC_LIBRARY_DIRS})

ADD_EXECUTABLE(Poisson ${SOURCE})

TARGET_LINK_LIBRARIES(Poisson 
    ${PETSC_LIBRARIES} 
    ${AMGX_LIBRARIES} 
    ${CUDA_LIBRARIES})

ADD_CUSTOM_COMMAND(
    TARGET Poisson
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory configs
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/configs configs/)
