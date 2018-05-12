# =====================================================================
# \file CMakeLists.txt
# \brief for cmake
# \author Pi-Yueh Chuang (pychuang@gwu.edu)
# \date 2017-05-24
# =====================================================================


CMAKE_MINIMUM_REQUIRED(VERSION 3.0.2)


# =====================================================================
# Cache variables
# =====================================================================

# set default c++ compilers
IF (NOT DEFINED ENV{CXX})
    SET(ENV{CXX} mpicxx)
ENDIF()

# set default c compilers
IF (NOT DEFINED ENV{CC})
    SET(ENV{CC} mpicc)
ENDIF()

# set default flags for release and debug mode
SET(CMAKE_CXX_FLAGS_DEBUG "-g" CACHE STRING "")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "")
SET(CMAKE_C_FLAGS_DEBUG "-g" CACHE STRING "")
SET(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "")
SET(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined" CACHE STRING "")

# set default build type
SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING "")

SET(PETSC_DIR "$ENV{PETSC_DIR}" CACHE PATH "The path to PETSc.")
SET(PETSC_ARCH "$ENV{PETSC_ARCH}" CACHE STRING "The specific build of PETSc under PETSC_DIR.")
SET(CUDA_DIR "$ENV{CUDA_DIR}" CACHE PATH "The path to CUDA.")
SET(AMGX_DIR "$ENV{AMGX_DIR}" CACHE PATH "The path to AmgX.")


# =====================================================================
# Project settings
# =====================================================================

MESSAGE("")
MESSAGE("====================================")
MESSAGE("Begin configuration")
MESSAGE("====================================")

# start generator
PROJECT(AmgxWrapper C CXX)

# use GNU standard installation folder heirarchy
INCLUDE(GNUInstallDirs)

# set the sub folders in build folder
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})

# Not many vendors are supported by CMAKE_CXX_STANDARD variable before cmake 3.9
# so we decided to check C++11 flag manually
INCLUDE(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
IF(COMPILER_SUPPORTS_CXX11)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
ELSE()
    MESSAGE(FATAL_ERROR
        "Compiler ${CMAKE_CXX_COMPILER} does not support flag -std=c++11. \
        Please specify C++11 flag manually through CMAKE_CXX_FLAGS.")
ENDIF()

# other cmake setting
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_COLOR_MAKEFILE ON)

# add options
OPTION(BUILD_SHARED_LIBS "Whether to build shared library." ON)

# add cmake-modules to CMake module searching path
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake-modules")


# =====================================================================
# PETSc
# =====================================================================

MESSAGE("")
MESSAGE("-- Finding PETSc")

FIND_PACKAGE(PETSc REQUIRED)

IF(${PETSc_FOUND})
    IF (${PETSC_VERSION} VERSION_LESS "3.8")
        MESSAGE(FATAL_ERROR "ERROR: the version of PETSc should >= 3.8")
    ENDIF()
    MESSAGE("-- Finding PETSc - Success")
ENDIF()


# =====================================================================
# CUDA
# =====================================================================

MESSAGE("")
MESSAGE("-- Finding CUDA")

SET(CUDA_TOOLKIT_ROOT_DIR ${CUDA_DIR} CACHE PATH "" FORCE)

FIND_PACKAGE(CUDA REQUIRED)

IF(${CUDA_FOUND})
    SET(CUDA_LIBRARY_DIRS ${CUDA_TOOLKIT_ROOT_DIR}/lib64)
    MESSAGE("-- Finding CUDA - Success")
ENDIF()


# =====================================================================
# AmgX
# =====================================================================

MESSAGE("")
MESSAGE("-- Finding AmgX")

MESSAGE("-- Finding amgx_c.h")
FIND_PATH(AMGX_INCLUDE_DIRS 
          NAMES amgx_c.h  
          PATHS ${AMGX_DIR} $ENV{AMGX_DIR} 
          PATH_SUFFIXES "base/include" "examples")

IF(${AMGX_INCLUDE_DIRS} STREQUAL "AMGX_INCLUDE_DIRS-NOTFOUND")
    MESSAGE(FATAL_ERROR "Could not find amgx_c.h!!!")
ELSE()
    MESSAGE("-- Finding amgx_c.h - found in " ${AMGX_INCLUDE_DIRS})
ENDIF()


MESSAGE("-- Finding amgx_config.h")
FIND_FILE(AMGX_CONFIG_FOUND 
          NAMES amgx_config.h 
          PATHS ${AMGX_INCLUDE_DIRS} 
          NO_DEFAULT_PATH)

IF(${AMGX_CONFIG_FOUND} STREQUAL "AMGX_CONFIG_FOUND-NOTFOUND")
    MESSAGE(FATAL_ERROR "Could not find amgx_config.h in " ${AMGX_INCLUDE_DIRS})
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

MESSAGE("-- Finding AmgX - Success")


# =====================================================================
# Doxygen
# =====================================================================

MESSAGE("")
MESSAGE("-- Finding Doxygen")

FIND_PACKAGE(Doxygen)

IF (DOXYGEN_FOUND)
    MESSAGE("-- Finding Doxygen - Success")

    SET(DOXYGEN_INPUT ${CMAKE_SOURCE_DIR}/doc/Doxyfile)
    SET(DOXYGEN_OUTPUT ${PROJECT_BINARY_DIR}/doc)
    FILE(MAKE_DIRECTORY ${DOXYGEN_OUTPUT})
ELSE()
    MESSAGE("-- Finding Doxygen - Failed.")
ENDIF()

# =====================================================================
# Targets
# =====================================================================

SET(SRC ${CMAKE_SOURCE_DIR}/src)
SET(SOURCE
    ${SRC}/AmgXSolver.cpp
    ${SRC}/init.cpp
    ${SRC}/misc.cpp
    ${SRC}/setA.cpp
    ${SRC}/solve.cpp)

# target AmgXWrapper
ADD_LIBRARY(AmgXWrapper ${SOURCE})

TARGET_INCLUDE_DIRECTORIES(AmgXWrapper
    PRIVATE ${PETSC_INCLUDES}
    PRIVATE ${CUDA_INCLUDE_DIRS}
    PRIVATE ${AMGX_INCLUDE_DIRS}
    PUBLIC ${SRC})


TARGET_LINK_LIBRARIES(AmgXWrapper
    PRIVATE ${PETSC_LIBRARIES}
    PRIVATE ${AMGX_LIBRARIES}
    PRIVATE ${CUDA_LIBRARIES})

ADD_CUSTOM_COMMAND(
    TARGET AmgXWrapper
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/include
    COMMAND ${CMAKE_COMMAND} -E copy ${SRC}/AmgXSolver.hpp ${PROJECT_BINARY_DIR}/include)

# target doc
IF (DOXYGEN_FOUND)
    ADD_CUSTOM_TARGET(doc
        COMMAND ${CMAKE_COMMAND} -E echo "Building API Documentation..."
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT}
        COMMAND ${CMAKE_COMMAND} -E echo "Done."
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
ENDIF()

# installation
INSTALL(TARGETS AmgXWrapper
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

INSTALL(FILES ${SRC}/AmgXSolver.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
INSTALL(DIRECTORY ${PROJECT_BINARY_DIR}/doc DESTINATION ${CMAKE_INSTALL_PREFIX})

# =====================================================================
# Print information
# =====================================================================

MESSAGE("")
MESSAGE("====================================")
MESSAGE("Config Information:")
MESSAGE("====================================")
MESSAGE("")

MESSAGE("Build type: " ${CMAKE_BUILD_TYPE})
MESSAGE("Build shared library: " ${BUILD_SHARED_LIBS})
MESSAGE("Doxygen available: " ${DOXYGEN_FOUND})
MESSAGE("Installation path: " ${CMAKE_INSTALL_PREFIX})
MESSAGE("")

MESSAGE("Compilers:")
MESSAGE("\tCXX compiler: " ${CMAKE_CXX_COMPILER})
MESSAGE("\tC compiler: " ${CMAKE_C_COMPILER})
MESSAGE("")

MESSAGE("Flags:")
MESSAGE("\tCXX: " ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}})
MESSAGE("")

MESSAGE("PETSc:")
MESSAGE("\tPETSc dir: " ${PETSC_DIR})
MESSAGE("\tHeader path: " ${PETSC_INCLUDES})
MESSAGE("\tLibraries: " ${PETSC_LIBRARIES})
MESSAGE("")

MESSAGE("CUDA:")
MESSAGE("\tCUDA dir: " ${CUDA_TOOLKIT_ROOT_DIR})
MESSAGE("\tVersion: " ${CUDA_VERSION_STRING})
MESSAGE("\tHeader path: " ${CUDA_INCLUDE_DIRS})
MESSAGE("\tLibrary path: " ${CUDA_LIBRARY_DIRS})
MESSAGE("\tLibraries: " ${CUDA_LIBRARIES})
MESSAGE("")

MESSAGE("AmgX:")
MESSAGE("\tAmgX dir: " ${AMGX_DIR})
MESSAGE("\tHeader path: " ${AMGX_INCLUDE_DIRS})
MESSAGE("\tLibraries: " ${AMGX_LIBRARIES})
MESSAGE("")
