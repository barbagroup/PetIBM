# ==============================================================================
# \file CMakeLists.txt
# \brief for cmake
# \author Pi-Yueh Chuang (pychuang@gwu.edu)
# \date 2022-05-13
# ==============================================================================

cmake_minimum_required(VERSION 3.18.0)

# ==============================================================================
# project settings
# ==============================================================================
message(STATUS "Begin PetIBM configuration")

project(PetIBM
    VERSION 0.5.4
    HOMEPAGE_URL https://github.com/barbagroup/PetIBM
    LANGUAGES CXX
)

include(GNUInstallDirs)  # use GNU standard installation folder hierarchy
include(FetchContent)  # for downloading dependencies if needed
include(CMakePackageConfigHelpers)  # generate configs for downstream projects

# set default build type
if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY VALUE "Release")
endif()

# add options (i.e., BOOL-type CACHE variables)
option(BUILD_SHARED_LIBS "Whether to build shared library." ON)
option(PETIBM_USE_AMGX "Whether to build with GPU solver from AmgX." OFF)
option(PETIBM_ENABLE_TESTS "Whether to build tests." OFF)
option(PETIBM_BUILD_YAMLCPP "Whether to download and build yaml-cpp." OFF)
option(PETIBM_BUILD_AMGXWRAPPER "Whether to download and build AmgXWrapper." OFF)

# search-paths for dependencies (defaults determined from env vars)
set(PETSC_DIR "$ENV{PETSC_DIR}" CACHE PATH "The path to PETSc.")
set(PETSC_ARCH "$ENV{PETSC_ARCH}" CACHE STRING "The specific build of PETSc under PETSC_DIR.")
set(AMGXWRAPPER_DIR "$ENV{AMGXWRAPPER_DIR}" CACHE PATH "The path to AmgXWrapper.")
set(YAMLCPP_DIR "$ENV{YAMLCPP_DIR}" CACHE PATH "The path to yaml-cpp.")
set(SYMENGINE_DIR "$ENV{SYMENGINE_DIR}" CACHE PATH "The path to symengine.")

# search-paths when PETIBM_BUILD_AMGXWRAPPER is ON/TRUE
set(AMGX_DIR "$ENV{AMGX_DIR}" CACHE PATH "The path to AmgX.")
set(CUDA_DIR "$ENV{CUDA_DIR}" CACHE PATH "The path to CUDA.")

# hardcoded configuration
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(FETCHCONTENT_QUIET ON CACHE BOOL "" FORCE)

# combine PETSC_DIR and PETSC_ARCH into a single variable
set(PETSC_DIR "${PETSC_DIR}/${PETSC_ARCH}")

# configure extra paths to find cmake/pkg-config modules/configs
set(CMAKE_PREFIX_PATH ${YAMLCPP_DIR} ${CMAKE_PREFIX_PATH})
set(CMAKE_PREFIX_PATH ${SYMENGINE_DIR} ${CMAKE_PREFIX_PATH})
set(CMAKE_PREFIX_PATH ${PETSC_DIR} ${CMAKE_PREFIX_PATH})
set(CMAKE_PREFIX_PATH ${AMGXWRAPPER_DIR} ${CMAKE_PREFIX_PATH})
set(CMAKE_MODULE_PATH ${YAMLCPP_DIR} ${CMAKE_MODULE_PATH})
set(CMAKE_MODULE_PATH ${SYMENGINE_DIR} ${CMAKE_MODULE_PATH})
set(CMAKE_MODULE_PATH ${PETSC_DIR} ${CMAKE_MODULE_PATH})
set(CMAKE_MODULE_PATH ${AMGXWRAPPER_DIR} ${CMAKE_MODULE_PATH})

# ==============================================================================
# finding pkg-config and MPI
# ==============================================================================
find_package(PkgConfig REQUIRED)
find_package(MPI REQUIRED)
find_package(Doxygen)

# ==============================================================================
# finding petsc
# ==============================================================================
message(CHECK_START "Finding PETSc")
pkg_search_module(PETSC REQUIRED QUIET IMPORTED_TARGET petsc>=3.16)

message(CHECK_START "Checking PETSc version")
if (PETSC_VERSION VERSION_GREATER_EQUAL "3.17.0")
    message(FATAL_ERROR "PETSc >= 3.17 is not supported.")
endif()
message(CHECK_PASS "found PETSc ${PETSC_VERSION}")
message(CHECK_PASS "found at ${PETSC_PREFIX}")

# ==============================================================================
# finding yaml-cpp
# ==============================================================================
if (PETIBM_BUILD_YAMLCPP)  # download and build if an user requested
    message(CHECK_START "Downloading yaml-cpp 420c98")
    set(YAML_CPP_INSTALL ON CACHE BOOL "" FORCE)
    set(YAML_CPP_BUILD_TEST OFF CACHE BOOL "" FORCE)
    set(YAML_CPP_COMMIT 420c98231094b1cd2e5de3a714c4e3ee9b4f1118)
    FetchContent_Declare(
        _yamlcpp
        URL https://github.com/jbeder/yaml-cpp/tarball/${YAML_CPP_COMMIT}
    )
    message(CHECK_PASS "done")
    list(APPEND DOWNLOADED_DEPS "_yamlcpp")
    set(YAMLCPP_DIR ${CMAKE_INSTALL_PREFIX} CACHE PATH "" FORCE)
else()  # use user-provided yaml-cpp; use imported target to unify later code
    message(CHECK_START "Finding yaml-cpp")
    pkg_search_module(YAMLCPP REQUIRED  QUIET IMPORTED_TARGET yaml-cpp>=0.7.0)
    add_library(yaml-cpp UNKNOWN IMPORTED)
    set_target_properties(yaml-cpp PROPERTIES 
        NO_SONAME ON
        IMPORTED_LOCATION ${YAMLCPP_LINK_LIBRARIES}
        INTERFACE_INCLUDE_DIRECTORIES ${YAMLCPP_INCLUDE_DIRS}
    )
    message(CHECK_PASS "found yaml-cpp ${YAMLCPP_VERSION}")
endif()

# ==============================================================================
# finding symengine
# ==============================================================================
message(CHECK_START "Finding symengine")
find_package(SymEngine 0.9 REQUIRED)
message(CHECK_PASS "found symengine ${SymEngine_VERSION}")

# the imported target from SymEngine is defective. Patching it.
set_target_properties(symengine PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${SYMENGINE_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${SYMENGINE_LIBRARIES}"
)

# ==============================================================================
# finding amgxwrapper (amgxwrapper should include cuda and amgx)
# ==============================================================================
if (PETIBM_USE_AMGX)
    if (PETIBM_BUILD_AMGXWRAPPER)
        message(CHECK_START "Downloading amgxwrapper 1.6.1")
        FetchContent_Declare(
            _amgxwrapper
            GIT_REPOSITORY https://github.com/barbagroup/AmgXWrapper
            GIT_TAG v1.6.1
            GIT_SHALLOW ON
        )
        message(CHECK_PASS "done")
    list(APPEND DOWNLOADED_DEPS "_amgxwrapper")
    set(AMGXWRAPPER_DIR ${CMAKE_INSTALL_PREFIX} CACHE PATH "" FORCE)
    else()  # AmgXWrapper>=1.6.1 provides things CMake needs
        message(CHECK_START "Finding AmgXWrapper")
        find_package(amgxwrapper 1.6.1 REQUIRED)
        add_library(amgxwrapper ALIAS amgxwrapper::amgxwrapper)
        message(CHECK_PASS "found AmgXWrapper ${amgxwrapper_VERSION}")
    endif()
endif()

# ==============================================================================
# download google-test if tests are enabled (always download from GitHub)
# ==============================================================================
if (PETIBM_ENABLE_TESTS)
    message(CHECK_START "Downloading gtest 1.11.0")
    FetchContent_Declare(
        _googletest
        GIT_REPOSITORY https://github.com/google/googletest
        GIT_TAG release-1.11.0
        GIT_SHALLOW ON
    )
    message(CHECK_PASS "done")
    list(APPEND DOWNLOADED_DEPS "_googletest")
endif()

# ==============================================================================
# configure and populate targets for downloaded dependencies
# ==============================================================================
FetchContent_MakeAvailable(${DOWNLOADED_DEPS})

if (PETIBM_BUILD_YAMLCPP)
    target_compile_options(yaml-cpp PRIVATE -Wno-effc++)
endif()

# =============================================================================
# prepare runtime path list
# =============================================================================
# don't use CMAKE_INSTALL_RPATH because it may be overwritten by a parent
# project if PetIBM is included as a sub-project

# always include self' lib dir into rpath
list(APPEND PETIBM_RPATH ${CMAKE_INSTALL_FULL_LIBDIR})

# can't use imported target for PETSc because PkgConfig uses an interface
# target rather than an imported target
list(APPEND PETIBM_RPATH ${PETSC_LIBRARY_DIRS})

# include symengine lib dir into rpath
list(APPEND PETIBM_RPATH $<TARGET_LINKER_FILE_DIR:symengine>)

# only add yaml-cpp to rpath if using a pre-built package
# if using a downloaded package, TARGET_LINKER_FILE_DIR will be the build dir
# not the installation dir
list(APPEND PETIBM_RPATH
    $<$<NOT:$<BOOL:${PETIBM_BUILD_YAMLCPP}>>:$<TARGET_LINKER_FILE_DIR:yaml-cpp>>
)

# only add AmgXWrapper to rpath if using a pre-built package (line too long..)
string(APPEND PETIBM_RPATH  # note it's `string` append, not `list` append
    ";"
    "$<"
        "$<AND:"
            "$<BOOL:${PETIBM_USE_AMGX}>,"
            "$<NOT:$<BOOL:${PETIBM_BUILD_AMGXWRAPPER}>>"
        ">:"
        "$<TARGET_LINKER_FILE_DIR:amgxwrapper>"
    ">"
)

# =============================================================================
# targets
# =============================================================================
add_subdirectory(src)
add_subdirectory(applications)

if (PETIBM_ENABLE_TESTS)
    add_subdirectory(tests)
endif()

# =============================================================================
# examples
# =============================================================================
install(
    DIRECTORY examples
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/petibm
    PATTERN output EXCLUDE
)

# =============================================================================
# configuration files
# =============================================================================
# generate petibm-config.cmake for cmake's `find_package`
configure_package_config_file(
    petibm-config.cmake.in
    ${PROJECT_BINARY_DIR}/petibm-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/petibm
    PATH_VARS
        CMAKE_INSTALL_LIBDIR
        CMAKE_INSTALL_INCLUDEDIR
        PETSC_DIR
        YAMLCPP_DIR
        AMGXWRAPPER_DIR
)

# generate a petibm-config-version.cmake for cmake's `find_package`
write_basic_package_version_file(
    ${PROJECT_BINARY_DIR}/petibm-config-version.cmake
    COMPATIBILITY AnyNewerVersion
)

install(
    FILES
    ${PROJECT_BINARY_DIR}/petibm-config.cmake
    ${PROJECT_BINARY_DIR}/petibm-config-version.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/petibm
)
