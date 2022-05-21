add_executable(mesh
    cartesianmesh2d_dirichlet.cpp
    cartesianmesh2d_yperiodic.cpp
    cartesianmesh3d_dirichlet.cpp
    main.cpp
)

target_link_libraries(mesh
    PRIVATE MPI::MPI_CXX
    PRIVATE PkgConfig::PETSC
    PRIVATE yaml-cpp
    PRIVATE petibm
    PRIVATE gtest_main
)

cmake_path(GET CMAKE_CURRENT_BINARY_DIR PARENT_PATH TEST_DIR)
gtest_discover_tests(mesh WORKING_DIRECTORY ${TEST_DIR})
