add_executable(misc delta_test.cpp)

target_link_libraries(misc
    PRIVATE MPI::MPI_CXX
    PRIVATE PkgConfig::PETSC
    PRIVATE yaml-cpp
    PRIVATE petibm
    PRIVATE gtest_main
)

cmake_path(GET CMAKE_CURRENT_BINARY_DIR PARENT_PATH TEST_DIR)
gtest_discover_tests(misc WORKING_DIRECTORY ${TEST_DIR})
