add_executable(body singlebody_test.cpp)

target_link_libraries(body
    PRIVATE MPI::MPI_CXX
    PRIVATE PkgConfig::PETSC
    PRIVATE yaml-cpp
    PRIVATE petibm
    PRIVATE gtest_main
)

add_custom_command(
    TARGET body
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E
        copy ${CMAKE_CURRENT_SOURCE_DIR}/body2d.txt ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND ${CMAKE_COMMAND} -E
        copy ${CMAKE_CURRENT_SOURCE_DIR}/body3d.txt ${CMAKE_CURRENT_BINARY_DIR}
)

cmake_path(GET CMAKE_CURRENT_BINARY_DIR PARENT_PATH TEST_DIR)
gtest_discover_tests(body WORKING_DIRECTORY ${TEST_DIR})
