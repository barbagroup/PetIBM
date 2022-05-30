enable_testing()

include(GoogleTest)

add_subdirectory(body)
add_subdirectory(boundary)
add_subdirectory(mesh)
add_subdirectory(misc)
add_subdirectory(operators)

add_custom_target(check
    COMMAND ${CMAKE_CTEST_COMMAND}
    DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Run tests"
    DEPENDS body boundary mesh misc operators
)
