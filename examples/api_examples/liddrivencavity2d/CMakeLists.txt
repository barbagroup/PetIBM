cmake_minimum_required (VERSION 3.18)
project(liddrivencavity2d)

# path to petibm installation
set(PETIBM_DIR "$ENV{PETIBM_DIR}" CACHE PATH "The path to PetIBM installation.")
set(CMAKE_PREFIX_PATH ${PETIBM_DIR} ${CMAKE_PREFIX_PATH})

# needs explicitly use MPI
find_package(MPI REQUIRED)

# find petibm
message(CHECK_START "Finding petibm")
find_package(petibm 0.5.4 REQUIRED)
message(CHECK_PASS "done")

# add executable
add_executable(liddrivencavity2d main.cpp)

# dependencies to link against
target_link_libraries(liddrivencavity2d PRIVATE MPI::MPI_CXX petibm::petibm)

# copy the case configuration/post-processing files to the build folder
add_custom_command(
    TARGET liddrivencavity2d
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E
        copy ${CMAKE_CURRENT_SOURCE_DIR}/config.yaml ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND ${CMAKE_COMMAND} -E
        copy ${CMAKE_CURRENT_SOURCE_DIR}/postprocessing.xmf ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND ${CMAKE_COMMAND} -E
        copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/config ${CMAKE_CURRENT_BINARY_DIR}/config
)
