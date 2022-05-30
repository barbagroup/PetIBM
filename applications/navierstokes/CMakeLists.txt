add_executable(petibm-navierstokes navierstokes.h navierstokes.cpp main.cpp)

set_target_properties(petibm-navierstokes PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    INSTALL_RPATH "${PETIBM_RPATH}"
)

target_compile_options(petibm-navierstokes
    PRIVATE -Wno-deprecated-declarations
)

target_compile_definitions(petibm-navierstokes
    PRIVATE $<$<BOOL:${PETIBM_USE_AMGX}>:HAVE_AMGX>
)

target_link_libraries(petibm-navierstokes
    PRIVATE MPI::MPI_CXX
    PRIVATE PkgConfig::PETSC
    PRIVATE yaml-cpp
    PRIVATE petibm
)

target_link_options(petibm-navierstokes
    PRIVATE "-Wl,--as-needed"
    PRIVATE "-Wl,--no-allow-shlib-undefined"
    PRIVATE "-Wl,--no-undefined"
)

install(
    TARGETS petibm-navierstokes
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
