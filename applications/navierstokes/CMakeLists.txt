add_executable(navierstokes navierstokes.h navierstokes.cpp main.cpp)

set_target_properties(navierstokes PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    INSTALL_RPATH "${PETIBM_RPATH}"
)

target_compile_options(navierstokes
    PRIVATE -Wno-deprecated-declarations
)

target_compile_definitions(navierstokes
    PRIVATE $<$<BOOL:${PETIBM_USE_AMGX}>:HAVE_AMGX>
)

target_link_libraries(navierstokes
    PRIVATE MPI::MPI_CXX
    PRIVATE PkgConfig::PETSC
    PRIVATE yaml-cpp
    PRIVATE petibm
)

target_link_options(navierstokes
    PRIVATE "-Wl,--as-needed"
    PRIVATE "-Wl,--no-allow-shlib-undefined"
    PRIVATE "-Wl,--no-undefined"
)

install(
    TARGETS navierstokes
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
