add_executable(decoupledibpm decoupledibpm.h decoupledibpm.cpp main.cpp)

set_target_properties(decoupledibpm PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    INSTALL_RPATH "${PETIBM_RPATH}"
)

target_compile_options(decoupledibpm
    PRIVATE -Wno-deprecated-declarations
)

target_compile_definitions(decoupledibpm
    PRIVATE $<$<BOOL:${PETIBM_USE_AMGX}>:HAVE_AMGX>
)

target_link_libraries(decoupledibpm
    PRIVATE MPI::MPI_CXX
    PRIVATE PkgConfig::PETSC
    PRIVATE yaml-cpp
    PRIVATE petibm
    PRIVATE petibmapps
)

target_link_options(decoupledibpm
    PRIVATE "-Wl,--as-needed"
    PRIVATE "-Wl,--no-allow-shlib-undefined"
    PRIVATE "-Wl,--no-undefined"
)

install(
    TARGETS decoupledibpm
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
