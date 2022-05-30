add_executable(petibm-decoupledibpm decoupledibpm.h decoupledibpm.cpp main.cpp)

set_target_properties(petibm-decoupledibpm PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    INSTALL_RPATH "${PETIBM_RPATH}"
)

target_compile_options(petibm-decoupledibpm
    PRIVATE -Wno-deprecated-declarations
)

target_compile_definitions(petibm-decoupledibpm
    PRIVATE $<$<BOOL:${PETIBM_USE_AMGX}>:HAVE_AMGX>
)

target_link_libraries(petibm-decoupledibpm
    PRIVATE MPI::MPI_CXX
    PRIVATE PkgConfig::PETSC
    PRIVATE yaml-cpp
    PRIVATE petibm
    PRIVATE petibmapps
)

target_link_options(petibm-decoupledibpm
    PRIVATE "-Wl,--as-needed"
    PRIVATE "-Wl,--no-allow-shlib-undefined"
    PRIVATE "-Wl,--no-undefined"
)

install(
    TARGETS petibm-decoupledibpm
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
