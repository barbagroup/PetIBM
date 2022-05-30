add_executable(petibm-vorticity main.cpp)

set_target_properties(petibm-vorticity PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    INSTALL_RPATH "${PETIBM_RPATH}"
)

target_compile_options(petibm-vorticity
    PRIVATE -Wno-deprecated-declarations
)

target_compile_definitions(petibm-vorticity
    PRIVATE $<$<BOOL:${PETIBM_USE_AMGX}>:HAVE_AMGX>
)

target_link_libraries(petibm-vorticity
    PRIVATE MPI::MPI_CXX
    PRIVATE PkgConfig::PETSC
    PRIVATE yaml-cpp
    PRIVATE petibm
)

target_link_options(petibm-vorticity
    PRIVATE "-Wl,--as-needed"
    PRIVATE "-Wl,--no-allow-shlib-undefined"
    PRIVATE "-Wl,--no-undefined"
)

install(
    TARGETS petibm-vorticity
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
