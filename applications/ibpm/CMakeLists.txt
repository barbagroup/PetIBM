add_executable(ibpm ibpm.h ibpm.cpp main.cpp)

set_target_properties(ibpm PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    INSTALL_RPATH "${PETIBM_RPATH}"
)

target_compile_options(ibpm
    PRIVATE -Wno-deprecated-declarations
)

target_compile_definitions(ibpm
    PRIVATE $<$<BOOL:${PETIBM_USE_AMGX}>:HAVE_AMGX>
)

target_link_libraries(ibpm
    PRIVATE MPI::MPI_CXX
    PRIVATE PkgConfig::PETSC
    PRIVATE yaml-cpp
    PRIVATE petibm
    PRIVATE petibmapps
)

target_link_options(ibpm
    PRIVATE "-Wl,--as-needed"
    PRIVATE "-Wl,--no-allow-shlib-undefined"
    PRIVATE "-Wl,--no-undefined"
)

install(
    TARGETS ibpm
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
