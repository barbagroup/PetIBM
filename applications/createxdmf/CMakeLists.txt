add_executable(petibm-createxdmf main.cpp)

set_target_properties(petibm-createxdmf PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    INSTALL_RPATH "${PETIBM_RPATH}"
)

target_compile_options(petibm-createxdmf
    PRIVATE -Wno-deprecated-declarations
)

target_compile_definitions(petibm-createxdmf
    PRIVATE $<$<BOOL:${PETIBM_USE_AMGX}>:HAVE_AMGX>
)

target_link_libraries(petibm-createxdmf
    PRIVATE MPI::MPI_CXX
    PRIVATE PkgConfig::PETSC
    PRIVATE yaml-cpp
    PRIVATE petibm
)

target_link_options(petibm-createxdmf
    PRIVATE "-Wl,--as-needed"
    PRIVATE "-Wl,--no-allow-shlib-undefined"
    PRIVATE "-Wl,--no-undefined"
)

install(
    TARGETS petibm-createxdmf
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
