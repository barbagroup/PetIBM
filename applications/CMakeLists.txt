# library of developing applications
add_library(petibmapps
    navierstokes/navierstokes.cpp
    navierstokes/navierstokes.h
    ibpm/ibpm.cpp
    ibpm/ibpm.h
    decoupledibpm/decoupledibpm.cpp
    decoupledibpm/decoupledibpm.h
    rigidkinematics/rigidkinematics.cpp
    rigidkinematics/rigidkinematics.h
)

set_target_properties(petibmapps PROPERTIES
    POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS}
    INSTALL_RPATH "${PETIBM_RPATH}"
)

target_include_directories(petibmapps
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_compile_options(petibmapps
    PRIVATE -Wno-deprecated-declarations
)

target_compile_definitions(petibmapps
    PRIVATE $<$<BOOL:${PETIBM_USE_AMGX}>:HAVE_AMGX>
)

target_link_libraries(petibmapps
    PRIVATE MPI::MPI_CXX
    PRIVATE PkgConfig::PETSC
    PUBLIC yaml-cpp
    PUBLIC petibm
)

target_link_options(petibmapps
    PRIVATE "-Wl,--as-needed"
    PRIVATE "-Wl,--no-allow-shlib-undefined"
    PRIVATE "-Wl,--no-undefined"
)

install(
    TARGETS petibmapps
    EXPORT PetIBMAppsTarget
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(
    DIRECTORY
        ${CMAKE_CURRENT_SOURCE_DIR}/navierstokes
        ${CMAKE_CURRENT_SOURCE_DIR}/ibpm
        ${CMAKE_CURRENT_SOURCE_DIR}/decoupledibpm
        ${CMAKE_CURRENT_SOURCE_DIR}/rigidkinematics
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/petibm
    FILES_MATCHING PATTERN "*.h"
)

# executables
add_subdirectory(createxdmf)
add_subdirectory(writemesh)
add_subdirectory(vorticity)
add_subdirectory(navierstokes)
add_subdirectory(ibpm)
add_subdirectory(decoupledibpm)
