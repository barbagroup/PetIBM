set(SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/private/private.h
    ${CMAKE_CURRENT_SOURCE_DIR}/body/bodypack.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/body/singlebody.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/body/singlebodypoints.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/boundary/boundary.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/boundary/boundarysimple.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/boundary/singleboundary.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/boundary/singleboundarybase.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/boundary/singleboundaryconvective.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/boundary/singleboundarydirichlet.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/boundary/singleboundaryneumann.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/boundary/singleboundaryperiodic.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/io/io.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/linsolver/linsolverksp.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/linsolver/linsolver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mesh/cartesianmesh.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mesh/mesh.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/misc/delta.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/misc/lininterp.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/misc/misc.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/misc/probes.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/misc/type.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/operators/createbn.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/operators/createconvection.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/operators/createdelta.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/operators/creatediagmatrix.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/operators/createdivergence.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/operators/creategradient.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/operators/createlaplacian.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/parser/parser.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/solution/solution.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/solution/solutionsimple.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/timeintegration/timeintegration.cpp
    $<$<BOOL:${PETIBM_USE_AMGX}>:${CMAKE_CURRENT_SOURCE_DIR}/linsolver/linsolveramgx.cpp>
)

# add library target
add_library(petibm ${SRC})

set_target_properties(petibm PROPERTIES
    POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS}
    INSTALL_RPATH "${PETIBM_RPATH}"
)

target_include_directories(petibm
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_compile_options(petibm
    PRIVATE -Wno-deprecated-declarations
)

target_compile_definitions(petibm
    PRIVATE $<$<BOOL:${PETIBM_USE_AMGX}>:HAVE_AMGX>
)

target_link_libraries(petibm
    PRIVATE MPI::MPI_CXX
    PRIVATE PkgConfig::PETSC
    PUBLIC yaml-cpp
    PUBLIC $<$<BOOL:${PETIBM_USE_AMGX}>:AmgXWrapper>
)

target_link_options(petibm
    PRIVATE "-Wl,--as-needed"
    PRIVATE "-Wl,--no-allow-shlib-undefined"
    PRIVATE "-Wl,--no-undefined"
)

install(
    TARGETS petibm
    EXPORT PetIBMTarget
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
    DIRECTORY ${PROJECT_SOURCE_DIR}/include/petibm
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
