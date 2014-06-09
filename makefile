ALL: bin/PetIBM
DIRS       = src src/include src/solvers
YAMLOBJ    = external/yaml-cpp/src/*.o external/yaml-cpp/src/contrib/*.o
LIBS       = lib/libclasses.a lib/libyaml.a lib/libsolvers.a
SRC        = ${wildcard src/*.cpp}
OBJ        = ${SRC:.cpp=.o}
CLEANFILES = ${LIBS} ${addsuffix /*.o, ${DIRS}} ${wildcard bin/*}

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

PETSC_CC_INCLUDES += -I./src/include -I./src/solvers
PCC_FLAGS += -std=c++0x -Wextra -pedantic

bin/PetIBM: ${OBJ} ${LIBS}
	${CLINKER} $^ -o $@ ${PETSC_SYS_LIB}

lib/libclasses.a:
	cd src/include; ${MAKE}

lib/libyaml.a:
	cd external/yaml-cpp; ${MAKE}

lib/libsolvers.a:
	cd src/solvers; ${MAKE}

check2d:
	${MPIEXEC} -n 4 bin/PetIBM -caseFolder cases/2d/test

memcheck2d:
	valgrind --tool=memcheck --leak-check=full --show-reachable=yes bin/PetIBM -caseFolder cases/2d/memtest

cavity:
	${MPIEXEC} -n 4 bin/PetIBM -caseFolder cases/2d/cavityRe100 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

body2d:
	${MPIEXEC} -n 4 bin/PetIBM -caseFolder cases/2d/bodyTest -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

cylinder:
	${MPIEXEC} -n 4 bin/PetIBM -caseFolder cases/2d/cylinderRe40 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

cylinderPeriodicDomain:
	${MPIEXEC} -n 4 bin/PetIBM -caseFolder cases/2d/cylinderPeriodicDomain -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

memcheck3d:
	valgrind --tool=memcheck --leak-check=full --show-reachable=yes bin/PetIBM -caseFolder cases/3d/memtest

cavityX:
	${MPIEXEC} -n 4 bin/PetIBM -caseFolder cases/3d/cavityX -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

cavityY:
	${MPIEXEC} -n 4 bin/PetIBM -caseFolder cases/3d/cavityY -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

cavityZ:
	${MPIEXEC} -n 4 bin/PetIBM -caseFolder cases/3d/cavityZ -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

body3d:
	${MPIEXEC} -n 4 bin/PetIBM -caseFolder cases/3d/bodyTest -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

vars:
	@echo CLINKER: ${CLINKER}
	@echo CXX: ${CXX}
	@echo PCC_LINKER: ${PCC_LINKER}
	@echo PCC_LINKER_FLAGS: ${PCC_LINKER_FLAGS}
	@echo CFLAGS: ${CFLAGS}
	@echo CC_FLAGS: ${CC_FLAGS}
	@echo PCC_FLAGS: ${PCC_FLAGS}
	@echo CPPFLAGS: ${CPPFLAGS}
	@echo CPP_FLAGS: ${CPP_FLAGS}
	@echo CCPPFLAGS: ${CCPPFLAGS}
	@echo PETSC_CC_INCLUDES: ${PETSC_CC_INCLUDES}
	@echo RM: ${RM}
	@echo MV: ${MV}
	@echo MAKE: ${MAKE}
	@echo MFLAGS: ${MFLAGS}
	@echo OMAKE: ${OMAKE}
	@echo AR: ${AR}
	@echo ARFLAGS: ${ARFLAGS}
	@echo RANLIB: ${RANLIB}
	@echo SRC: ${SRC}
	@echo OBJ: ${OBJ}
	@echo CLEANFILES: ${CLEANFILES}
	@echo FIND: ${FIND}

cleanoutput:
	find . -name '*.d' -exec rm -rf {} \;
	find ./cases -name '*.txt' -exec rm -rf {} \;
	find ./cases -name '0*' -prune -exec rm -rf {} \;
	find ./cases -name 'output' -prune -exec rm -rf {} \;

cleanall: cleanout
	${RM} ${YAMLOBJ}

.PHONY: ${LIBS} check4 memcheck vars cleanall
