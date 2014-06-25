ALL: bin/PetIBM2d bin/PetIBM3d
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

src/PetIBM2d.o:
	${PCC} -o src/PetIBM2d.o -D DIMENSIONS=2 -c ${PCC_FLAGS} ${CFLAGS} ${CCPPFLAGS} src/PetIBM.cpp

src/PetIBM3d.o:
	${PCC} -o src/PetIBM3d.o -D DIMENSIONS=3 -c ${PCC_FLAGS} ${CFLAGS} ${CCPPFLAGS} src/PetIBM.cpp

bin/PetIBM2d: src/PetIBM2d.o ${LIBS}
	${CLINKER} $^ -o $@ ${PETSC_SYS_LIB}

bin/PetIBM3d: src/PetIBM3d.o ${LIBS}
	${CLINKER} $^ -o $@ ${PETSC_SYS_LIB}

lib/libclasses.a:
	cd src/include; ${MAKE}

lib/libyaml.a:
	cd external/yaml-cpp; ${MAKE}

lib/libsolvers.a:
	cd src/solvers; ${MAKE}

check2d:
	${MPIEXEC} -n 4 bin/PetIBM2d -caseFolder cases/2d/test

memcheck2dSerial:
	${MPIEXEC} -n 1 valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes bin/PetIBM2d -caseFolder cases/2d/memtest

memcheck2dParallel:
	${MPIEXEC} -n 2 valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes bin/PetIBM2d -caseFolder cases/2d/memtest

memcheck2dBodySerial:
	${MPIEXEC} -n 1 valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes bin/PetIBM2d -caseFolder cases/2d/memtestBody

cavityRe100Serial:
	${MPIEXEC} -n 1 bin/PetIBM2d -caseFolder cases/2d/lidDrivenCavity/Re100 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

cavityRe100Parallel:
	${MPIEXEC} -n 4 bin/PetIBM2d -caseFolder cases/2d/lidDrivenCavity/Re100 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

cavityRe1000:
	${MPIEXEC} -n 4 bin/PetIBM2d -caseFolder cases/2d/lidDrivenCavity/Re1000 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

body2dSerial:
	${MPIEXEC} -n 1 bin/PetIBM2d -caseFolder cases/2d/bodyTest -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

body2dParallel:
	${MPIEXEC} -n 4 bin/PetIBM2d -caseFolder cases/2d/bodyTest -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

cylinderRe40:
	${MPIEXEC} -n 4 bin/PetIBM2d -caseFolder cases/2d/cylinder/Re40 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

cylinderRe150:
	${MPIEXEC} -n 4 bin/PetIBM2d -caseFolder cases/2d/cylinder/Re150 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

cylinderRe550:
	${MPIEXEC} -n 4 bin/PetIBM2d -caseFolder cases/2d/cylinder/Re550 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

cylinderRe3000:
	${MPIEXEC} -n 4 bin/PetIBM2d -caseFolder cases/2d/cylinder/Re3000 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

cylinderPeriodicDomain:
	${MPIEXEC} -n 4 bin/PetIBM2d -caseFolder cases/2d/cylinderPeriodicDomain -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

memcheck3dSerial:
	${MPIEXEC} -n 1 valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes bin/PetIBM3d -caseFolder cases/3d/memtest

memcheck3dParallel:
	${MPIEXEC} -n 2 valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes bin/PetIBM3d -caseFolder cases/3d/memtest

cavityX:
	${MPIEXEC} -n 4 bin/PetIBM3d -caseFolder cases/3d/cavityX -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

cavityY:
	${MPIEXEC} -n 4 bin/PetIBM3d -caseFolder cases/3d/cavityY -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

cavityZ:
	${MPIEXEC} -n 4 bin/PetIBM3d -caseFolder cases/3d/cavityZ -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

body3dSerial:
	${MPIEXEC} -n 1 bin/PetIBM3d -caseFolder cases/3d/bodyTest -sys2_pc_gamg_agg_nsmooths 1 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

body3dParallel:
	${MPIEXEC} -n 4 bin/PetIBM3d -caseFolder cases/3d/bodyTest -sys2_pc_gamg_agg_nsmooths 1 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

bodyAngle:
	${MPIEXEC} -n 4 bin/PetIBM3d -caseFolder cases/3d/bodyAngle -sys2_pc_gamg_agg_nsmooths 1 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

vortexShedding:
	${MPIEXEC} -n 4 bin/PetIBM3d -caseFolder cases/3d/vortexShedding -sys2_pc_gamg_agg_nsmooths 1 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

flatPlateRe200:
	${MPIEXEC} -n 4 bin/PetIBM3d -caseFolder cases/3d/flatPlateRe200 -sys2_pc_gamg_agg_nsmooths 1 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

cylinder3d:
	${MPIEXEC} -n 4 bin/PetIBM3d -caseFolder cases/3d/cylinder/Re40 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

cylinderRe200:
	${MPIEXEC} -n 4 bin/PetIBM3d -caseFolder cases/3d/cylinder/Re200 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

vars:
	@echo CLINKER: ${CLINKER}
	@echo CXX: ${CXX}
	@echo PCC: ${PCC}
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

cleanall: clean cleanoutput
	${RM} ${YAMLOBJ}

.PHONY: ${LIBS} check4 memcheck vars cleanall
