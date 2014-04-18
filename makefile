ALL: bin/PetIBM
DIRS       = src src/include src/solvers external/yaml-cpp/src external/yaml-cpp/src/contrib
LIBS       = lib/libclasses.a lib/libyaml.a lib/libsolvers.a
SRC        = ${wildcard src/*.cpp}
OBJ        = ${SRC:.cpp=.o}
CLEANFILES = ${LIBS} ${addsuffix /*.o, ${DIRS}} ${wildcard bin/*}

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

PETSC_CC_INCLUDES += -I./src/include -I./src/solvers

bin/PetIBM: ${OBJ} ${LIBS}
	${CLINKER} $^ -o $@ ${PETSC_SYS_LIB}

lib/libclasses.a:
	cd src/include; ${MAKE}

lib/libyaml.a:
	cd external/yaml-cpp; ${MAKE}

lib/libsolvers.a:
	cd src/solvers; ${MAKE}

check4:
	${MPIEXEC} -n 4 bin/PetIBM -caseFolder cases/cavityRe100

vars:
	@echo CLINKER: ${CLINKER}
	@echo CXX: ${CXX}
	@echo CFLAGS: ${CFLAGS}
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

.PHONY: ${LIBS} check4 vars
