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

memcheck:
	valgrind --tool=memcheck --leak-check=full --show-reachable=yes bin/PetIBM -caseFolder cases/test

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

cleanall: clean
	${RM} ${YAMLOBJ}

.PHONY: ${LIBS} check4 memcheck vars cleanall
