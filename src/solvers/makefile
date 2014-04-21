ALL: libsolvers.a
SRC        = ${wildcard *.cpp}
OBJ        = ${SRC:.cpp=.o}
CLEANFILES = ${OBJ}

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

PETSC_CC_INCLUDES += -I../../external/yaml-cpp/include -I../include

libsolvers.a: ${OBJ}
	${AR} ${ARFLAGS} $@ $^
	${RANLIB} $@
	${MV} $@ ../../lib/
