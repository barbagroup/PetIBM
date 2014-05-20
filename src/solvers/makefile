ALL: libsolvers.a
SRC        = ${wildcard *.cpp}
OBJ        = ${SRC:.cpp=.o}
CLEANFILES = ${OBJ}

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

PETSC_CC_INCLUDES += -I../../external/yaml-cpp/include -I../include
PCC_FLAGS += -std=c++0x -Wextra -pedantic

libsolvers.a: ${OBJ}
	${AR} ${ARFLAGS} $@ $^
	${RANLIB} $@
	${MV} $@ ../../lib/
