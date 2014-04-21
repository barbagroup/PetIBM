ALL: libyaml.a
DIRS       = src src/contrib
SRC        = ${wildcard src/*.cpp} ${wildcard src/contrib/*.cpp}
OBJ        = ${SRC:.cpp=.o}
CLEANFILES = ${OBJ}

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

PETSC_CC_INCLUDES += -I./include

libyaml.a: ${OBJ}
	${AR} ${ARFLAGS} $@ $^
	${RANLIB} $@
	${MV} $@ ../../lib/
