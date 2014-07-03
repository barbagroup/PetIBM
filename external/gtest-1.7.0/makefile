ALL: libgtest.a
GTEST_HEADERS = include/gtest/*.h include/gtest/internal/*.h
GTEST_SRCS = src/*.cc src/*.h $(GTEST_HEADERS)
GTEST_OBJS = src/gtest-all.o
CLEANFILES = ${GTEST_OBJS}

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

PETSC_CC_INCLUDES += -I ./include -I ./

libgtest.a: ${GTEST_OBJS}
	${AR} ${ARFLAGS} $@ $^
	${RANLIB} $@
	${MV} $@ ../../lib/
