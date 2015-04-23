# file: makefile
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: Generates static library libclasses.a.


TARGET = libclasses.a

SUFFIX = .cpp
SRCS = $(wildcard *$(SUFFIX))
OBJ = $(SRCS:$(SUFFIX)=.o)
CLEANFILES = $(OBJ)

ALL: $(TARGET)

include $(PETSC_DIR)/conf/variables
include $(PETSC_DIR)/conf/rules

PETSC_CC_INCLUDES += -I $(YAML)/include -I $(BOOST_INCLUDE)

$(TARGET): $(OBJ)
	$(AR) $(ARFLAGS) $@ $^
	$(RANLIB) $@
	$(MV) $@ $(LIB_DIR)

.PHONY: ALL $(TARGET)