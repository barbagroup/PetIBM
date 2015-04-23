# file: makefile
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: Generates static library for gtest-1.7.0.


TARGET = libgtest.a

SUFFIX = .cc
SRC_DIR = src
SRCS = $(wildcard $(SRC_DIR)/*$(SUFFIX))
OBJ = $(SRC_DIR)/gtest-all.o
CLEANFILES = $(OBJ)

ALL: $(TARGET)

include $(PETSC_DIR)/conf/variables
include $(PETSC_DIR)/conf/rules

PETSC_CC_INCLUDES += -I ./include -I ./

$(TARGET): $(OBJ)
	$(AR) $(ARFLAGS) $@ $^
	$(RANLIB) $@
	$(MV) $@ $(LIB_DIR)

.PHONY: ALL $(TARGET)