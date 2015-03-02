# file: makefile
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: Generates static library libsolvers.a.


TARGET = libsolvers.a

SUFFIX = .cpp
SRCS = $(wildcard *$(SUFFIX))
OBJ = $(SRCS:$(SUFFIX)=.o)
LIB_DIR = $(PETIBM_DIR)/lib
CLEANFILES = $(OBJ)

ALL: $(TARGET)

include $(PETSC_DIR)/conf/variables
include $(PETSC_DIR)/conf/rules

PETSC_CC_INCLUDES += -I../../external/yaml-cpp/include -I../include
PCC_FLAGS += -std=c++0x -Wall -Wextra -pedantic -MMD
CXX_FLAGS += -std=c++0x -Wall -Wextra -pedantic -MMD

$(TARGET): $(OBJ)
	$(AR) $(ARFLAGS) $@ $^
	$(RANLIB) $@
	$(MV) $@ $(LIB_DIR)

.PHONY: ALL $(TARGET)