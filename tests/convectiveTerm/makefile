# file: makefile
# author: Olivier Mesnard (mesnardo@gwu.edu)
# brief: Builds the unit-test for the explicit convective term.


CONVECTIVE_DIR = $(abspath $(dir $(firstword $(MAKEFILE_LIST))))
BIN_DIR = $(CONVECTIVE_DIR)/bin
BUILD_DIR = $(CONVECTIVE_DIR)/build
DATA_DIR = $(CONVECTIVE_DIR)/data

CONVECTIVE2D = $(BIN_DIR)/convectiveTerm2d
CONVECTIVE3D = $(BIN_DIR)/convectiveTerm3d

LIB_DIR = $(PETIBM_DIR)/lib
LIBS = $(addprefix $(LIB_DIR)/, libclasses.a libsolvers.a)
EXT_LIBS = $(addprefix $(LIB_DIR)/, libyaml.a)

.PHONY: ALL variables runTest cleanTest

ALL: $(CONVECTIVE2D) $(CONVECTIVE3D)

include $(PETSC_DIR)/conf/variables
include $(PETSC_DIR)/conf/rules

PETSC_CC_INCLUDES += -I $(PETIBM_DIR)/src/include -I $(PETIBM_DIR)/src/solvers

PCC_FLAGS += -std=c++0x -Wextra -pedantic
CXX_FLAGS += -std=c++0x -Wextra -pedantic


$(CONVECTIVE2D): $(BUILD_DIR)/convectiveTermTest2d.o $(BUILD_DIR)/ConvectiveTerm2d.o $(LIBS) $(EXT_LIBS)
	@echo "\n$@ - Linking ..."
	@mkdir -p $(BIN_DIR)
	$(CLINKER) $^ -o $@ $(PETSC_SYS_LIB)

$(CONVECTIVE3D): $(BUILD_DIR)/convectiveTermTest3d.o $(BUILD_DIR)/ConvectiveTerm3d.o $(LIBS) $(EXT_LIBS)
	@echo "\n$@ - Linking ..."
	@mkdir -p $(BIN_DIR)
	$(CLINKER) $^ -o $@ $(PETSC_SYS_LIB)

$(BUILD_DIR)/convectiveTermTest2d.o: $(CONVECTIVE_DIR)/convectiveTermTest.cpp
	@echo "\n$@ - Compiling ..."
	@mkdir -p $(BUILD_DIR)
	$(PETSC_COMPILE) -D DIMENSIONS=2 $^ -o $@

$(BUILD_DIR)/convectiveTermTest3d.o: $(CONVECTIVE_DIR)/convectiveTermTest.cpp
	@echo "\n$@ - Compiling ..."
	@mkdir -p $(BUILD_DIR)
	$(PETSC_COMPILE) -D DIMENSIONS=3 $^ -o $@

$(BUILD_DIR)/ConvectiveTerm2d.o: $(CONVECTIVE_DIR)/ConvectiveTerm.cpp
	@echo "\n$@ - Compiling ..."
	@mkdir -p $(BUILD_DIR)
	$(PETSC_COMPILE) -D DIMENSIONS=2 $^ -o $@

$(BUILD_DIR)/ConvectiveTerm3d.o: $(CONVECTIVE_DIR)/ConvectiveTerm.cpp
	@echo "\n$@ - Compiling ..."
	@mkdir -p $(BUILD_DIR)
	$(PETSC_COMPILE) -D DIMENSIONS=3 $^ -o $@


cleanTest:
	$(RM) -rf $(BIN_DIR)
	$(RM) -rf $(BUILD_DIR)
	$(RM) -rf $(DATA_DIR)


runTest: runConvective2d runConvective3d

runConvective2d:
	@mkdir -p $(DATA_DIR)
	$(RM) -f $(DATA_DIR)/relativeErrors2d.dat
	$(MPIEXEC) -n 1 $(CONVECTIVE2D) -caseFolder $(PETIBM_DIR)/cases/2d/tests/6 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	$(MPIEXEC) -n 1 $(CONVECTIVE2D) -caseFolder $(PETIBM_DIR)/cases/2d/tests/12 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	$(MPIEXEC) -n 1 $(CONVECTIVE2D) -caseFolder $(PETIBM_DIR)/cases/2d/tests/24 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	$(MPIEXEC) -n 1 $(CONVECTIVE2D) -caseFolder $(PETIBM_DIR)/cases/2d/tests/48 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

runConvective3d:
	@mkdir -p $(DATA_DIR)
	$(RM) -f $(DATA_DIR)/relativeErrors3d.dat
	$(MPIEXEC) -n 1 $(CONVECTIVE3D) -caseFolder $(PETIBM_DIR)/cases/3d/tests/6 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	$(MPIEXEC) -n 1 $(CONVECTIVE3D) -caseFolder $(PETIBM_DIR)/cases/3d/tests/12 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	$(MPIEXEC) -n 1 $(CONVECTIVE3D) -caseFolder $(PETIBM_DIR)/cases/3d/tests/24 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	$(MPIEXEC) -n 1 $(CONVECTIVE3D) -caseFolder $(PETIBM_DIR)/cases/3d/tests/48 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	
memoryCheck2dSerial:
	$(MPIEXEC) -n 1 valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes $(CONVECTIVE2D) -caseFolder $(PETIBM_DIR)/cases/2d/tests/6

memoryCheck2dParallel:
	$(MPIEXEC) -n 2 valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes $(CONVECTIVE2D) -caseFolder $(PETIBM_DIR)/cases/2d/tests/6

memoryCheck3dSerial:
	$(MPIEXEC) -n 1 valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes $(CONVECTIVE3D) -caseFolder $(PETIBM_DIR)/cases/3d/tests/6

memoryCheck3dParallel:
	$(MPIEXEC) -n 2 valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes $(CONVECTIVE3D) -caseFolder $(PETIBM_DIR)/cases/3d/tests/6

variables:
	@echo PETSC_DIR: $(PETSC_DIR)
	@echo PETIBM_DIR: $(PETIBM_DIR)
	@echo CONVECTIVE_DIR: $(CONVECTIVE_DIR)
	@echo CONVECTIVE2D: $(CONVECTIVE2D)
	@echo CONVECTIVE3D: $(CONVECTIVE3D)
	@echo BIN_DIR: $(BIN_DIR)
	@echo BUILD_DIR: $(BUILD_DIR)
	@echo DATA_DIR: $(DATA_DIR)
	@echo LIB_DIR: $(LIB_DIR)
	@echo LIBS: $(LIBS)
	@echo EXT_LIBS: $(EXT_LIBS)