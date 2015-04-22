# file: makefile
# author: Olivier Mesnard (mesnardo@gwu.edu)
# brief: Builds the unit-test for the explicit diffusive term.


DIFFUSIVE_DIR = $(abspath $(dir $(firstword $(MAKEFILE_LIST))))
BIN_DIR = $(DIFFUSIVE_DIR)/bin
BUILD_DIR = $(DIFFUSIVE_DIR)/build
DATA_DIR = $(DIFFUSIVE_DIR)/data

DIFFUSIVE2D = $(BIN_DIR)/diffusiveTerm2d
DIFFUSIVE3D = $(BIN_DIR)/diffusiveTerm3d

LIB_DIR = $(PETIBM_DIR)/lib
LIBS = $(addprefix $(LIB_DIR)/, libclasses.a libsolvers.a)
EXT_LIBS = $(addprefix $(LIB_DIR)/, libyaml.a)

.PHONY: ALL variables runTest cleanTest

ALL: $(DIFFUSIVE2D) $(DIFFUSIVE3D)

include $(PETSC_DIR)/conf/variables
include $(PETSC_DIR)/conf/rules

PETSC_CC_INCLUDES += -I $(PETIBM_DIR)/src/include -I $(PETIBM_DIR)/src/solvers

PCC_FLAGS += -std=c++0x -Wextra -pedantic
CXX_FLAGS += -std=c++0x -Wextra -pedantic


$(DIFFUSIVE2D): $(BUILD_DIR)/diffusiveTermTest2d.o $(BUILD_DIR)/DiffusiveTerm2d.o $(LIBS) $(EXT_LIBS)
	@echo "\n$@ - Linking ..."
	@mkdir -p $(BIN_DIR)
	$(CLINKER) $^ -o $@ $(PETSC_SYS_LIB)

$(DIFFUSIVE3D): $(BUILD_DIR)/diffusiveTermTest3d.o $(BUILD_DIR)/DiffusiveTerm3d.o $(LIBS) $(EXT_LIBS)
	@echo "\n$@ - Linking ..."
	@mkdir -p $(BIN_DIR)
	$(CLINKER) $^ -o $@ $(PETSC_SYS_LIB)

$(BUILD_DIR)/diffusiveTermTest2d.o: $(DIFFUSIVE_DIR)/diffusiveTermTest.cpp
	@echo "\n$@ - Compiling ..."
	@mkdir -p $(BUILD_DIR)
	$(PETSC_COMPILE) -D DIMENSIONS=2 $^ -o $@

$(BUILD_DIR)/diffusiveTermTest3d.o: $(DIFFUSIVE_DIR)/diffusiveTermTest.cpp
	@echo "\n$@ - Compiling ..."
	@mkdir -p $(BUILD_DIR)
	$(PETSC_COMPILE) -D DIMENSIONS=3 $^ -o $@

$(BUILD_DIR)/DiffusiveTerm2d.o: $(DIFFUSIVE_DIR)/DiffusiveTerm.cpp
	@echo "\n$@ - Compiling ..."
	@mkdir -p $(BUILD_DIR)
	$(PETSC_COMPILE) -D DIMENSIONS=2 $^ -o $@

$(BUILD_DIR)/DiffusiveTerm3d.o: $(DIFFUSIVE_DIR)/DiffusiveTerm.cpp
	@echo "\n$@ - Compiling ..."
	@mkdir -p $(BUILD_DIR)
	$(PETSC_COMPILE) -D DIMENSIONS=3 $^ -o $@


cleanTest:
	$(RM) -rf $(BIN_DIR)
	$(RM) -rf $(BUILD_DIR)
	$(RM) -rf $(DATA_DIR)


runTest: runDiffusive2d runDiffusive3d

runDiffusive2d:
	@mkdir -p $(DATA_DIR)
	$(RM) -f $(DATA_DIR)/relativeErrors2d.dat
	$(MPIEXEC) -n 1 $(DIFFUSIVE2D) -caseFolder $(PETIBM_DIR)/cases/2d/tests/6 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	$(MPIEXEC) -n 1 $(DIFFUSIVE2D) -caseFolder $(PETIBM_DIR)/cases/2d/tests/12 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	$(MPIEXEC) -n 1 $(DIFFUSIVE2D) -caseFolder $(PETIBM_DIR)/cases/2d/tests/24 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	$(MPIEXEC) -n 1 $(DIFFUSIVE2D) -caseFolder $(PETIBM_DIR)/cases/2d/tests/48 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

runDiffusive3d:
	@mkdir -p $(DATA_DIR)
	$(RM) -f $(DATA_DIR)/relativeErrors3d.dat
	$(MPIEXEC) -n 1 $(DIFFUSIVE3D) -caseFolder $(PETIBM_DIR)/cases/3d/tests/6 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	$(MPIEXEC) -n 1 $(DIFFUSIVE3D) -caseFolder $(PETIBM_DIR)/cases/3d/tests/12 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	$(MPIEXEC) -n 1 $(DIFFUSIVE3D) -caseFolder $(PETIBM_DIR)/cases/3d/tests/24 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	$(MPIEXEC) -n 1 $(DIFFUSIVE3D) -caseFolder $(PETIBM_DIR)/cases/3d/tests/48 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	
memoryCheck2dSerial:
	$(MPIEXEC) -n 1 valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes $(DIFFUSIVE2D) -caseFolder $(PETIBM_DIR)/cases/2d/tests/6

memoryCheck2dParallel:
	$(MPIEXEC) -n 2 valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes $(DIFFUSIVE2D) -caseFolder $(PETIBM_DIR)/cases/2d/tests/6

memoryCheck3dSerial:
	$(MPIEXEC) -n 1 valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes $(DIFFUSIVE3D) -caseFolder $(PETIBM_DIR)/cases/3d/tests/6

memoryCheck3dParallel:
	$(MPIEXEC) -n 2 valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes $(DIFFUSIVE3D) -caseFolder $(PETIBM_DIR)/cases/3d/tests/6

variables:
	@echo PETSC_DIR: $(PETSC_DIR)
	@echo PETIBM_DIR: $(PETIBM_DIR)
	@echo DIFFUSIVE_DIR: $(DIFFUSIVE_DIR)
	@echo DIFFUSIVE2D: $(DIFFUSIVE2D)
	@echo DIFFUSIVE3D: $(DIFFUSIVE3D)
	@echo BIN_DIR: $(BIN_DIR)
	@echo BUILD_DIR: $(BUILD_DIR)
	@echo DATA_DIR: $(DATA_DIR)
	@echo LIB_DIR: $(LIB_DIR)
	@echo LIBS: $(LIBS)
	@echo EXT_LIBS: $(EXT_LIBS)