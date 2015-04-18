# file: makefile
# author: Olivier Mesnard (mesnardo@gwu.edu)
# brief: Builds the unit-test for the explicit convection term.


CONVECTION_DIR = $(abspath $(dir $(firstword $(MAKEFILE_LIST))))

CONVECTION2D = $(CONVECTION_DIR)/convection2d
CONVECTION3D = $(CONVECTION_DIR)/convection3d

LIB_DIR = $(PETIBM_DIR)/lib
LIBS = $(addprefix $(LIB_DIR)/, libclasses.a libsolvers.a)
EXT_LIBS = $(addprefix $(LIB_DIR)/, libyaml.a)

.PHONY: ALL variables cleanConvection

ALL: $(CONVECTION2D) $(CONVECTION3D)

include $(PETSC_DIR)/conf/variables
include $(PETSC_DIR)/conf/rules

PETSC_CC_INCLUDES += -I $(PETIBM_DIR)/src/include -I $(PETIBM_DIR)/src/solvers

PCC_FLAGS += -std=c++0x -Wextra -pedantic
CXX_FLAGS += -std=c++0x -Wextra -pedantic

$(CONVECTION2D): $(CONVECTION_DIR)/main2d.o $(CONVECTION_DIR)/ConvectionTerm2d.o $(LIBS) $(EXT_LIBS)
	@echo "\n$@ - Linking ..."
	$(CLINKER) $^ -o $@ $(PETSC_SYS_LIB)

$(CONVECTION3D): $(CONVECTION_DIR)/main3d.o $(CONVECTION_DIR)/ConvectionTerm3d.o $(LIBS) $(EXT_LIBS)
	@echo "\n$@ - Linking ..."
	$(CLINKER) $^ -o $@ $(PETSC_SYS_LIB)

$(CONVECTION_DIR)/main2d.o: $(CONVECTION_DIR)/main.cpp
	@echo "\n$@ - Compiling ..."
	$(PETSC_COMPILE) -D DIMENSIONS=2 $^ -o $@

$(CONVECTION_DIR)/main3d.o: $(CONVECTION_DIR)/main.cpp
	@echo "\n$@ - Compiling ..."
	$(PETSC_COMPILE) -D DIMENSIONS=3 $^ -o $@

$(CONVECTION_DIR)/ConvectionTerm2d.o: $(CONVECTION_DIR)/ConvectionTerm.cpp
	@echo "\n$@ - Compiling ..."
	$(PETSC_COMPILE) -D DIMENSIONS=2 $^ -o $@

$(CONVECTION_DIR)/ConvectionTerm3d.o: $(CONVECTION_DIR)/ConvectionTerm.cpp
	@echo "\n$@ - Compiling ..."
	$(PETSC_COMPILE) -D DIMENSIONS=3 $^ -o $@

cleanConvection:
	$(RM) -f $(CONVECTION2D) $(CONVECTION3D) 
	$(RM) -f $(CONVECTION_DIR)/*.o
	$(RM) -f $(CONVECTION_DIR)/*.dat

runConvection2d:
	$(MPIEXEC) -n 1 $(CONVECTION2D) -caseFolder $(PETIBM_DIR)/cases/2d/tests/6 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	$(MPIEXEC) -n 1 $(CONVECTION2D) -caseFolder $(PETIBM_DIR)/cases/2d/tests/12 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	$(MPIEXEC) -n 1 $(CONVECTION2D) -caseFolder $(PETIBM_DIR)/cases/2d/tests/24 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	$(MPIEXEC) -n 1 $(CONVECTION2D) -caseFolder $(PETIBM_DIR)/cases/2d/tests/48 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

runConvection3d:
	$(MPIEXEC) -n 1 $(CONVECTION3D) -caseFolder $(PETIBM_DIR)/cases/3d/tests/6 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	$(MPIEXEC) -n 1 $(CONVECTION3D) -caseFolder $(PETIBM_DIR)/cases/3d/tests/12 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	$(MPIEXEC) -n 1 $(CONVECTION3D) -caseFolder $(PETIBM_DIR)/cases/3d/tests/24 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	$(MPIEXEC) -n 1 $(CONVECTION3D) -caseFolder $(PETIBM_DIR)/cases/3d/tests/48 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	
memoryCheck2dSerial:
	$(MPIEXEC) -n 1 valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes $(CONVECTION2D) -caseFolder $(PETIBM_DIR)/cases/2d/tests/6

memoryCheck2dParallel:
	$(MPIEXEC) -n 2 valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes $(CONVECTION2D) -caseFolder $(PETIBM_DIR)/cases/2d/tests/6

memoryCheck3dSerial:
	$(MPIEXEC) -n 1 valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes $(CONVECTION3D) -caseFolder $(PETIBM_DIR)/cases/3d/tests/6

memoryCheck3dParallel:
	$(MPIEXEC) -n 2 valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes $(CONVECTION3D) -caseFolder $(PETIBM_DIR)/cases/3d/tests/6

variables:
	@echo PETSC_DIR: $(PETSC_DIR)
	@echo PETIBM_DIR: $(PETIBM_DIR)
	@echo DIFFUSION_DIR: $(CONVECTION_DIR)
	@echo DIFFUSION2D: $(CONVECTION2D)
	@echo DIFFUSION3D: $(CONVECTION3D)
	@echo LIB_DIR: $(LIB_DIR)
	@echo LIBS: $(LIBS)
	@echo EXT_LIBS: $(EXT_LIBS)