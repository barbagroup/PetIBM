# file: makefile
# author: Olivier Mesnard (mesnardo@gwu.edu)
# brief: Builds the unit-test for the explicit diffusion term.


DIFFUSION_DIR = $(abspath $(dir $(firstword $(MAKEFILE_LIST))))

DIFFUSION2D = $(DIFFUSION_DIR)/diffusion2d
DIFFUSION3D = $(DIFFUSION_DIR)/diffusion3d

LIB_DIR = $(PETIBM_DIR)/lib
LIBS = $(addprefix $(LIB_DIR)/, libclasses.a libsolvers.a)
EXT_LIBS = $(addprefix $(LIB_DIR)/, libyaml.a)

.PHONY: ALL variables cleanDiffusion

ALL: $(DIFFUSION2D) $(DIFFUSION3D)

include $(PETSC_DIR)/conf/variables
include $(PETSC_DIR)/conf/rules

PETSC_CC_INCLUDES += -I $(PETIBM_DIR)/src/include -I $(PETIBM_DIR)/src/solvers

PCC_FLAGS += -std=c++0x -Wextra -pedantic
CXX_FLAGS += -std=c++0x -Wextra -pedantic

$(DIFFUSION2D): $(DIFFUSION_DIR)/main2d.o $(DIFFUSION_DIR)/DiffusionTerm2d.o $(LIBS) $(EXT_LIBS)
	@echo "\n$@ - Linking ..."
	$(CLINKER) $^ -o $@ $(PETSC_SYS_LIB)

$(DIFFUSION3D): $(DIFFUSION_DIR)/main3d.o $(DIFFUSION_DIR)/DiffusionTerm3d.o $(LIBS) $(EXT_LIBS)
	@echo "\n$@ - Linking ..."
	$(CLINKER) $^ -o $@ $(PETSC_SYS_LIB)

$(DIFFUSION_DIR)/main2d.o: $(DIFFUSION_DIR)/main.cpp
	@echo "\n$@ - Compiling ..."
	$(PETSC_COMPILE) -D DIMENSIONS=2 $^ -o $@

$(DIFFUSION_DIR)/main3d.o: $(DIFFUSION_DIR)/main.cpp
	@echo "\n$@ - Compiling ..."
	$(PETSC_COMPILE) -D DIMENSIONS=3 $^ -o $@

$(DIFFUSION_DIR)/DiffusionTerm2d.o: $(DIFFUSION_DIR)/DiffusionTerm.cpp
	@echo "\n$@ - Compiling ..."
	$(PETSC_COMPILE) -D DIMENSIONS=2 $^ -o $@

$(DIFFUSION_DIR)/DiffusionTerm3d.o: $(DIFFUSION_DIR)/DiffusionTerm.cpp
	@echo "\n$@ - Compiling ..."
	$(PETSC_COMPILE) -D DIMENSIONS=3 $^ -o $@

cleanDiffusion:
	$(RM) -f $(DIFFUSION2D) $(DIFFUSION3D) 
	$(RM) -f $(DIFFUSION_DIR)/*.o
	$(RM) -f $(DIFFUSION_DIR)/*.dat

runDiffusion2d:
	$(MPIEXEC) -n 1 $(DIFFUSION2D) -caseFolder $(PETIBM_DIR)/cases/2d/tests/6 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	$(MPIEXEC) -n 1 $(DIFFUSION2D) -caseFolder $(PETIBM_DIR)/cases/2d/tests/12 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	$(MPIEXEC) -n 1 $(DIFFUSION2D) -caseFolder $(PETIBM_DIR)/cases/2d/tests/24 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	$(MPIEXEC) -n 1 $(DIFFUSION2D) -caseFolder $(PETIBM_DIR)/cases/2d/tests/48 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

runDiffusion3d:
	$(MPIEXEC) -n 1 $(DIFFUSION3D) -caseFolder $(PETIBM_DIR)/cases/3d/tests/6 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	$(MPIEXEC) -n 1 $(DIFFUSION3D) -caseFolder $(PETIBM_DIR)/cases/3d/tests/12 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	$(MPIEXEC) -n 1 $(DIFFUSION3D) -caseFolder $(PETIBM_DIR)/cases/3d/tests/24 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	$(MPIEXEC) -n 1 $(DIFFUSION3D) -caseFolder $(PETIBM_DIR)/cases/3d/tests/48 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	
memoryCheck2dSerial:
	$(MPIEXEC) -n 1 valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes $(DIFFUSION2D) -caseFolder $(PETIBM_DIR)/cases/2d/tests/6

memoryCheck2dParallel:
	$(MPIEXEC) -n 2 valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes $(DIFFUSION2D) -caseFolder $(PETIBM_DIR)/cases/2d/tests/6

memoryCheck3dSerial:
	$(MPIEXEC) -n 1 valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes $(DIFFUSION3D) -caseFolder $(PETIBM_DIR)/cases/3d/tests/6

memoryCheck3dParallel:
	$(MPIEXEC) -n 2 valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes $(DIFFUSION3D) -caseFolder $(PETIBM_DIR)/cases/3d/tests/6

variables:
	@echo PETSC_DIR: $(PETSC_DIR)
	@echo PETIBM_DIR: $(PETIBM_DIR)
	@echo DIFFUSION_DIR: $(DIFFUSION_DIR)
	@echo DIFFUSION2D: $(DIFFUSION2D)
	@echo DIFFUSION3D: $(DIFFUSION3D)
	@echo LIB_DIR: $(LIB_DIR)
	@echo LIBS: $(LIBS)
	@echo EXT_LIBS: $(EXT_LIBS)