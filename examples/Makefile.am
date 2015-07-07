include $(PETSC_DIR)/conf/variables

PETIBM2D = $(top_builddir)/src/petibm2d
PETIBM3D = $(top_builddir)/src/petibm3d

examples:
	@echo "Copying examples folders..."
	@if test "$(top_srcdir)" != "$(top_builddir)" ; then \
		cp -r $(top_srcdir)/examples/2d $(PWD) ; \
		cp -r $(top_srcdir)/examples/3d $(PWD) ; \
	fi ;


# Two-dimensional examples

cavity2dRe100Serial:
	@echo "2D lid-driven cavity flow at Re=100 on 1 process."
	@$(MPIEXEC) -n 1 $(PETIBM2D) \
							-caseFolder ./2d/lidDrivenCavity/Re100 \
							-sys2_pc_type gamg \
							-sys2_pc_gamg_type agg \
							-sys2_pc_gamg_agg_nsmooths 1

cavity2dRe100Parallel:
	@echo "2D lid-driven cavity flow at Re=100 on 4 processes."
	@$(MPIEXEC) -n 4 $(PETIBM2D) \
							-caseFolder ./2d/lidDrivenCavity/Re100 \
							-sys2_pc_type gamg \
							-sys2_pc_gamg_type agg \
							-sys2_pc_gamg_agg_nsmooths 1

cavity2dRe100NonUniform:
	@echo "2D lid-driven cavity flow at Re=100 on 1 process."
	@$(MPIEXEC) -n 1 $(PETIBM2D) \
							-caseFolder ./2d/lidDrivenCavity/Re100NonUniform \
							-sys2_pc_type gamg \
							-sys2_pc_gamg_type agg \
							-sys2_pc_gamg_agg_nsmooths 1

cavity2dRe1000:
	@echo "2D lid-driven cavity flow at Re=1000 on 4 processes."
	@$(MPIEXEC) -n 4 $(PETIBM2D) \
							-caseFolder ./2d/lidDrivenCavity/Re1000 \
							-sys2_pc_type gamg \
							-sys2_pc_gamg_type agg \
							-sys2_pc_gamg_agg_nsmooths 1

cavity2dRe3200:
	@echo "2D lid-driven cavity flow at Re=3200 on 4 processes."
	@$(MPIEXEC) -n 4 $(PETIBM2D) \
							-caseFolder ./2d/lidDrivenCavity/Re3200 \
							-sys2_pc_type gamg \
							-sys2_pc_gamg_type agg \
							-sys2_pc_gamg_agg_nsmooths 1

cavity2dRe5000:
	@echo "2D lid-driven cavity flow at Re=5000 on 4 processes."
	@$(MPIEXEC) -n 4 $(PETIBM2D) \
							-caseFolder ./2d/lidDrivenCavity/Re5000 \
							-sys2_pc_type gamg \
							-sys2_pc_gamg_type agg \
							-sys2_pc_gamg_agg_nsmooths 1

cylinder2dRe40:
	@echo "2D cylinder at Re=40 on 2 processes."
	@$(MPIEXEC) -n 2 $(PETIBM2D) \
							-caseFolder ./2d/cylinder/Re40 \
							-sys2_pc_type gamg \
							-sys2_pc_gamg_type agg \
							-sys2_pc_gamg_agg_nsmooths 1

cylinder2dRe40Periodic:
	@echo "2D cylinder at Re=40 with periodic boundary conditions on 4 processes."
	@$(MPIEXEC) -n 4 $(PETIBM2D) \
							-caseFolder ./2d/cylinder/Re40PeriodicDomain \
							-sys2_pc_type gamg \
							-sys2_pc_gamg_type agg \
							-sys2_pc_gamg_agg_nsmooths 1

cylinder2dRe150:
	@echo "2D cylinder at Re=150 on 4 processes."
	@$(MPIEXEC) -n 4 $(PETIBM2D) \
							-caseFolder ./2d/cylinder/Re150 \
							-sys2_pc_type gamg \
							-sys2_pc_gamg_type agg \
							-sys2_pc_gamg_agg_nsmooths 1

cylinder2dRe250:
	@echo "2D cylinder at Re=250 on 4 processes."
	@$(MPIEXEC) -n 4 $(PETIBM2D) \
							-caseFolder ./2d/cylinder/Re250 \
							-sys2_pc_type gamg \
							-sys2_pc_gamg_type agg \
							-sys2_pc_gamg_agg_nsmooths 1

cylinder2dRe550:
	@echo "2D cylinder at Re=550 on 4 processes."
	@$(MPIEXEC) -n 4 $(PETIBM2D) \
							-caseFolder ./2d/cylinder/Re550 \
							-sys2_pc_type gamg \
							-sys2_pc_gamg_type agg \
							-sys2_pc_gamg_agg_nsmooths 1

cylinder2dRe3000:
	@echo "2D cylinder at Re=3000 on 4 processes."
	@$(MPIEXEC) -n 4 $(PETIBM2D) \
							-caseFolder ./2d/cylinder/Re3000 \
							-sys2_pc_type gamg \
							-sys2_pc_gamg_type agg \
							-sys2_pc_gamg_agg_nsmooths 1

memoryCheck2dSerial:
	@echo "Checking memory leak on a 2d-case on 1 process..."
	@$(MPIEXEC) -n 1 valgrind \
							--tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes \
							$(PETIBM2D) -caseFolder ./2d/memoryTest

memoryCheck2dParallel:
	@echo "Checking memory leak on a 2d-case on 2 process..."
	@$(MPIEXEC) -n 2 valgrind \
							--tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes \
							$(PETIBM2D) -caseFolder ./2d/memoryTest

memoryCheck2dBodySerial:
	@echo "Checking memory leak on a 2d-case with an immersed body on 1 process..."
	@$(MPIEXEC) -n 1 valgrind \
							--tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes \
							$(PETIBM2D) -caseFolder ./2d/memoryTestBody


# Three-dimensional examples

cavity3dRe100PeriodicX:
	@echo "3D lid-driven cavity flow at Re=100 with periodic boundary conditions in the x-direction on 4 processes."
	@$(MPIEXEC) -n 4 $(PETIBM3D) \
							-caseFolder ./3d/lidDrivenCavity/Re100PeriodicX \
							-sys2_pc_type gamg \
							-sys2_pc_gamg_type agg \
							-sys2_pc_gamg_agg_nsmooths 1

cavity3dRe100PeriodicY:
	@echo "3D lid-driven cavity flow at Re=100 with periodic boundary conditions in the y-direction on 4 processes."
	@$(MPIEXEC) -n 4 $(PETIBM3D) \
							-caseFolder ./3d/lidDrivenCavity/Re100PeriodicY \
							-sys2_pc_type gamg \
							-sys2_pc_gamg_type agg \
							-sys2_pc_gamg_agg_nsmooths 1

cavity3dRe100PeriodicZ:
	@echo "3D lid-driven cavity flow at Re=100 with periodic boundary conditions in the z-direction on 4 processes."
	@$(MPIEXEC) -n 4 $(PETIBM3D) \
							-caseFolder ./3d/lidDrivenCavity/Re100PeriodicZ \
							-sys2_pc_type gamg \
							-sys2_pc_gamg_type agg \
							-sys2_pc_gamg_agg_nsmooths 1

cylinder3dRe40:
	@echo "3D cylinder at Re=40 on 4 processes."
	@$(MPIEXEC) -n 4 $(PETIBM3D) \
							-caseFolder ./3d/cylinder/Re40 \
							-sys2_pc_type gamg \
							-sys2_pc_gamg_type agg \
							-sys2_pc_gamg_agg_nsmooths 1

memoryCheck3dSerial:
	@echo "Checking memory leak on a 3d-case on 1 process..."
	@$(MPIEXEC) -n 1 valgrind \
							--tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes \
							$(PETIBM3D) -caseFolder ./3d/memoryTest

memoryCheck3dParallel:
	@echo "Checking memory leak on a 3d-case on 2 process..."
	@$(MPIEXEC) -n 2 valgrind \
							--tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes \
							$(PETIBM3D) -caseFolder ./2d/memoryTest							