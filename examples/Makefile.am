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

lidDrivenCavity2dRe100Serial:
	@echo "==================================================="
	@echo "2d lid-driven cavity flow at Re=100 using 1 process"
	@echo "==================================================="
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/lidDrivenCavity/Re100

lidDrivenCavity2dRe100Parallel:
	@echo "====================================================="
	@echo "2d lid-driven cavity flow at Re=100 using 4 processes"
	@echo "====================================================="
	@$(MPIEXEC) -n 4 $(PETIBM2D) -directory ./2d/lidDrivenCavity/Re100

lidDrivenCavity2dRe100NonUniform:
	@echo "===================================================="
	@echo "2d lid-driven cavity flow at Re=100 unsing 1 process"
	@echo "===================================================="
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/lidDrivenCavity/Re100NonUniform

lidDrivenCavity2dRe1000:
	@echo "======================================================"
	@echo "2d lid-driven cavity flow at Re=1000 using 4 processes"
	@echo "======================================================"
	@$(MPIEXEC) -n 4 $(PETIBM2D) -directory ./2d/lidDrivenCavity/Re1000

lidDrivenCavity2dRe3200:
	@echo "======================================================"
	@echo "2d lid-driven cavity flow at Re=3200 using 4 processes"
	@echo "======================================================"
	@$(MPIEXEC) -n 4 $(PETIBM2D) -directory ./2d/lidDrivenCavity/Re3200

lidDrivenCavity2dRe5000:
	@echo "======================================================"
	@echo "2d lid-driven cavity flow at Re=5000 using 4 processes"
	@echo "======================================================"
	@$(MPIEXEC) -n 4 $(PETIBM2D) -directory ./2d/lidDrivenCavity/Re5000

cylinder2dRe40:
	@echo "================================================"
	@echo "2d flow over cylinder at Re=40 using 2 processes"
	@echo "================================================"
	@$(MPIEXEC) -n 2 $(PETIBM2D) -directory ./2d/cylinder/Re40

cylinder2dRe40Periodic:
	@echo "=============================================================="
	@echo "2d flow over cylinder at Re=40 (y- periodic) using 4 processes"
	@echo "=============================================================="
	@$(MPIEXEC) -n 4 $(PETIBM2D) -directory ./2d/cylinder/Re40PeriodicDomain

cylinder2dRe150:
	@echo "================================================="
	@echo "2d flow over cylinder at Re=150 using 4 processes"
	@echo "================================================="
	@$(MPIEXEC) -n 4 $(PETIBM2D) -directory ./2d/cylinder/Re150

cylinder2dRe250:
	@echo "================================================="
	@echo "2d flow over cylinder at Re=250 using 4 processes"
	@echo "================================================="
	@$(MPIEXEC) -n 4 $(PETIBM2D) -directory ./2d/cylinder/Re250

cylinder2dRe550:
	@echo "================================================="
	@echo "2d flow over cylinder at Re=550 using 4 processes"
	@echo "================================================="
	@$(MPIEXEC) -n 4 $(PETIBM2D) -directory ./2d/cylinder/Re550

cylinder2dRe3000:
	@echo "=================================================="
	@echo "2d flow over cylinder at Re=3000 using 4 processes"
	@echo "=================================================="
	@$(MPIEXEC) -n 4 $(PETIBM2D) -directory ./2d/cylinder/Re3000

memoryCheck2dSerial:
	@echo "============================================"
	@echo "Check memory leak on 2d-case using 1 process"
	@echo "============================================"
	@$(MPIEXEC) -n 1 valgrind \
							--tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes \
							$(PETIBM2D) -directory ./2d/memoryTest

memoryCheck2dParallel:
	@echo "=============================================="
	@echo "Check memory leak on 2d-case using 2 processes"
	@echo "=============================================="
	@$(MPIEXEC) -n 2 valgrind \
							--tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes \
							$(PETIBM2D) -directory ./2d/memoryTest

memoryCheck2dBodySerial:
	@echo "==============================================================="
	@echo "Check memory leak on 2d-case with immersed-body using 1 process"
	@echo "==============================================================="
	@$(MPIEXEC) -n 1 valgrind \
							--tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes \
							$(PETIBM2D) -directory ./2d/memoryTestBody

convergence2d:
	@echo "================================="
	@echo "Two-dimensional convergence tests"
	@echo "================================="
	@$(MAKE) lidDrivenCavity2d20
	@$(MAKE) lidDrivenCavity2d30
	@$(MAKE) taylorGreenVortex2d20
	@$(MAKE) taylorGreenVortex2d25
	@$(MAKE) taylorGreenVortex2d30

lidDrivenCavity2d20:
	@echo "================================="
	@echo "2d lid-driven cavity flow (20x20)"
	@echo "================================="
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/lidDrivenCavity20/20
	@echo "================================="
	@echo "2d lid-driven cavity flow (60x60)"
	@echo "================================="
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/lidDrivenCavity20/60
	@echo "==================================="
	@echo "2d lid-driven cavity flow (180x180)"
	@echo "==================================="
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/lidDrivenCavity20/180
	@echo "==================================="
	@echo "2d lid-driven cavity flow (540x540)"
	@echo "==================================="
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/lidDrivenCavity20/540

lidDrivenCavity2d30:
	@echo "================================="
	@echo "2d lid-driven cavity flow (30x30)"
	@echo "================================="
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/lidDrivenCavity30/30
	@echo "================================="
	@echo "2d lid-driven cavity flow (90x90)"
	@echo "================================="
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/lidDrivenCavity30/90
	@echo "==================================="
	@echo "2d lid-driven cavity flow (270x270)"
	@echo "==================================="
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/lidDrivenCavity30/270
	@echo "==================================="
	@echo "2d lid-driven cavity flow (810x810)"
	@echo "==================================="
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/lidDrivenCavity30/810

taylorGreenVortex2d20:
	@echo "=============================="
	@echo "2d Taylor-Green vortex (20x20)"
	@echo "=============================="
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/taylorGreenVortex20/20
	@echo "=============================="
	@echo "2d Taylor-Green vortex (60x60)"
	@echo "=============================="
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/taylorGreenVortex20/60
	@echo "================================"
	@echo "2d Taylor-Green vortex (180x180)"
	@echo "================================"
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/taylorGreenVortex20/180
	@echo "================================"
	@echo "2d Taylor-Green vortex (540x540)"
	@echo "================================"
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/taylorGreenVortex20/540

taylorGreenVortex2d25:
	@echo "=============================="
	@echo "2d Taylor-Green vortex (25x25)"
	@echo "=============================="
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/taylorGreenVortex25/25
	@echo "=============================="
	@echo "2d Taylor-Green vortex (75x75)"
	@echo "=============================="
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/taylorGreenVortex25/75
	@echo "================================"
	@echo "2d Taylor-Green vortex (225x225)"
	@echo "================================"
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/taylorGreenVortex25/225
	@echo "================================"
	@echo "2d Taylor-Green vortex (675x675)"
	@echo "================================"
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/taylorGreenVortex25/675

taylorGreenVortex2d30:
	@echo "=============================="
	@echo "2d Taylor-Green vortex (30x30)"
	@echo "=============================="
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/taylorGreenVortex30/30
	@echo "=============================="
	@echo "2d Taylor-Green vortex (90x90)"
	@echo "=============================="
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/taylorGreenVortex30/90
	@echo "================================"
	@echo "2d Taylor-Green vortex (270x270)"
	@echo "================================"
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/taylorGreenVortex30/270
	@echo "================================"
	@echo "2d Taylor-Green vortex (810x810)"
	@echo "================================"
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/taylorGreenVortex30/810


# Three-dimensional examples

lidDrivenCavity3dRe100PeriodicX:
	@echo "=================================================================="
	@echo "3d lid-driven cavity flow at Re=100 (x-periodic) using 4 processes"
	@echo "=================================================================="
	@$(MPIEXEC) -n 4 $(PETIBM3D) -directory ./3d/lidDrivenCavity/Re100PeriodicX

lidDrivenCavity3dRe100PeriodicY:
	@echo "=================================================================="
	@echo "3d lid-driven cavity flow at Re=100 (y-periodic) using 4 processes"
	@echo "=================================================================="
	@$(MPIEXEC) -n 4 $(PETIBM3D) -directory ./3d/lidDrivenCavity/Re100PeriodicY

lidDrivenCavity3dRe100PeriodicZ:
	@echo "=================================================================="
	@echo "3d lid-driven cavity flow at Re=100 (z-periodic) using 4 processes"
	@echo "=================================================================="
	@$(MPIEXEC) -n 4 $(PETIBM3D) -directory ./3d/lidDrivenCavity/Re100PeriodicZ

cylinder3dRe40:
	@echo "================================================"
	@echo "3d flow over cylinder at Re=40 using 4 processes"
	@echo "================================================"
	@$(MPIEXEC) -n 4 $(PETIBM3D) -directory ./3d/cylinder/Re40

memoryCheck3dSerial:
	@echo "============================================"
	@echo "Check memory leak on 3d-case using 1 process"
	@echo "============================================"
	@$(MPIEXEC) -n 1 valgrind \
							--tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes \
							$(PETIBM3D) -directory ./3d/memoryTest

memoryCheck3dParallel:
	@echo "=============================================="
	@echo "Check memory leak on 3d-case using 2 processes"
	@echo "=============================================="
	@$(MPIEXEC) -n 2 valgrind \
							--tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes \
							$(PETIBM3D) -directory ./3d/memoryTest							