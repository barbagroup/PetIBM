include $(PETSC_DIR)/lib/petsc/conf/variables

PETIBM2D = $(top_builddir)/src/petibm2d
PETIBM3D = $(top_builddir)/src/petibm3d

examples:
	@echo "Copying examples folders..."
	@if test "$(top_srcdir)" != "$(top_builddir)" ; then \
		cp -r $(top_srcdir)/examples/2d $(PWD) ; \
		cp -r $(top_srcdir)/examples/3d $(PWD) ; \
		cp -r $(top_srcdir)/examples/GPU $(PWD) ; \
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
	@echo "==================================================="
	@echo "2d lid-driven cavity flow at Re=100 using 1 process"
	@echo "==================================================="
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
							--tool=memcheck --leak-check=full --track-origins=yes \
							$(PETIBM2D) -directory ./2d/memoryTest

memoryCheck2dParallel:
	@echo "=============================================="
	@echo "Check memory leak on 2d-case using 2 processes"
	@echo "=============================================="
	@$(MPIEXEC) -n 2 valgrind \
							--tool=memcheck --leak-check=full --track-origins=yes \
							$(PETIBM2D) -directory ./2d/memoryTest

memoryCheck2dBodySerial:
	@echo "==============================================================="
	@echo "Check memory leak on 2d-case with immersed-body using 1 process"
	@echo "==============================================================="
	@$(MPIEXEC) -n 1 valgrind \
							--tool=memcheck --leak-check=full --track-origins=yes \
							$(PETIBM2D) -directory ./2d/memoryTestBody

convergence2d:
	@echo "================================="
	@echo "Two-dimensional convergence tests"
	@echo "================================="
	@$(MAKE) lidDrivenCavity2d20
	@$(MAKE) lidDrivenCavity2d30
	@$(MAKE) decayingVortices2d20
	@$(MAKE) decayingVortices2d30
	@$(MAKE) translatingVortices2d20
	@$(MAKE) translatingVortices2d30

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
	@echo "====================================="
	@echo "2d lid-driven cavity flow (1620x1620)"
	@echo "====================================="
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/lidDrivenCavity20/1620

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
	@echo "====================================="
	@echo "2d lid-driven cavity flow (2430x2430)"
	@echo "====================================="
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/lidDrivenCavity30/2430

decayingVortices2d20:
	@echo "=============================="
	@echo "2d Taylor-Green vortex (20x20)"
	@echo "=============================="
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/decayingVortices20/20
	@echo "=============================="
	@echo "2d Taylor-Green vortex (60x60)"
	@echo "=============================="
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/decayingVortices20/60
	@echo "================================"
	@echo "2d Taylor-Green vortex (180x180)"
	@echo "================================"
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/decayingVortices20/180
	@echo "================================"
	@echo "2d Taylor-Green vortex (540x540)"
	@echo "================================"
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/decayingVortices20/540

decayingVortices2d30:
	@echo "=============================="
	@echo "2d Taylor-Green vortex (30x30)"
	@echo "=============================="
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/decayingVortices30/30
	@echo "=============================="
	@echo "2d Taylor-Green vortex (90x90)"
	@echo "=============================="
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/decayingVortices30/90
	@echo "================================"
	@echo "2d Taylor-Green vortex (270x270)"
	@echo "================================"
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/decayingVortices30/270
	@echo "================================"
	@echo "2d Taylor-Green vortex (810x810)"
	@echo "================================"
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/decayingVortices30/810

translatingVortices2d20:
	@echo "==============================="
	@echo "2d translating vortices (20x20)"
	@echo "==============================="
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/translatingVortices20/20
	@echo "==============================="
	@echo "2d translating vortices (60x60)"
	@echo "==============================="
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/translatingVortices20/60
	@echo "================================="
	@echo "2d translating vortices (180x180)"
	@echo "================================="
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/translatingVortices20/180
	@echo "================================"
	@echo "2d translatingVortices (540x540)"
	@echo "================================"
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/translatingVortices20/540

translatingVortices2d30:
	@echo "==============================="
	@echo "2d translating vortices (30x30)"
	@echo "==============================="
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/translatingVortices30/30
	@echo "==============================="
	@echo "2d translating vortices (90x90)"
	@echo "==============================="
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/translatingVortices30/90
	@echo "================================="
	@echo "2d translating vortices (270x270)"
	@echo "================================="
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/translatingVortices30/270
	@echo "================================"
	@echo "2d translatingVortices (810x810)"
	@echo "================================"
	@$(MPIEXEC) -n 1 $(PETIBM2D) -directory ./2d/convergence/translatingVortices30/810

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

sphere3dRe50:
	@echo "================================================"
	@echo "3d flow over a sphere at Re=50 using 8 processes"
	@echo "================================================"
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/sphere/Re50

sphere3dRe100:
	@echo "================================================="
	@echo "3d flow over a sphere at Re=100 using 8 processes"
	@echo "================================================="
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/sphere/Re100

sphere3dRe150:
	@echo "================================================="
	@echo "3d flow over a sphere at Re=150 using 8 processes"
	@echo "================================================="
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/sphere/Re150

sphere3dRe200:
	@echo "================================================="
	@echo "3d flow over a sphere at Re=200 using 8 processes"
	@echo "================================================="
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/sphere/Re200

flatPlate3dRe100AR2dx04:
	@echo "==============================================================="
	@echo "3d flow over an inclined flat-plate (AR=2) at Re=100 using 8 processes"
	@echo "grid-spacing in uniform region: dx=0.04"
	@echo "==============================================================="
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.04/flatPlateRe100AoA0AR2
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.04/flatPlateRe100AoA10AR2
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.04/flatPlateRe100AoA20AR2
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.04/flatPlateRe100AoA30AR2
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.04/flatPlateRe100AoA40AR2
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.04/flatPlateRe100AoA50AR2
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.04/flatPlateRe100AoA60AR2
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.04/flatPlateRe100AoA70AR2
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.04/flatPlateRe100AoA80AR2
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.04/flatPlateRe100AoA90AR2

flatPlate3dRe100AR2dx03:
	@echo "==============================================================="
	@echo "3d flow over an inclined flat-plate (AR=2) at Re=100 using 8 processes"
	@echo "grid-spacing in uniform region: dx=0.03"
	@echo "==============================================================="
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.03/flatPlateRe100AoA0AR2
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.03/flatPlateRe100AoA10AR2
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.03/flatPlateRe100AoA20AR2
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.03/flatPlateRe100AoA30AR2
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.03/flatPlateRe100AoA40AR2
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.03/flatPlateRe100AoA50AR2
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.03/flatPlateRe100AoA60AR2
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.03/flatPlateRe100AoA70AR2
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.03/flatPlateRe100AoA80AR2
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.03/flatPlateRe100AoA90AR2

flatPlate3dRe100AR2dx025:
	@echo "==============================================================="
	@echo "3d flow over an inclined flat-plate (AR=2) at Re=100 using 8 processes"
	@echo "grid-spacing in uniform region: dx=0.025"
	@echo "==============================================================="
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.025/flatPlateRe100AoA0AR2
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.025/flatPlateRe100AoA10AR2
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.025/flatPlateRe100AoA20AR2
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.025/flatPlateRe100AoA30AR2
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.025/flatPlateRe100AoA40AR2
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.025/flatPlateRe100AoA50AR2
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.025/flatPlateRe100AoA60AR2
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.025/flatPlateRe100AoA70AR2
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.025/flatPlateRe100AoA80AR2
	@$(MPIEXEC) -n 8 $(PETIBM3D) -directory ./3d/flatePlate/resolution0.025/flatPlateRe100AoA90AR2

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
							--tool=memcheck --leak-check=full --track-origins=yes \
							$(PETIBM3D) -directory ./3d/memoryTest

memoryCheck3dParallel:
	@echo "=============================================="
	@echo "Check memory leak on 3d-case using 2 processes"
	@echo "=============================================="
	@$(MPIEXEC) -n 2 valgrind \
							--tool=memcheck --leak-check=full --track-origins=yes \
							$(PETIBM3D) -directory ./3d/memoryTest							

# GPU examples

cylinder2dRe40GPU:
	@echo "================================================================"
	@echo "2d flow over cylinder at Re=40 using 2 processes and AmgX solver"
	@echo "================================================================"
	@$(MPIEXEC) -n 2 $(PETIBM2D) -directory ./GPU/2d/cylinder/Re40

cylinder2dRe550GPU:
	@echo "================================================================="
	@echo "2d flow over cylinder at Re=550 using 4 processes and AmgX solver"
	@echo "================================================================="
	@$(MPIEXEC) -n 4 $(PETIBM2D) -directory ./GPU/2d/cylinder/Re550

cylinder2dRe3000GPU:
	@echo "=================================================================="
	@echo "2d flow over cylinder at Re=3000 using 4 processes and AmgX solver"
	@echo "=================================================================="
	@$(MPIEXEC) -n 4 $(PETIBM2D) -directory ./GPU/2d/cylinder/Re3000
