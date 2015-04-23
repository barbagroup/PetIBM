# file: run_cases.sh
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Runs the 2d Taylor-Green vortex on 4 grids 
#              with constant grid-refinement ratio 
#              and plot the spatial-convergence

PETIBM2D="$PETIBM_DIR/bin/PetIBM2d"
MPIEXEC="$PETSC_DIR/arch-linux2-c-opt/bin/mpiexec"
N=1

# grid: 25x25
$MPIEXEC -n $N $PETIBM2D -caseFolder ./25 \
           -sys2_pc_type gamg \
           -sys2_pc_gamg_type agg \
           -sys2_pc_gamg_agg_nsmooths 1

# grid: 75x75
$MPIEXEC -n $N $PETIBM2D -caseFolder ./75 \
           -sys2_pc_type gamg \
           -sys2_pc_gamg_type agg \
           -sys2_pc_gamg_agg_nsmooths 1

# grid: 225x225
$MPIEXEC -n $N $PETIBM2D -caseFolder ./225 \
           -sys2_pc_type gamg \
           -sys2_pc_gamg_type agg \
           -sys2_pc_gamg_agg_nsmooths 1

# grid: 675x675
$MPIEXEC -n $N $PETIBM2D -caseFolder ./675 \
           -sys2_pc_type gamg \
           -sys2_pc_gamg_type agg \
           -sys2_pc_gamg_agg_nsmooths 1

# plot the spatial convergence
python $PETIBM_DIR/scripts/python/validation/taylorGreenVortexConvergence.py --no-show
