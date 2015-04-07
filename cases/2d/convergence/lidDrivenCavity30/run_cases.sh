# file: run_cases.sh
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Runs the 2d lid-driven cavity flow on 4 grids 
#              with constant grid-refinement ratio 
#              and plot the spatial-convergence

PETIBM2D="$PETIBM_DIR/bin/PetIBM2d"
MPIEXEC="$PETSC_DIR/arch-linux2-c-opt/bin/mpiexec"
N=1

# grid: 30x30
$MPIEXEC -n $N $PETIBM2D -caseFolder ./30 \
           -sys2_pc_type gamg \
           -sys2_pc_gamg_type agg \
           -sys2_pc_gamg_agg_nsmooths 1

# grid: 90x90
$MPIEXEC -n $N $PETIBM2D -caseFolder ./90 \
           -sys2_pc_type gamg \
           -sys2_pc_gamg_type agg \
           -sys2_pc_gamg_agg_nsmooths 1

# grid: 270x270
$MPIEXEC -n $N $PETIBM2D -caseFolder ./270 \
           -sys2_pc_type gamg \
           -sys2_pc_gamg_type agg \
           -sys2_pc_gamg_agg_nsmooths 1

# grid: 810x810
$MPIEXEC -n $N $PETIBM2D -caseFolder ./810 \
           -sys2_pc_type gamg \
           -sys2_pc_gamg_type agg \
           -sys2_pc_gamg_agg_nsmooths 1

# plot the spatial convergence
python $PETIBM_DIR/scripts/python/validation/cavityConvergence.py --no-show
