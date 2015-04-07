# file: run_cases.sh
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Runs the 2d lid-driven cavity flow on 4 grids 
#              with constant grid-refinement ratio 
#              and plot the spatial-convergence

PETIBM2D="$PETIBM_DIR/bin/PetIBM2d"
MPIEXEC="$PETSC_DIR/arch-linux2-c-opt/bin/mpiexec"
N=1

# grid: 20x20
$MPIEXEC -n $N $PETIBM2D -caseFolder ./20 \
           -sys2_pc_type gamg \
           -sys2_pc_gamg_type agg \
           -sys2_pc_gamg_agg_nsmooths 1

# grid: 60x60
$MPIEXEC -n $N $PETIBM2D -caseFolder ./60 \
           -sys2_pc_type gamg \
           -sys2_pc_gamg_type agg \
           -sys2_pc_gamg_agg_nsmooths 1

# grid: 180x180
$MPIEXEC -n $N $PETIBM2D -caseFolder ./180 \
           -sys2_pc_type gamg \
           -sys2_pc_gamg_type agg \
           -sys2_pc_gamg_agg_nsmooths 1

# grid: 540x540
$MPIEXEC -n $N $PETIBM2D -caseFolder ./540 \
           -sys2_pc_type gamg \
           -sys2_pc_gamg_type agg \
           -sys2_pc_gamg_agg_nsmooths 1

# plot the spatial convergence
python $PETIBM_DIR/scripts/python/validation/cavityConvergence.py --no-show
