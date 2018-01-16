# 2D lid-driven cavity flow at Re=100

Run the example:

```
mpiexec -np 1 petibm-navierstokes -options_left -log_view ascii:stdout.txt
```

The simulation completes within a few seconds when using 1 CPU process
(Intel(R) Core(TM) i7-3770 CPU @ 3.40GHz).

Plot the centerline velocity components and compare with the numerical results
from Ghia et al. (1982):

```
python scripts/plotCenterlineVelocities.py
```

The plot is saved in the subfolder `figures` on the simulation directory.
