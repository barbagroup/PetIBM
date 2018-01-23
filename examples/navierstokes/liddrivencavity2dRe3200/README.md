# 2D lid-driven cavity flow at Re=3200

Run the example:

```
mpiexec -np 4 petibm-navierstokes -options_left -log_view ascii:stdout.txt
```

The simulation completes in about 6 minutes when using 4 CPU processes
(Intel(R) Core(TM) i7-3770 CPU @ 3.40GHz).

Plot the centerline velocity components and compare with the numerical results
from Ghia et al. (1982):

```
python scripts/plotCenterlineVelocities.py
```

The plot is saved in the subfolder `figures` of the simulation directory.
