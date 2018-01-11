# 2D flow around cylinder (Re=40)

Run the example using 2 CPU processes:

```
mpiexec -np 2 petibm-tairacolonius -options_left -log_view ascii:stdout.txt
```

The simulation takes less than 5 minutes to complete when using:
- 2 CPU processes (Intel(R) Core(TM) i7-3770 CPU @ 3.40GHz)

Plot the instantaneous drag coefficients and compares with numerical results
from Koumoutsakos and Leonard (1995):

```
python scripts/plotDragCoefficient.py
```

Compute the vorticity field at saved time steps:

```
petibm-vorticity
```

Plot the contour of the vorticity field at last saved time step:

```
python scripts/plotVorticity.py
```

The plots are saved as PNG files in the subfolder `figures` of the simulation
directory.
