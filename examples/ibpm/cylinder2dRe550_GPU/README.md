# 2D flow around cylinder (Re=550)

Run the example using 2 CPU processes and 1 GPU device:

```
export CUDA_VISIBLE_DEVICES=0
mpiexec -np 2 petibm-tairacolonius -options_left -log_view ascii:stdout.txt
```

The simulation completes in less than 5 minutes when using:
- 2 CPU processes (Intel(R) Core(TM) i7-3770 CPU @ 3.40GHz),
- 1 NVIDIA K40 GPU device.

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
