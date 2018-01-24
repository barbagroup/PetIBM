# 2D flow around cylinder (Re=3000)

Run the example using 4 CPU processes and 2 GPUs:

```
export CUDA_VISIBLE_DEVICES=<idx1>,<idx2>
mpiexec -np 4 petibm-decoupledibpm -options_left -log_view ascii:stdout.txt
```

The simulation completes in less than 20 minutes to complete when using:
- 4 CPU processes (Intel(R) Core(TM) i7-3770 CPU @ 3.40GHz),
- 2 NVIDIA K20 GPU devices.

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

The plots are saved as PNG files in the sub-folder `figures` of the simulation
directory.
