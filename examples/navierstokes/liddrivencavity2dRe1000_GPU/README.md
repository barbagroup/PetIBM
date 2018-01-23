# 2D lid-driven cavity flow at Re=1000 (GPU)

Run the example using 1 GPU (index 0):

```
export CUDA_VISIBLE_DEVICES=0
mpiexec -np 2 petibm-navierstokes -options_left -log_view ascii:stdout.txt
```

The simulation completes in less than 2 minutes when using:
- 2 CPU processes (Intel(R) Core(TM) i7-3770 CPU @ 3.40GHz),
- 1 NVIDIA K40 GPU device.

Plot the centerline velocity components and compare with the numerical results
from Ghia et al. (1982):

```
python scripts/plotCenterlineVelocities.py
```

The plot is saved in the subfolder `figures` of the simulation directory.
