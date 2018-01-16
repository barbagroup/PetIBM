# 2D flow around cylinder (Re=100)

Run the example using 4 CPU processes and 2 GPUs:

```
export CUDA_VISIBLE_DEVICES=1,3
mpiexec -np 4 petibm-decoupledibpm -options_left -log_view ascii:stdout.txt
```

The simulation completes in less than an hour when using:
- 4 CPU processes (Intel(R) Core(TM) i7-3770 CPU @ 3.40GHz),
- 2 NVIDIA K20 GPU devices.

Plot the instantaneous force coefficients:

```
python scripts/plotDragCoefficient.py
```

We obtained a mean drag coefficient of 1.3422 and a mean lift coefficient of
0.0022 (min: -0.3473, max: 0.3474).

The plot is saved in the subfolder `figures` of the simulation directory.

Compute the vorticity field at saved time steps:

```
petibm-vorticity
```

Create XDMF files to visualize the data with VisIt:

```
petibm-createxdmf
```

The XDMF files `p.xmf` (pressure), `u.xmf` (x-velocity), `v.xmf` (y-velocity),
and `wz.xmf` (vorticity) are saved in the simulation directory and can be used
as input files to visualize data with VisIt.
