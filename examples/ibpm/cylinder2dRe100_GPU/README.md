# 2D flow around cylinder (Re=100)

Run the example using 4 CPU processes and 2 GPUs:

```
export CUDA_VISIBLE_DEVICES=<idx1>,<idx2>
mpiexec -np 4 petibm-ibpm -options_left -log_view ascii:view.log
```

The simulation completes in less than 2 hours when using:
- 4 CPU processes (Intel(R) Core(TM) i7-3770 CPU @ 3.40GHz),
- 2 NVIDIA K20 GPU devices.

Plot the instantaneous force coefficients:

```
python scripts/plotDragCoefficient.py
```

We obtained a mean drag coefficient of 1.3336 and a mean lift coefficient of
0.0032 (min: -0.3488, max: 0.3487).

The plot is saved in the sub-folder `figures` of the simulation directory.

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
