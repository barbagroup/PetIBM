# 2D flow around inline oscillating cylinder (Re=100)

To compile this example:

```
make oscillatingcylinder
```

Run the example using 4 CPU processes and 1 GPU:

```
export CUDA_VISIBLE_DEVICES=<idx1>
mpiexec -np 4 petibm-decoupledibpm -options_left -log_view ascii:stdout.txt
```

The simulation completes in less than 30 minutes when using:
- 4 CPU processes (Intel(R) Core(TM) i7-3770 CPU @ 3.40GHz),
- 1 NVIDIA K40 GPU device.

Plot the instantaneous force coefficients:

```
python scripts/plotDragCoefficient.py
```

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
