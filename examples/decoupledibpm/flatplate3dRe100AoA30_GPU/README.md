# 3D flow around inclined flat plate (Re=100, AR=2, AoA=30deg)

Run the example using 4 CPU processes and 1 GPU:

```
export CUDA_VISIBLE_DEVICES=<idx1>
mpiexec -np 4 petibm-decoupledibpm -options_left -log_view ascii:stdout.txt
```

The simulation completes in less than 20 minutes when using:
- 4 CPU processes (Intel(R) Core(TM) i7-3770 CPU @ 3.40GHz),
- 1 NVIDIA K40 GPU device.

Plot the instantaneous force coefficients:

```
python scripts/plotDragCoefficient.py
```

We obtained a mean drag coefficient of 0.7393 and a mean lift coefficient of
0.7059.

The plot is saved in the sub-folder `figures` of the simulation directory.

Compute the vorticity field at saved time steps:

```
petibm-vorticity
```

Create XDMF files to visualize the data with VisIt:

```
petibm-createxdmf
```

The XDMF files `p.xmf` for the pressure, `u.xmf`, `v.xmf`, and `w.xmf` for the
velocity components, and `wx.xmf`, `wy.xmf`, and `wz.xmf` for the vorticity
components are saved in the simulation directory and can be used as input files
to visualize data with VisIt.
