# 2D flow around 2 vertically aligned cylinders (Re=100)

To generate the boundary coordinates for the two cylinders:

```
python scripts/createBodies.py
```

To run the example using 4 CPU processes and 1 GPU:

```
export CUDA_VISIBLE_DEVICES=<idx1>
mpiexec -np 4 petibm-decoupledibpm -options_left -log_view ascii:stdout.txt
```

The simulation completes in less than 20 minutes when using:
- 4 CPU processes (Intel(R) Core(TM) i7-3770 CPU @ 3.40GHz),
- 1 NVIDIA K40 GPU device.

To plot the instantaneous force coefficients:

```
python scripts/plotDragCoefficient.py
```

Here are the force coefficients we obtained when averaging between 125 and 200 time units of flow simulation:

- Body 1:
    + Cd = 1.7603
    + Cl = -0.0023 ([min, max] = [-0.4888, 0.4888])
- Body 2:
    + Cd = 1.7604
    + Cl = 0.0022 ([min, max] = [-0.4886, 0.4886])

The plot is saved in the sub-folder `figures` of the simulation directory.

To compute the vorticity field at saved time steps:

```
petibm-vorticity
```

To create XDMF files to visualize the data with VisIt:

```
petibm-createxdmf
```

The XDMF files `p.xmf` (pressure), `u.xmf` (x-velocity), `v.xmf` (y-velocity),
and `wz.xmf` (vorticity) are saved in the simulation directory and can be used
as input files to visualize data with VisIt.
