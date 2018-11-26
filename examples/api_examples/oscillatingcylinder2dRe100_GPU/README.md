# 2D flow around inline oscillating cylinder (Re=100, KC=5)

If `Makefile` is present in this folder, you can compile the example and create the executable `oscillatingcylinder` with

```bash
make oscillatingcylinder
```

You can run the example (for example, with 4 MPI processes and 1 GPU device) with

```bash
export CUDA_VISIBLE_DEVICES=<idx1>
mpiexec -np 4 oscillatingcylinder -options_left -log_view ascii:stdout.txt
```

The simulation completes in less than 40 minutes when using:

- 2 CPU processes (Intel(R) Core(TM) i7-3770 CPU @ 3.40GHz),
- 2 NVIDIA K20 GPU devices.

Plot the drag coefficient over 4 oscillation cycles:

```bash
python scripts/plotDragCoefficient.py
```

The Python script will display the Matplotlib figure and save it into the sub-folder `figures` in the simulation directory.

Compute the vorticity field at saved time steps:

```bash
petibm-vorticity
```

Create XDMF files to visualize the data with VisIt:

```bash
petibm-createxdmf
```

The XDMF files `p.xmf` (pressure), `u.xmf` (x-velocity), `v.xmf` (y-velocity),
and `wz.xmf` (vorticity) are saved in the simulation directory and can be used
as input files to visualize data with VisIt.
