# 2D flow around inline oscillating cylinder (Re=100, KC=5)

## Build the example

If a `Makefile` is present in this folder, you can compile the example and create the executable `oscillatingcylinder` with

```shell
make oscillatingcylinder
```

## Run the example

You can run the example (for example, using 4 MPI processes and 1 GPU device) with

```shell
export CUDA_VISIBLE_DEVICES=<idx1>
mpiexec -np 4 oscillatingcylinder -options_left -log_view ascii:stdout.txt
```

The simulation should complete in less than 40 minutes when using:

* 2 CPU processes (Intel(R) Core(TM) i7-3770 CPU @ 3.40GHz),
* 2 NVIDIA K20 GPU devices.

## Post-processing

Plot the drag coefficient over 4 oscillation cycles:

```shell
python scripts/plotDragCoefficient.py
```

The Python script will display the Matplotlib figure and save it into the sub-folder `figures` in the simulation directory.

Compute the vorticity field at saved time steps:

```shell
petibm-vorticity
```

Create XDMF files to visualize the data with VisIt:

```shell
petibm-createxdmf
```

The XDMF files `p.xmf` (pressure), `u.xmf` (x-velocity), `v.xmf` (y-velocity),
and `wz.xmf` (vorticity) are saved in the simulation directory and can be used
as input files to visualize data with VisIt.
