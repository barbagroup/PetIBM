# 2D flow around inline oscillating cylinder (Re=100)

If `Makefile` is present in this folder, you can compile the example with

```
make oscillatingcylinder
```

This will create the executable `oscillatingcylinder` and you can run the example (with 4 MPI processes and 1 GPU device) with

```
export CUDA_VISIBLE_DEVICES=<idx1>
mpiexec -np 4 ./oscillatingcylinder -options_left -log_view ascii:stdout.txt
```

Otherwise, there is the source file `main.cpp` and it is up to the user to compile it against the PetIBM library.

`Makefile.am` and `Makefile.in` are only included for simplicity but useless.

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
