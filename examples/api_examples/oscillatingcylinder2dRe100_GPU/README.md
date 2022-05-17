# 2D flow around inline oscillating cylinder (Re=100, KC=5)

## Build the example

Build procedure is the same to the basic API example, `liddrivencavity2d`. Here we just give the commands:

```shell
$ mkdir build
$ cd build
$ cmake -DPETIBM_DIR=<petibm installation path> ../
$ make all -j <number of CPUs>
```

The executable `oscillatingcylinder` will be available in the `build` folder (i.e., the current folder).

If you're interested, in the `CMakeLists.txt`, we also link against `libpetibmapps.so` using the imported target `petibm::petibmapps`.
We do so because `oscillatingcylinder` is a solver derived from the `RigidKinematicsSolver` solver class.
It is the major difference between `oscillatingcylinder` and the `liddrivencavity2d` example.

## Run the example

You can run the example (for example, using 4 MPI processes and 1 GPU device) with

```shell
export CUDA_VISIBLE_DEVICES=<idx1>
mpiexec -np 4 oscillatingcylinder -options_left -log_view ascii:view.log
```
`<idx1>` denotes the IDs of the two GPUs to be used (if there are more than two GPUs on your systems.)

The simulation should complete in less than 40 minutes when using:

* 2 CPU processes (Intel(R) Core(TM) i7-3770 CPU @ 3.40GHz),
* 2 NVIDIA K20 GPU devices.

Note, if MPI was installed using `conda`/`mamba`, then you may need to add the flag `--mca opal_cuda_support 1`, i.e., `mpiexec -np 4 --mca opal_cuda_support 1 ...`.
This is due to that the Anaconda's OpenMPI package disables the CUDA support by default.

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
