# Example: poisson
------------------------

This example solves a Poisson equation in 2D and 3D. We use standard central
difference with uniform grid within each direction. The exact solutions are
known, so `poisson` can also serve as a kind of test.

The Poisson equation is
<center>Δu = -8π cos(2πx) cos(2πy)</center>
for 2D. And for 3D, it is
<center>Δu = -12π cos(2πx) cos(2πy) cos(2πz)</center>

The exact solutions are
<center>Δu = cos(2πx) cos(2πy)</center>
for 2D, and
<center>Δu = cos(2πx) cos(2πy) cos(2πz)</center>
for 3D.

The boundary condition is all-Neumann BC, except that we pin a point as a reference
point (i.e., apply Dirichlet BC to that point) to avoid singular matrix. We
choose the point represented by the first row in matrix A as our reference point.

You can solve other Poisson equations with known exact solutions by modifying
the hard-coded equation in functions `generateRHS` and `generateExt`. The
computational domain is hard-coded as `Lx=Ly=Lz=1`. You are free to change them.

Most of code is regarding to setting up the Poisson problem with PETSc. If you
are already familiar with PETSc, you should be able to quickly find out how to
solve PETSc linear systems with `AmgXWrapper`.

## Build
---------

We use CMake to build this example.

### Step 1: create a directory for build
-----------------------------------------

Create a folder to accommodate the compiled binary executable and go into it.
For example:

```bash
$ mkdir ${HOME}/build
$ cd ${HOME}/build
```

### Step 2: generate Makefile with CMake
-----------------------------------------

```bash
$ cmake \
    -DPETSC_DIR=${PATH_TO_PETSC_INSTALLATION} \
    -DPETSC_ARCH=${NAME_OF_PETSC_BUILD_DESIRED} \
    -DCUDA_DIR=${PATH_TO_CUDA6.5} \
    -DAMGX_DIR=${PATH_TO_AMGX} \
    ${PATH_TO_THE_FOLDER_OF_EXAMPLE_poisson}
```

### Step 3: run make
---------------------

Run `make` to build the binary executable.

```bash
$ make
```

After compilation, the binary executable file `example` will be in the
folder `bin`, and the samples of solver configuration files are in the folder
`configs`.


## Running `poisson`
--------------------------

Assume you are now in the `build` folder we used in the previous section. To
see all command line arguments for `poisson`, use

```bash
$ bin/poisson -print
```

Note, if you use `-help` instead of `-print`, you will get all PETSc arguments.

To run a test with PETSc and with 4 CPU cores to solve the 2D Poisson, for example,

```bash
$ mpiexec -n 4 bin/poisson \
    -caseName ${CASE_NAME} \
    -mode PETSc \
    -cfgFileName ${PATH_TO_KSP_SETTING_FILE} \
    -Nx ${NUMBER_OF_CELL_IN_X_DIRECTION} \
    -Ny ${NUMBER_OF_CELL_IN_Y_DIRECTION}
```

1. `${CASE_NAME}` is the name of this run. It's only purpose is to be printed on
    the monitor and can serve as a keyword when parsing the results in post-processing.
2. `${PATH_TO_KSP_SETTING_FILE}` is the file containing the setting to KSP solver.

Note, for argument `-mode`, available options are `PETSc` and `AmgX_GPU`.
CPU version of AmgX solver is not supported in the wrapper and this example.

To run a test with AmgX and with 4 CPU cores plus GPUs to solve the 3D Poisson,
for example,

```bash
$ mpiexec -n 4 bin/poisson \
    -caseName ${CASE_NAME} \
    -mode AmgX_GPU \
    -cfgFileName ${PATH_TO_KSP_SETTING_FILE} \
    -Nx ${NUMBER_OF_CELL_IN_X_DIRECTION} \
    -Ny ${NUMBER_OF_CELL_IN_Y_DIRECTION} \
    -Nz ${NUMBER_OF_CELL_IN_Z_DIRECTION}
```

The number of GPUs used will depend on the available GPUs on the machines. To
limit the number of GPUs used, use the environmental variable `CUDA_VISIBLE_DEVICES`.
For example, to only use the first two GPUs on the machine which has four GPUs, do

```bash
$ export CUDA_VISIBLE_DEVICES=0,1
```

before running the program.
