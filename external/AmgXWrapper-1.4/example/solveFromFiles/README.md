# Example: solveFromFiles
------------------------

This example reads binary files of PETSc Mat and Vecs and solve the system. It
also shows how to directly include the source files of AmgXWrapper in a PETSc
application and compile them.

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
    ${PATH_TO_THE_FOLDER_OF_EXAMPLE_solveFromFiles}
```

### Step 3: run make
---------------------

Run `make` to build the binary executable.

```bash
$ make
```

After compilation, the binary executable file `solveFromFiles` will be in the
folder `bin`, and the example of solver configuration files are in the folder
`configs`.


## Running `solveFromFiles`
--------------------------

Assume you are now in the `build` folder we used in the previous section. To
see all command line arguments for `solveFromFiles`, use

```bash
$ bin/solveFromFiles -print
```

Note, if you use `-help` instead of `-print`, you will get all PETSc arguments.

To run a test with PETSc and with 4 CPU cores, for example,

```bash
$ mpiexec -n 4 bin/solveFromFiles \
    -caseName ${CASE_NAME} \
    -mode PETSc \
    -cfgFileName ${PATH_TO_KSP_SETTING_FILE} \
    -matrixFileName ${BINARY_FILE_OF_COEFFICIENT_MATRIX} \
    -rhsFileName ${BINARY_FILE_OF_RHS_VECTOR} \
    -exactFileName ${BINARY_FILE_OF_EXACT_SOLUTION_VECTOR}
```

1. `${CASE_NAME}` is the name of this run. It's only purpose is to be printed on
    the monitor and can serve as a keyword when parsing the results in post-processing.
2. `${PATH_TO_KSP_SETTING_FILE}` is the file containing the setting to KSP solver.
3. `${BINARY_FILE_OF_COEFFICIENT_MATRIX}` is the binary-format file of the coefficient
    matrix, which is a PETSc Mat.
4. `${BINARY_FILE_OF_RHS_VECTOR}` is the binary-format file of the right-hand-side
    vector, which is a PETSc Vec.
5. `${BINARY_FILE_OF_EXACT_SOLUTION_VECTOR}` is the binary-format file of the
    exact solution vector. It's used to compute errors. It's not necessry to
    be exact solutions, if you are not using `solveFromFiles` for benchmarks
    and tests and just want to solve a system.

Note, for argument `-mode`, available options are `PETSc` and `AmgX_GPU`.
CPU version of AmgX solver is not supported in the wrapper and this example.
