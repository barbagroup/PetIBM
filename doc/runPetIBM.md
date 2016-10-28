Once the installation completed, the PetIBM executables (`petibm2d` and `petibm3d`) should be located into the sub-folder `bin` of your build directory.
They will used to run cases in two and three dimensions.
For simplicity, we define `BUILD_DIR` to be the build directory of PetIBM and create the variables `PETIBM2D = $BUILD_DIR/bin/petibm2d` and `PETIBM2D = $BUILD_DIR/bin/petibm3d`


## Run PetIBM using one processor

To simulate a 2d flow in serial, run the following command from your PetIBM build directory:

    $PETIBM2D -directory path/to/simulation/directory

A bunch of examples are provided in PetIBM, located in the folder `examples` of the PetIBM directory.
From your PetIBM build directory (`cd $BUILD_DIR`), go into the sub-folder `examples` and run:

    make examples

This rule copies the examples provided in PetIBM to the current build directory.

For example, to compute the 2d flow over a circular cylinder at Reynolds number 40, specify the folder in which the input files are present:

    $PETIBM2D -directory $BUILD_DIR/examples/2d/cylinder/Re40


## Run PetIBM in Parallel

To run the code in parallel, use `mpiexec` from you build directory:

    mpiexec -n 4 $PETIBM2D -directory $BUILD_DIR/examples/2d/cylinder/Re40

The above command runs the code using 4 MPI processes (specified via the `-n` flag). The folder containing `mpiexec` and other MPI compilers should be provided while installing PETSc using the command line option `--with-mpi-dir`.

If you installed MPI using the `--download-mpich` flag in the configure step of your PETSc installation, then the correct path to the executable is `$PETSC_DIR/$PETSC_ARCH/bin/mpiexec`:

    $PETSC_DIR/$PETSC_ARCH/bin/mpiexec -n 4 $BUILD_DIR/bin/petibm2d -directory $BUILD_DIR/examples/2d/cylinder/Re40


## Run using a `make` rule

If you want to run the code from within a `make` rule, include the following statements at the beginning of your makefile:

    include $(PETSC_DIR)/conf/variables
    include $(PETSC_DIR)/conf/rules

And use `$(MPIEXEC)` to run in parallel. This variable stores the correct path of the `mpiexec` that PETSc was configured with. Write your own `make` rule:

    cylinder2dRe40:
        $(MPIEXEC) -n 4 your/build/directory/bin/petibm2d -caseFolder your/simulation/directory

And call it from the command-line:

    make cylinder2dRe40
