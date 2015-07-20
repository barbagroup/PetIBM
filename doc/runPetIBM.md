After installation, the PetIBM executables (`petibm2d` and `petibm3d`) are copied into the sub-folder `bin` of your build directory. They will used to run cases in two and three dimensions.


## Run in Serial

To simulate a 2d flow in serial, run the following command from your PetIBM build directory:

    > bin/petibm2d -caseFolder path/to/folder

For example, to compute the 2d flow over a circular cylinder at Reynolds number 40, specify the folder in which the input files for that particular case are present:

    > bin/petibm2d -caseFolder examples/2d/cylinder/Re40


## Run in Parallel

To run the code in parallel, use `mpiexec` from you build directory:

    > mpiexec -n 4 bin/petibm2d -caseFolder examples/2d/cylinder/Re40

The above command runs the code using 4 MPI processes (specified via the `-n` flag). The folder containing `mpiexec` and other MPI compilers should be provided while installing PETSc using the command line option `--with-mpi-dir`.

If you installed MPI using the `--download-mpich` flag in the configure step of your PETSc installation, then the correct path to the executable is `$PETSC_DIR/$PETSC_ARCH/bin/mpiexec`:

    > $PETSC_DIR/$PETSC_ARCH/bin/mpiexec -n 4 bin/petibm2d -caseFolder examples/2d/cylinder/Re40

Visit the [PETSc Installation](http://www.mcs.anl.gov/petsc/documentation/installation.html#mpi) page for more information on how to install PETSc with MPI.


## Run using a `make` rule

If you want to run the code from within a `make` rule, include the following statements at the beginning of your makefile:

    include $(PETSC_DIR)/conf/variables
    include $(PETSC_DIR)/conf/rules

And use `$(MPIEXEC)` to run in parallel. This variable stores the correct path of the `mpiexec` that PETSc was configured with. Write your own `make` rule:

    cylinder2dRe40:
        $(MPIEXEC) -n 4 your/build/directory/bin/petibm2d -caseFolder your/simulation/directory

And call it from the command-line:

    > make cylinder2dRe40


## Command line options

Only one command-line parameter is mandatory, which is the flag `-caseFolder` followed by the directory of the simulation, which contains the input files (`.yaml` files).

To run the case of flow in a lid-driven square cavity at Reynolds number 100 with two processes, use the following command within your build directory:

    > mpiexec -n 2 bin/petibm2d -caseFolder examples/2d/lidDrivenCavity/Re100

Any of the PETSc command line options can be passed to the program. To set the options related to the Krylov solver used to calculate the intermediate velocity, use the prefix `sys1_`. For example, to use a relative tolerance of `1.0E-06` to test for convergence:

    > mpiexec -n 2 bin/petibm2d -caseFolder examples/2d/lidDrivenCavity/Re100 -sys1_ksp_rtol 1.0E-06
    > mpiexec -n 4 bin/petibm2d -caseFolder examples/2d/lidDrivenCavity/Re100 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1

Options for the linear solver that computes the pressure and body forces can be specified with the prefix `sys2_`. The following options must be passed to use the smoothed aggregation preconditioner for this step:
