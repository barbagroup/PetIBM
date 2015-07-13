After installation, two executables `PetIBM2d` and `PetIBM3d` are generated in the folder `bin/`. These are used to simulate 2-D and 3-D flows respectively.

## Run in Serial

To simulate a 2-D flow in serial, run the following command from the PetIBM root directory:

    > bin/PetIBM2d -caseFolder path/to/folder

For example, to compute flow over a circular cylinder at Reynolds number 40, specify the folder in which the input files for that particular case are present:

    > bin/PetIBM2d -caseFolder cases/2d/cylinder/Re40

## Run in Parallel

To run the code in parallel, use `mpiexec`:

    > mpiexec -n 4 bin/PetIBM2d -caseFolder cases/2d/cylinder/Re40

The above command runs the code using 4 MPI processes (specified via the `-n` flag). The folder containing `mpiexec` and other MPI compilers should be provided while installing PETSc using the command line option `--with-mpi-dir`.

If you installed MPI using the `--download-mpich` flag in the configure step of your PETSc installation, then the correct path to the executable is `${PETSC_DIR}/${PETSC_ARCH}/bin/mpiexec`:

    > ${PETSC_DIR}/${PETSC_ARCH}/bin/mpiexec -n 4 bin/PetIBM2d -caseFolder cases/2d/cylinder/Re40

Visit the [PETSc Installation](http://www.mcs.anl.gov/petsc/documentation/installation.html#mpi) page for more information on how to install PETSc with MPI.

## Run using a `make` rule

If you want to run the code from within a makefile, include the following statements at the beginning of the file:

    include ${PETSC_DIR}/conf/variables
    include ${PETSC_DIR}/conf/rules

And use `${MPIEXEC}` to run in parallel. This variable stores the correct path of the `mpiexec` that PETSc was configured with. Write your own `make` rule:

    cylinderRe40:
        ${MPIEXEC} -n 4 bin/PetIBM2d -caseFolder cases/2d/cylinder/Re40

And call it from the command line:

    > make cylinderRe40

## Command line options

Only one command line option is mandatory, and that is the `-caseFolder` option. This gives the location of the input files for the simulation. To run the case of flow in a lid-driven square cavity at Reynolds number 100 with two processes, use the following command:

    > mpiexec -n 2 bin/PetIBM2d -caseFolder cases/2d/lidDrivenCavity/Re100

Any of the PETSc command line options can be passed to the program. To set the options related to the Krylov solver used to calculate the intermediate velocity, use the prefix `sys1_`. For example, to use a relative tolerance of 1e-6 to test for convergence:

    > mpiexec -n 2 bin/PetIBM2d -caseFolder cases/2d/lidDrivenCavity/Re100 -sys1_ksp_rtol 1e-6

Options for the linear solver that computes the pressure and body forces can be specified with the prefix `sys2_`. The following options must be passed to use the smoothed aggregation preconditioner for this step:

    > mpiexec -n 4 bin/PetIBM2d -caseFolder cases/2d/lidDrivenCavity/Re100 -sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1