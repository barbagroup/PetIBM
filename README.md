PetIBM - A PETSc-based Immersed Boundary Method code
====================================================

Installation instructions
-------------------------

### Dependencies

Ensure that the following dependencies are installed before compiling PetIBM:

* GNU C++ Compiler(`g++`) version 4.6 or above
* Git distributed version control system (`git`)
* PETSc version 3.5.0 or above (use branch `petsc-3.4-compatible` to run with PETSc 3.4)
* [Boost](http://www.boost.org) version 1.55.0 or above

PetIBM has been tested and run on Ubuntu 12.04, Ubuntu 14.10 and Mac OS X 10.9.

#### GNU C++ Compiler (`g++`)

On Ubuntu, install `g++` using the following command:

    > sudo apt-get install g++

Check the version of G++ installed:

    > g++ --version

Other development and version control tools can be installed by following the instructions under Step 1 in the 
[CompilingEasyHowTo](https://help.ubuntu.com/community/CompilingEasyHowTo) page on the Ubuntu Community Help Wiki. 
Software developers will find it useful to install all of them.

`g++` can be installed on Mac OS X by installing XCode through the App Store.

#### Git (`git`)

Install `git`. On Ubuntu, this can be done via the Terminal using the following command:

    > sudo apt-get install git-core

You can download and install `git` on OS X from [here](http://git-scm.com/download/mac).

#### PETSc

Install PETSc-3.5.2 in the optimized mode by running the following:

Get and unpack PETSc:

    > cd $HOME/sfw
    > mkdir petsc && cd petsc
    > wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.5.2.tar.gz
    > tar xvfz petsc-lite-3.5.2.tar.gz
    > mv petsc-3.5.2 3.5.2
    > cd 3.5.2

Configure and compile PETSc:
    
    > export PETSC_DIR=$HOME/sfw/petsc/3.5.2
    > export PETSC_ARCH=linux-opt
    > ./configure --PETSC_ARCH=$PETSC_ARCH \
                  --with-debugging=0 \
                  --COPTFLAGS=-O3 --CXXOPTFLAGS=-O3 --FOPTFLAGS=-O3 \
                  --download-fblaslapack \
                  --download-mpich
    > make all
    > make test

All production runs must be compiled and run using the optimized mode.

The debug mode is recommended during development. To install PETSc in debug mode, use the following flags in the 
configure step:

    > ./configure --PETSC_ARCH=$PETSC_ARCH \
                  --download-fblaslapack \
                  --download-mpich

When running the code on an external cluster, make sure that you compile PETSc with the MPI that has been configured to 
work with the cluster. Hence, **do not use the `--download-mpich` flag**, but instead point to the folder with the MPI 
compilers and executables using the `--with-mpi-dir` flag. If BLAS and LAPACK are already installed in the system, you 
can point to the libraries using the `--with-blas-lib` and `--with-lapack-lib` flags.

If you are installing PETSc-3.4, use the option `--download-f-blas-lapack` instead of `--download-fblaslapack`.

[Detailed instructions](http://www.mcs.anl.gov/petsc/documentation/installation.html) with more options to customize 
your installation can be found on the PETSc website. Run `./configure --help` in the PETSc root directory to list all 
the available configure flags.

The PETSc Users Manual and the Manual Pages can be found on their 
[documentation page](http://www.mcs.anl.gov/petsc/documentation/index.html).


### Get PetIBM

Create a local copy of the PetIBM repository:

    > cd $HOME/sfw
    > mkdir petibm && cd petibm
    > git clone https://github.com/barbagroup/PetIBM.git

### Set environment variables

Before compiling PetIBM, set the environment variables `PETSC_DIR`, `PETSC_ARCH` and `PETIBM_DIR`.
**These will be displayed at the end of the configure step during PETSc installation**. 
You can either set them through the command line or insert the following lines in your `~/.bashrc` or `~/.bash_profile` files:

    > export PETSC_DIR=$HOME/sfw/petsc/3.5.2
    > export PETSC_ARCH=linux-opt
    > export PETIBM_DIR=$HOME/sfw/petibm/PetIBM

If you are using C shell, use the `setenv` command.
Then, restart the terminal or source the file:

    > source $HOME/.bashrc

### Configure and compile a PetIBM build:

    > mkdir petibm-linux-opt && cd petibm-linux-opt
    > $PETIBM_DIR/configure --prefix=$HOME/sfw/petibm/petibm-linux-opt \
                            CC=$PETSC_DIR/$PETSC_ARCH/bin/mpicc \
                            CXX=$PETSC_DIR/$PETSC_ARCH/bin/mpicxx
    > make all
    > make check

### Run a test case

The following case simulates flow over an impulsively started 2-D circular cylinder at Reynolds number 40. It runs with 
2 MPI processes.

    > cd $HOME/sfw/petibm/petibm-linux-opt/examples
    > make examples
    > make cylinder2dRe40

Postprocess the results:

We also provide a series of Python scripts to post-process the numerical solution saved. Those scripts are located in the folder `$PETIBM_DIR/scripts/python`.

For instance, the Python script `plotFields2d.py` generates the contour of the velocity components, the pressure and the vorticity at saved time-steps. The following command-line displays the list of options available:

    > python $(PETIBM_DIR)/scripts/python/plotFields2d.py --help

The contours, for the 2-D circular cylinder at Reynolds number 40, are generated with:

    > python $(PETIBM_DIR)/scripts/python/plotFields2d.py --case 2d/cylinder/Re40 \
                                                          --bottom-left -2.0 -2.0 \
                                                          --top-right 4.0 2.0 \

To generate .vtk files, readable by open-source visualization tools such as [Paraview](http://www.paraview.org/) or [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/), use the Python script: `generateVTKFiles.py`.

### Contact

Please e-mail [Anush Krishnan](mailto:k.anush@gmail.com), [Olivier Mesnard](mailto:mesnardo@gwu.edu) if you have any questions, suggestions or feedback.
