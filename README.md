PetIBM - A PETSc-based Immersed Boundary Method code
====================================================

Installation instructions
-------------------------

### Dependencies

Ensure that the following dependencies are installed before compiling PetIBM:

* GNU C++ Compiler(`g++`) version 4.4 or above
* Git distributed version control system (`git`)
* PETSc version 3.4.0 or above

PetIBM has been tested and run on Ubuntu 12.04 and Mac OS X 10.9.

#### GNU C++ Compiler (`g++`)

On Ubuntu, install `g++` using the following command:

    > sudo apt-get install g++

Check the version of G++ installed:

    > g++ --version

Other development and version control tools can be installed by following the instructions under Step 1 in the [CompilingEasyHowTo](https://help.ubuntu.com/community/CompilingEasyHowTo) page on the Ubuntu Community Help Wiki. Software developers will find it useful to install all of them.

`g++` can be installed on Mac OS X by installing XCode through the App Store.

#### Git (`git`)

Install `git`. On Ubuntu, this can be done via the Terminal using the following command:

    > sudo apt-get install git-core

You can download and install `git` on OS X from [here](http://git-scm.com/download/mac).

#### PETSc

Install PETSc-3.4.4 in debug mode by running the following:

    > cd $HOME
    > mkdir src
    > cd src
    > wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.4.4.tar.gz
    > tar -xvzf petsc-lite-3.4.4.tar.gz
    > cd petsc-3.4.4
    > ./configure --download-f-blas-lapack --download-mpich
    > make all test

To install PETSc in the optimized mode, use the following flags in the configure step:

    > ./configure --with-debugging=0 --COPTFLAGS=-O3 --FOPTFLAGS=-O3 --CXXOPTFLAGS=-O3 --download-f-blas-lapack --download-mpich

Production runs should be compiled and run using the optimized mode. Use the debug mode during development.

When running the code on an external cluster, make sure that you compile PETSc with the MPI that has been configured to work with the cluster. Hence, *do not use the* `--download-mpich` *flag*, and instead specify the MPI program and compilers on the cluster using the `--with-mpiexec` and `--with-mpi-dir` flags. Run `./configure --help` to list all the available configure options.

Detailed [installation instructions](http://www.mcs.anl.gov/petsc/documentation/installation.html) can be found on the PETSc website.

The PETSc Users Manual and the Manual Pages can be found on their [documentation page](http://www.mcs.anl.gov/petsc/documentation/index.html).

### Set environment variables

Before compiling PetIBM, set the environment variables `PETSC_DIR` and `PETSC_ARCH`. *These will be displayed at the end of the configure step during PETSc installation*. You can either set them through the command line or insert the following in your `~/.bashrc` or `~/.bash_profile` files:

    > export PETSC_DIR=$HOME/src/petsc-3.4.4
    > export PETSC_ARCH=arch-linux2-c-debug

Replace the values above with the version and mode specific to your installation. Use the `setenv` command to do the same in C shell. Restart the terminal or source the file:

    > source $HOME/.bashrc


### Compiling PetIBM

Create a local copy of the PetIBM repository:

    > cd $HOME/src
    > git clone https://github.com/barbagroup/PetIBM.git

Compile all the files:
    
    > cd PetIBM
    > make

Run a test case
---------------

The following case simulates flow over an impulsively started 2-D circular cylinder at Reynolds number 40. It runs with 2 MPI processes.

    > cd $HOME/src/PetIBM    
    > make cylinderRe40

Postprocess the results:

    > python scripts/python/plot.py -folder cases/2d/cylinder/Re40 -xmin -2 -xmax 4 -ymin -2 -ymax 2

Contour plots of the x and y components of velocity, the pressure and the vorticity will be written to the folder `cases/2d/cylinder/Re40/output`.

Contact
-------

Please e-mail [Anush Krishnan](mailto:k.anush@gmail.com) if you have any questions, suggestions or feedback.