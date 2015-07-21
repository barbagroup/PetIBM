# PetIBM - A PETSc-based Immersed Boundary Method code

[![Build Status](https://travis-ci.org/barbagroup/PetIBM.png?branch=develop)](https://travis-ci.org/barbagroup/PetIBM)

PetIBM solves the incompressible Navier-Stokes equations in two and three dimensions using the immersed boundary projection method from [Taira and Colonius](http://colonius.caltech.edu/pdfs/TairaColonius2007.pdf) (2007) and is implemented using PETSc, the Portable, Extensible Toolkit for Scientific Computation.

---

## Installation instructions

#### Dependencies

Ensure that the following dependencies are installed before compiling PetIBM:

* GNU C++ Compiler(`g++`) version 4.6 or above
* [PETSc](http://www.mcs.anl.gov/petsc/) version 3.5.0 or above (use branch `petsc-3.4-compatible` to run with PETSc 3.4)
* [Boost](http://www.boost.org) version 1.55.0 or above (no build is required)

PetIBM has been tested and run on Ubuntu 12.04, Ubuntu 14.10 and Mac OS X 10.8.

---

#### GNU C++ Compiler (`g++`)

On Ubuntu, install `g++` using the following command:

    > sudo apt-get install g++

Check the version of G++ installed:

    > g++ --version

Other development and version control tools can be installed by following the instructions under Step 1 in the
[CompilingEasyHowTo](https://help.ubuntu.com/community/CompilingEasyHowTo) page on the Ubuntu Community Help Wiki.
Software developers will find it useful to install all of them.

`g++` can be installed on Mac OS X by installing XCode through the App Store.

---

#### Get PETSc

Install PETSc-3.5.2 in debugging mode and/or optimized mode. The debug mode is recommended during development, while the optimized mode should be used for production runs.

Get and unpack PETSc:

    > cd $HOME/sfw
    > mkdir -p petsc/3.5.2
    > wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.5.2.tar.gz
    > tar -xvf petsc-lite-3.5.2.tar.gz -C petsc/3.5.2 --strip-components=1
    > cd petsc/3.5.2

Configure and build an debugging version of PETSc:

    > export PETSC_DIR=$HOME/sfw/petsc/3.5.2
    > export PETSC_ARCH=linux-dbg
    > ./configure --PETSC_ARCH=$PETSC_ARCH \
                  --download-fblaslapack \
                  --download-mpich
    > make all
    > make test

Configure and build an optimized version of PETSc:

    > export PETSC_DIR=$HOME/sfw/petsc/3.5.2
    > export PETSC_ARCH=linux-opt
    > ./configure --PETSC_ARCH=$PETSC_ARCH \
                  --with-debugging=0 \
                  --COPTFLAGS=-O3 --CXXOPTFLAGS=-O3 --FOPTFLAGS=-O3 \
                  --download-fblaslapack \
                  --download-mpich
    > make all
    > make test

When running the code on an external cluster, make sure that you compile PETSc with the MPI that has been configured to work with the cluster. Hence, **do not use the `--download-mpich` flag**, but instead point to the folder with the MPI compilers and executables using the `--with-mpi-dir` flag. If BLAS and LAPACK are already installed in the system, you can point to the libraries using the `--with-blas-lib` and `--with-lapack-lib` flags.

If you are installing PETSc-3.4, use the option `--download-f-blas-lapack` instead of `--download-fblaslapack`.

[Detailed instructions](http://www.mcs.anl.gov/petsc/documentation/installation.html) with more options to customize your installation can be found on the PETSc website. Run `./configure --help` in the PETSc root directory to list all the available configure flags.

The PETSc Users Manual and the Manual Pages can be found on their
[documentation page](http://www.mcs.anl.gov/petsc/documentation/index.html).

---

#### Get Boost

The parser used in PetIBM to go through the input files requires some Boost header files.
If the Boost library is not already installed on your machine, you may type the following command-lines:

    > cd $HOME/sfw
    > mkdir -p boost/1.57.0
    > wget http://sourceforge.net/projects/boost/files/boost/1.57.0/boost_1_57_0.tar.gz
    > tar -xvf boost_1_57_0.tar.gz -C boost/1.57.0 --strip-components=1

That's it! We only need some header files from the Boost library.

---

#### Get PetIBM

Create a local copy of the PetIBM repository:

    > cd $HOME/sfw
    > mkdir petibm && cd petibm
    > git clone https://github.com/barbagroup/PetIBM.git
    > export PETIBM_DIR=$HOME/sfw/petibm/PetIBM

You can set the environment variable `PETIBM_DIR` to your `$HOME/.bashrc` or `$HOME/.bash_profile` files by adding the line:

    > export PETIBM_DIR=$HOME/sfw/petibm/PetIBM

 and restart you terminal or  source the file:

    > source $HOME/.bashrc

Note: if you are using C shell, use the `setenv` command instead of `export`.


Configure and build PetIBM using the optimized PETSc build:

    > export PETSC_DIR=$HOME/sfw/petsc/3.5.2
    > export PETSC_ARCH=linux-opt
    > mkdir petibm-linux-opt && cd petibm-linux-opt
    > $PETIBM_DIR/configure --prefix=$HOME/sfw/petibm/petibm-linux-opt \
                            --with-boost=$HOME/sfw/boost/1.57.0 \
                            CC=$PETSC_DIR/$PETSC_ARCH/bin/mpicc \
                            CXX=$PETSC_DIR/$PETSC_ARCH/bin/mpicxx
    > make all
    > make check
    > make install

The command `make install` copies the PetIBM executables (`petibm2d` and `petibm3d`) in the sub-folder `bin` of the directory given by the flag `--prefix`.

You may also want to build a debugging version:

    > export PETSC_DIR=$HOME/sfw/petsc/3.5.2
    > export PETSC_ARCH=linux-dbg
    > mkdir petibm-linux-dbg && cd petibm-linux-dbg
    > $PETIBM_DIR/configure --prefix=$HOME/sfw/petibm/petibm-linux-dbg \
                            --with-boost=$HOME/sfw/boost/1.57.0 \
                            CC=$PETSC_DIR/$PETSC_ARCH/bin/mpicc \
                            CXX=$PETSC_DIR/$PETSC_ARCH/bin/mpicxx
    > make all
    > make check
    > make install

---

## Run a case

The sub-folder `examples` in your build directory contains several 'ready-to-run' cases.
For instance, to simulate (using 2 MPI processes) a 2d flow over an impulsively started circular cylinder at Reynolds number 40, use the follwing command-lines:

    > cd $HOME/sfw/petibm/petibm-linux-opt/examples
    > make examples
    > make cylinder2dRe40

#### Post-process the results

We also provide a series of Python scripts to post-process the numerical solution saved. Those scripts are located in the folder `$PETIBM_DIR/scripts/python`.

Those scripts have been created and tested using [Anaconda](https://store.continuum.io/cshop/anaconda/), a very useful Python distribution for scientific computing.

For instance, the Python script `plotFields2d.py` generates the contour of the velocity components, the pressure and the vorticity at saved time-steps. The following command-line displays the list of options available:

    > python $(PETIBM_DIR)/scripts/python/plotFields2d.py --help

To generate the field variable contours of the 2d cylinder case at every time-step saved in an box `[-2.0, 4.0]x[-2.0, 2.0]`, use the command-line:

    > python $(PETIBM_DIR)/scripts/python/plotFields2d.py --case 2d/cylinder/Re40 \
                                                          --bottom-left -2.0 -2.0 \
                                                          --top-right 4.0 2.0 \

Images are saved in the sub-folder `images` of the case directory.

The script `generateVTKFiles.py` allow you to generate files containing the field variables that are readable by open-source visualization tools such as [Paraview](http://www.paraview.org) or [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit).

---

## Contact

Please e-mail [Anush Krishnan](mailto:k.anush@gmail.com), [Olivier Mesnard](mailto:mesnardo@gwu.edu) if you have any questions, suggestions or feedback.

To report bugs, please use the GitHub issue tracking system.
We are also open to pull-requests.
