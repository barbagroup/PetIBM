#### Dependencies

Ensure that the following dependencies are installed before compiling PetIBM:

* GNU C++ Compiler(`g++`) version 4.6 or above
* [PETSc](http://www.mcs.anl.gov/petsc/) version 3.7

PetIBM has been tested on:
* Ubuntu 16.04 with g++-5.4, and PETSc-3.7.4;
* Mac OS X El Capitan with g++-4.9, and PETSc-3.7.3.

---

#### GNU C++ Compiler (`g++`)

On Ubuntu, install `g++` using the following command:

    sudo apt-get install g++

Check the version of G++ installed:

    g++ --version

On Mac OS X, `g++` can be installed on Mac OS X via XCode or [Homebrew](brew.sh).

---

#### Get PETSc

Install PETSc-3.7 in debugging mode and/or optimized mode.
The debug mode is recommended during development, while the optimized mode should be used for production runs.

Get and unpack PETSc:

    cd $HOME/sfw
    mkdir -p petsc/3.7.4
    wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.4.tar.gz
    tar -xvf petsc-lite-3.7.4.tar.gz -C petsc/3.7.4 --strip-components=1
    cd petsc/3.7.4

Configure and build an debugging version of PETSc:

    export PETSC_DIR=$HOME/sfw/petsc/3.7.4
    export PETSC_ARCH=linux-dbg
    ./configure --PETSC_ARCH=$PETSC_ARCH \
        --with-cc=gcc --with-cxx=g++ --with-fc=gfortran \
        --COPTFLAGS="-O0" --CXXOPTFLAGS="-O0" --FOPTFLAGS="-O0" \
        --with-debugging=1 \
        --download-fblaslapack \
        --download-mpich
    make all
    make test

If you do not have a Fortran compiler ready, you use the following configuration:

    ./configure --PETSC_ARCH=$PETSC_ARCH \
        --with-cc=gcc --with-cxx=g++ \
        --COPTFLAGS="-O0" --CXXOPTFLAGS="-O0" \
        --with-debugging=1 \
        --download-f2cblaslapack \
        --download-mpich

Configure and build an optimized version of PETSc:

    export PETSC_DIR=$HOME/sfw/petsc/3.7.4
    export PETSC_ARCH=linux-opt
    ./configure --PETSC_ARCH=$PETSC_ARCH \
        --with-cc=gcc --with-cxx=g++ --with-fc=gfortran \
        --COPTFLAGS="-O3" --CXXOPTFLAGS="-O3" --FOPTFLAGS="-O3" \
        --with-debugging=0 \
        --download-fblaslapack \
        --download-mpich
    make all
    make test

When running the code on an external cluster, make sure that you compile PETSc with the MPI that has been configured to work with the cluster. Hence, **do not use the `--download-mpich` flag**, but instead point to the folder with the MPI compilers and executables using the `--with-mpi-dir` flag. If BLAS and LAPACK are already installed in the system, you can point to the libraries using the `--with-blas-lib` and `--with-lapack-lib` flags.

[Detailed instructions](http://www.mcs.anl.gov/petsc/documentation/installation.html) with more options to customize your installation can be found on the PETSc website. Run `./configure --help` in the PETSc root directory to list all the available configure flags.

The PETSc Users Manual and the Manual Pages can be found on their
[documentation page](http://www.mcs.anl.gov/petsc/documentation/index.html).

---

#### Get PetIBM

Create a local copy of the PetIBM repository:

    cd $HOME/sfw
    mkdir petibm && cd petibm
    git clone https://github.com/barbagroup/PetIBM.git
    export PETIBM_DIR=$HOME/sfw/petibm/PetIBM

You can set the environment variable `PETIBM_DIR` to your `$HOME/.bashrc` or `$HOME/.bash_profile` files by adding the line:

    export PETIBM_DIR=$HOME/sfw/petibm/PetIBM

and restart you terminal or  source the file:

    source $HOME/.bashrc

Note: if you are using C shell, use the `setenv` command instead of `export`.

Configure and build PetIBM using the debugging PETSc build:

    export PETSC_DIR=$HOME/sfw/petsc/3.7.4
    export PETSC_ARCH=linux-dbg
    mkdir petibm-linux-dbg && cd petibm-linux-dbg
    $PETIBM_DIR/configure \
        --prefix=$HOME/sfw/petibm/petibm-linux-dbg \
        CXX=$PETSC_DIR/$PETSC_ARCH/bin/mpicxx \
        CXXFLAGS="-g -O0 -std=c++11"
    make all
    make check
    make install

You may also want to build a optimized version:

    export PETSC_DIR=$HOME/sfw/petsc/3.7.4
    export PETSC_ARCH=linux-opt
    mkdir petibm-linux-opt && cd petibm-linux-opt
    $PETIBM_DIR/configure \
        --prefix=$HOME/sfw/petibm/petibm-linux-opt \
        CXX=$PETSC_DIR/$PETSC_ARCH/bin/mpicxx \
        CXXFLAGS="-g -O3 -std=c++11"
    make all
    make check
    make install


## Optional package: AmgXWrapper

We offer the possibility to solve the iterative systems (velocity and/or Poisson) on CUDA-capable devices using [AmgX](https://developer.nvidia.com/amgx) (version `1.2.0-build108`).

For this task, we use the package [`AmgXWrapper`](https://github.com/barbagroup/AmgXWrapper) (version `v1.0-beta` bundled with PetIBM).

AmgX depends on OpenMPI-1.8 and CUDA-6.5 and is available to CUDA Registered Developers [here](https://developer.nvidia.com/amgx).

For example, to build an optimized version of PetIBM with AmgXWrapper:

    export PETSC_DIR=$HOME/sfw/petsc/3.7.4
    export PETSC_ARCH=linux-opt
    mkdir petibm-linux-opt && cd petibm-linux-opt
    $PETIBM_DIR/configure \
        --prefix=$HOME/sfw/petibm/petibm-linux-opt \
        CXX=$PETSC_DIR/$PETSC_ARCH/bin/mpicxx \
        CXXFLAGS="-g -O3 -std=c++11" \
        --with-amgx="$HOME/sfw/amgx" \
        --with-cuda="/usr/local/cuda-6.5"
    make all
    make check
    make install


## Contributing and reporting bugs

To report bugs, please use the GitHub issue tracking system.
We are also open to pull-requests.
