# Installation

## Dependencies (current PetIBM version)

**Required**:

* GNU C++ compiler g++ (only 11.2 have been tested w/ the latest PetIBM)
* [PETSc](https://www.mcs.anl.gov/petsc/) (3.16+) with HDF5 enabled
* MPI: OpenMPI, MPICH, or Intel MPI
* [yaml-cpp](https://github.com/jbeder/yaml-cpp) (0.6.2)
* [gtest](https://github.com/google/googletest) (1.7.0)

**Optional for GPU linear solvers**:

* [AmgX](https://github.com/NVIDIA/AMGX) (v2.1.0)
* [AmgXWrapper](https://github.com/barbagroup/AmgXWrapper) (1.5)

**Optional for pre- and post-processing Python scripts**:

* Python (3.6+)
* NumPy (1.12.1)
* h5py (2.7.0)
* Matplotlib (2.0.2)

**Note**:

* MPI can be either installed during [PETSc configuration](#petsc) or installed explicitly by users.
* [yaml-cpp](https://github.com/jbeder/yaml-cpp),
  [gtest](https://github.com/google/googletest), and
  [AmgXWrapper](https://github.com/barbagroup/AmgXWrapper) can be automatically installed during PetIBM configuration or explicitly installed by users in advance.
* The current version of PetIBM has been tested with:
    * Arch Linux with g++-11.2.0 and PETSc 3.16.5
* PetIBM versions prior commit a657821 (or version 0.5.2, included) have been tested with:
    * Ubuntu 16.04 with g++-5.4, and PETSc-3.11.2
    * MacOS Sierra with g++-6.0, and PETSc-3.8.2
    * Arch Linux with g++-7.2, and PETSc-3.8.2
* Older versions of PetIBM have also been tested on the following HPC systems:
    * [GW ColonialOne](https://colonialone.gwu.edu/)
    * [Titan](https://www.olcf.ornl.gov/titan/) at ORNL

---

## GNU C++ Compiler (`g++`)

On Ubuntu, install `g++` using the following command:

    sudo apt-get install g++

Check the version of G++ installed:

    g++ --version

On Mac OS X, `g++` can be installed on Mac OS X via XCode or [Homebrew](brew.sh).

---

## PETSc

PetIBM relies on the data structures and parallel routines of the PETSc library.

Here, we provide the command-line instructions to install PETSc-3.11.3 in debugging mode and/or optimized mode.
The debug mode is recommended during development, while the optimized mode should be used for production runs.

Get and unpack PETSc:

```shell
cd $HOME/sfw
mkdir -p petsc/3.11.3
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.11.3.tar.gz
tar -xvf petsc-lite-3.11.3.tar.gz -C petsc/3.11.3 --strip-components=1
cd petsc/3.11.3
```

Configure and build an debugging version of PETSc:

```shell
export PETSC_DIR=$HOME/sfw/petsc/3.11.3
export PETSC_ARCH=linux-dbg
./configure --PETSC_DIR=$PETSC_DIR --PETSC_ARCH=$PETSC_ARCH \
    --with-cc=gcc \
    --with-cxx=g++ \
    --with-fc=gfortran \
    --COPTFLAGS="-O0" \
    --CXXOPTFLAGS="-O0" \
    --FOPTFLAGS="-O0" \
    --with-debugging=1 \
    --download-fblaslapack \
    --download-openmpi \
    --download-hdf5 \
    --download-hypre \
    --download-ptscotch \
    --download-metis \
    --download-parmetis \
    --download-superlu_dist
make all
make test
```

If you do not have a Fortran compiler ready, you use the following configuration:

```shell
./configure --PETSC_ARCH=$PETSC_ARCH \
    --with-cc=gcc \
    --with-cxx=g++ \
    --COPTFLAGS="-O0" \
    --CXXOPTFLAGS="-O0" \
    --with-debugging=1 \
    --download-f2cblaslapack \
    --download-openmpi \
    --download-hdf5 \
    --download-hypre \
    --download-ptscotch \
    --download-metis \
    --download-parmetis \
    --download-superlu_dist
```

Configure and build an optimized version of PETSc:

```shell
export PETSC_DIR=$HOME/sfw/petsc/3.11.3
export PETSC_ARCH=linux-opt
./configure --PETSC_DIR=$PETSC_DIR --PETSC_ARCH=$PETSC_ARCH \
    --with-cc=gcc \
    --with-cxx=g++ \
    --with-fc=gfortran \
    --COPTFLAGS="-O3" \
    --CXXOPTFLAGS="-O3" \
    --FOPTFLAGS="-O3" \
    --with-debugging=0 \
    --download-fblaslapack \
    --download-openmpi \
    --download-hdf5 \
    --download-hypre \
    --download-ptscotch \
    --download-metis \
    --download-parmetis \
    --download-superlu_dist
make all
make test
```

When running the code on an external cluster, make sure that you compile PETSc with the MPI that has been configured to work with the cluster.
Hence, **do not use the `--download-openmpi` flag**, but instead point to the folder with the MPI compilers and executables using the `--with-mpi-dir` flag.
If BLAS and LAPACK are already installed in the system, you can point to the libraries using the `--with-blas-lib` and `--with-lapack-lib` flags.

[Detailed instructions](http://www.mcs.anl.gov/petsc/documentation/installation.html) with more options to customize your installation can be found on the PETSc website.
Run `./configure --help` in the PETSc root directory to list all the available configure flags.

The PETSc Users Manual and the Manual Pages can be found on their
[documentation page](http://www.mcs.anl.gov/petsc/documentation/index.html).

---

## yaml-cpp

[yaml-cpp](https://github.com/jbeder/yaml-cpp) is a YAML parser in C++, used in PetIBM to parse the input configuration file.

The PetIBM configuration script gives the possibility to use a previously installed version of yaml-cpp (with the configuration argument `--with-yamlcpp-dir=<path>`, or `--with-yamlcpp-include=<path>` and `--with-yamlcpp-lib=<path>`) or to download and install yaml-cpp-0.6.2 (with `--enable-yamlcpp`).

---

## gtest

We use Google's C++ test framework [gtest](https://github.com/google/googletest) to run unit-tests after compiling PetIBM.
When configuring PetIBM, you can either use a previously installed version of gtest and provide the path of the directory of the gtest installation with `--with-gtest-dir=<path>` (or using `--with-gtest-include=<path>` and `--with-gtest-lib=<path>`) or request to download and install gtest-1.7.0 with `--enable-gtest`.

---

## PetIBM

Create a local copy of the PetIBM repository:

```shell
cd $HOME/sfw
mkdir petibm && cd petibm
git clone https://github.com/barbagroup/PetIBM.git
export PETIBM_DIR=$HOME/sfw/petibm/PetIBM
```

You can set the environment variable `PETIBM_DIR` to your `$HOME/.bashrc` or `$HOME/.bash_profile` files by adding the line:

```shell
export PETIBM_DIR=$HOME/sfw/petibm/PetIBM
```

and restart you terminal or source the file:

```shell
source $HOME/.bashrc
```

Note: if you are using C shell, use the `setenv` command instead of `export`.

Configure and build PetIBM with PETSc in debugging mode:

```shell
PETSC_DIR=$HOME/sfw/petsc/3.11.3
PETSC_ARCH=linux-dbg
mkdir petibm-linux-dbg && cd petibm-linux-dbg
$PETIBM_DIR/configure \
    --prefix=$HOME/sfw/petibm/petibm-linux-dbg \
    CXX=$PETSC_DIR/$PETSC_ARCH/bin/mpicxx \
    CXXFLAGS="-g -O0 -Wall -Wno-deprecated -std=c++14" \
    --with-petsc-dir=$PETSC_DIR \
    --with-petsc-arch=$PETSC_ARCH \
    --with-yamlcpp-dir=<path> \
    --with-gtest-dir=<path>
make all
make check
make install
```

If yaml-cpp and/or gtest are not available, you can request installation at configuration time using:

```shell
$PETIBM_DIR/configure \
    --prefix=$HOME/sfw/petibm/petibm-linux-dbg \
    CXX=$PETSC_DIR/$PETSC_ARCH/bin/mpicxx \
    CXXFLAGS="-g -O0 -Wall -Wno-deprecated -std=c++14" \
    --enable-yamlcpp \
    --enable-gtest
```

For production runs, you may want to use the optimized build of PETSc:

```shell
PETSC_DIR=$HOME/sfw/petsc/3.11.3
PETSC_ARCH=linux-opt
mkdir petibm-linux-opt && cd petibm-linux-opt
$PETIBM_DIR/configure \
    --prefix=$HOME/sfw/petibm/petibm-linux-opt \
    CXX=$PETSC_DIR/$PETSC_ARCH/bin/mpicxx \
    CXXFLAGS="-O3 -Wall -Wno-deprecated -std=c++14" \
    --with-petsc-dir=$PETSC_DIR \
    --with-petsc-arch=$PETSC_ARCH \
    --with-yamlcpp-dir=<path> \
    --with-gtest-dir=<path>
make all
make check
make install
```

Et voila!
You can now add the `bin` directory (that contains the executables) to you PATH environment variable:

```shell
export PATH=$HOME/sfw/petibm/petibm-linux-opt/bin:$PATH
```

---

## PetIBM examples

We provide some examples!
Input files are located in `$PETIBM_DIR`, but can be copied to another directory with

```shell
make copy-examples EXAMPLES_DIR=<directory>
```

`EXAMPLES_DIR` is optional and the default directory is the folder `examples` in the top build directory.

---

## Optional: using NVIDIA AmgX to solve linear systems on multiple GPUs

PetIBM implements the possibility to solve linear systems on CUDA-capable GPU devices using the NVIDIA [AmgX](https://github.com/NVIDIA/AMGX) library.

To solve linear systems on multiple GPU devices with PetIBM, the user must install CUDA and AmgX before building PetIBM and use the configuration arguments `--with-cuda-dir=<path>` for CUDA and `--with-amgx-dir=<path>` (or `--with-amgx-include=<path>` and `--with-amgx-lib=<path>`) for AmgX.

To handle the data conversion between PETSc and AmgX, we use our in-house wrapper, called [AmgXWrapper](https://github.com/barbagroup/AmgXWrapper).
When configuring PetIBM, you can either use a previously installed version of AmgXWrapper and provide the path of the directory of the AmgXWrapper installation with `--with-amgxwrapper-dir=<path>` (or using `--with-amgxwrapper-include=<path>` and `--with-amgxwrapper-lib=<path>`) or request to download and install AmgXWrapper-1.4 with `--enable-amgxwrapper`.

For example, to build an optimized version of PetIBM with CUDA, AmgX, and AmgXWrapper:

```shell
PETSC_DIR=$HOME/sfw/petsc/3.11.3
PETSC_ARCH=linux-opt
mkdir petibm-linux-opt && cd petibm-linux-opt
$PETIBM_DIR/configure \
    --prefix=$HOME/sfw/petibm/petibm-linux-opt \
    CXX=$PETSC_DIR/$PETSC_ARCH/bin/mpicxx \
    CXXFLAGS="-O3 -Wall -Wno-deprecated -std=c++14" \
    --with-petsc-dir=$PETSC_DIR \
    --with-petsc-arch=$PETSC_ARCH \
    --with-yamlcpp-dir=<path> \
    --with-gtest-dir=<path> \
    --with-cuda-dir=<path> \
    --with-amgx-dir=<path> \
    --with-amgxwrapper-dir=<path>
make all
make check
make install
```

---

## Contributing and reporting bugs

To report bugs, please use the GitHub issue tracking system.
We are also open to pull-requests.
