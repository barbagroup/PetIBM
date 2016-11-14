# How to install and run PetIBM on Colonial One

Colonial One is a high-performance computing cluster at The George Washington University. Usage instructions are available on the [Colonial One Wiki](https://colonialone.gwu.edu/).

General instructions to build PETSc can be found [here](http://www.mcs.anl.gov/petsc/documentation/installation.html).


#### Build PETSc

Log in to Colonial One with your approved username and password.4

Get and unpack PETSc:

    cd $HOME/sfw
    mkdir -p petsc/3.7.4
    wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.4.tar.gz
    tar -xvf petsc-lite-3.7.4.tar.gz -C petsc/3.7.4 --strip-components=1
    cd petsc/3.7.4

Configure and build an debugging version of PETSc:

    export PETSC_DIR=$HOME/sfw/petsc/3.7.4
    export PETSC_ARCH=linux-dbg
    ./configure \
        --PETSC_ARCH=$PETSC_ARCH \
        --with-mpi-dir=/c1/apps/openmpi/1.8/gcc/4.9.2 \
        --with-gcc-dir=/c1/apps/gcc/4.9.2 \
        --download-fblaslapack \
        --with-shared-libraries \
        --with-debugging=1 \
        --COPTFLAGS="-O0" --CXXOPTFLAGS="-O0" --FOPTFLAGS="-O0"
    make all
    make test

Configure and build an optimized version of PETSc:

    export PETSC_DIR=$HOME/sfw/petsc/3.7.4
    export PETSC_ARCH=linux-opt
    ./configure \
        --PETSC_ARCH=$PETSC_ARCH \
        --with-mpi-dir=/c1/apps/openmpi/1.8/gcc/4.9.2 \
        --with-gcc-dir=/c1/apps/gcc/4.9.2 \
        --download-fblaslapack \
        --with-shared-libraries \
        --with-debugging=0 \
        --COPTFLAGS="-O3" --CXXOPTFLAGS="-O3" --FOPTFLAGS="-O3"
    make all
    make test

The debugging mode is recommended when developing the code. For production runs, you should the optimized mode.


#### Build PetIBM

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


Configure and build PetIBM using the optimized PETSc build:

    export PETSC_DIR=$HOME/sfw/petsc/3.7.4
    export PETSC_ARCH=linux-opt
    mkdir petibm-linux-opt && cd petibm-linux-opt
    $PETIBM_DIR/configure \
        --prefix=$HOME/sfw/petibm/petibm-linux-opt \
        CXX="/c1/apps/openmpi/1.8/gcc/4.9.2/bin/mpicxx" \
        CXXFLAGS="-g -O3 -std=c++11"
    make all
    make check
    make install

You may also want to build a debugging version:

    export PETSC_DIR=$HOME/sfw/petsc/3.7.4
    export PETSC_ARCH=linux-dbg
    mkdir petibm-linux-dbg && cd petibm-linux-dbg
    $PETIBM_DIR/configure \
        --prefix=$HOME/sfw/petibm/petibm-linux-dbg \
        CXX="/c1/apps/openmpi/1.8/gcc/4.9.2/bin/mpicxx" \
        CXXFLAGS="-g -O0 -std=c++11"
    make all
    make check
    make install


#### Submitting jobs

General instructions to submit jobs to Colonial One can be found [here](https://colonialone.gwu.edu/#Submitting_jobs_on_the_cluster).

Navigate to your simulation directory, for instance `$HOME/sfw/petibm/petibm-linux-opt/examples/2d/cylinder/Re40`, where we are going to create a submission script:

    cd $HOME/sfw/petibm/petibm-linux-opt/examples/2d/cylinder/Re40
    vim runPetIBM2d.sh

In the submission script, you should make sure to load the required modules before calling the executable. Here is what could look like a submission script if I want to simulation the 2d cylinder case at Reynolds number `40`, using `4` MPI processes on a single node, on the `short` partition, with a time-limit of `15` minutes:

```
    #!/bin/sh

    #SBATCH --job-name="cylinder2dRe40_n4"
    #SBATCH --output=log%j.out
    #SBATCH --error=log%j.err
    #SBATCH --partition=short
    #SBATCH --time=00:15:00
    #SBATCH -n 4

    module load openmpi/1.8/gcc/4.9.2
    mpiexec $HOME/sfw/petibm-linux-opt/bin/petibm2d [optional PETSc command-line arguments]
```

Then, you can submit your job to the `short` partition:

    sbatch runPetIBM2d.sh
