## Build and install
--------------------------

We suggest two different ways to build or install AmgXWrapper.

### 1. Compile the source code together with your application
-------------------------------------------------------------

Since AmgXWrapper doesn't have a lot of code or complex hierarchy, it's very 
easy to include the source code in the building process of your applications.
This is also the way we use AmgXWrapper in our own applications. For example,
see [PetIBM](https://github.com/barbagroup/PetIBM).

### 2. Build and install AmgXWrapper with CMake
------------------------------------------------

If you want to use AmgXWrapper as a library, we provide CMakeList.txt for you to
build and install it to your preferred location.

#### Step 1

At any location you like, make a temporary folder to held compiled files and 
then go into that folder. For example,

```bash
$ mkdir ${HOME}/build
$ cd ${HOME}/build
```

#### Step 2

Run CMake.

```bash
$ cmake \
    -D CMAKE_INSTALL_PREFIX=${PATH_TO_WHERE_YOU_WANT_TO_INSTALL_AmgXWrapper} \
    -D PETSC_DIR=${PATH_TO_PETSC} \
    -D PETSC_ARCH=${THE_BUILD_OF_PETSC_YOU_WANT_TO_USE} \
    -D CUDA_DIR=${PATH_TO_CUDA} \
    -D AMGX_DIR=${PATH_TO_AMGX} \
    ${PATH_TO_AmgXWrapper_SOURCE}
```

Other available or possible CMake arguments are:

* `CMAKE_C_COMPILER`: the default is `mpicc`. Given that the current version of 
  AmgX only support OpenMPI 1.8, we are not sure if other MPI implementations
  will work. But you can try. If you simply want to change underlying C compiler,
  you can change the setting with environmental variable `OMPI_MPICC`. See 
  [Compiling MPI applications](https://www.open-mpi.org/faq/?category=mpi-apps).

* `CMAKE_CXX_COMPILER`: the default is `mpicxx`. Given that the current version
  of AmgX only supports OpenMPI 1.8, we are not sure if other MPI implementations
  will work. But you can try. If you simply want to change underlying C++ compiler,
  you can change the setting with environmental variable `OMPI_MPICXX`. See 
  [Compiling MPI applications](https://www.open-mpi.org/faq/?category=mpi-apps).

* `CMAKE_BUILD_TYPE`: the default is `RELEASE`. `RELEASE` mode will use flags
  `-O3 -DNDEBUG`. You can use build type `DEBUG` for using flag `-g`.

* `CMAKE_CXX_FLAGS` and `CMAKE_C_FLAGS`: set your own flags here, if you want.

* `BUILD_SHARED_LIBS`: the default is `ON`. This will create a shared library 
  `libAmgXWrapper.so`. To create a static library (libAmgXWrapper.a) only, set
  this argument to `OFF`.

#### Step 3

Build the library.

```bash
$ make
```

#### Step 4 (optional)

If you have Doxygen installed, and CMake is able to find Doxygen, you can build
the API documentation with

```bash
$ make doc
```

Both html and latex versions will be built.

#### Step 5

Install the library, header, and documents to your preferred location.

```bash
$ make install
```

After installation, you can safely delete the temporary build folder. For example,
here we can delete `${HOME}/build`.

Under the installation location, `lib` (or `lib64`, depends on your system)
contains the library, libAmgXWrapper.so or libAmgXWrapper.a. The header file 
AmgXSolver.hpp will be in the folder `include`, and the documents (if available)
will be in the folder `doc`.
