# Example of using PetIBM API: a simple Navier-Stokes solver

This folder contains a `main.cpp` file which implements a simple Navier-Stokes solver with block LU decomposition type of projection method shown in Perot's paper (1993).
The aim here is to show a basic usage of PetIBM API for those who want to use PetIBM as a library or toolbox to build their own flow solvers.

In order to simplify the code in this example, the temporal integration is hard-coded with first-order implicit and explicit Euler methods for diffusion and convection terms, respectively.
Some flow and numerical parameters are also hard-coded, such as the viscosity and the time-step size.
Thus, this example can focus on the following concepts:

1. Reading and parsing a YAML text file that has simulation parameters in it.
2. Creating basic objects: mesh, boundary conditions (ghost points), solution, linear solvers.
3. Creating operators (matrices).
4. Using CMake to build/compile an application with PetIBM.

By following the comments in the `main.cpp`, readers should be able to understand the basic usage of PetIBM API.
For more details about classes, functions, and namespaces in PetIBM, please refer to the [API documentation](https://barbagroup.github.io/PetIBM/modules.html).

This example omits many elements that a feature-complete CFD solver should have, such as restarting, outputting transient solutions, flexible temporal integration schemes, etc.
To see how to implement these elements, please refer to the source code of the Navier-Stokes solver in the `applications` folder under PetIBM's main directory.

## Build this example

In the `CMakeLists.txt`, we assume this application does not know where PetIBM was installed, as if it is an application developed by users. Please check the comments in `CMakeLists.txt`.

To configure:

```shell
$ mkdir build
$ cd build
$ cmake -DPETIBM_DIR=<path to PetIBM installation> ../
```

Then build with:

```shell
$ make all -j <number of CPUs>
```

To keep `CMakeLists.txt` simple, we didn't write the installation procedure in it, so there is no `make install`.
The build folder is itself where the solver `liddrivencavity2d` is installed.

## Run the example

Assuming we're still in the `build` folder, to run the example on CPU with 1 MPI process:

```shell
$ ./liddrivencavity2d
```

The numerical solution will be saved in the sub-folder `output` of the simulation directory.

## Post-processing

The numerical solution is written in HDF5 format.
By the end of running the program, two new files should have been created: `solution/0002000.h5` contains the solution of the flow variables after 2000 time steps and `grid.h5` contains the coordinates of the grid points for the flow variables.

A HDF5 file just contains raw data, and most post-processing software don't know what the data mean.
Thus, we need to provide additional information to the post-processing software on how to correctly read the data.
(For example, what array in the file represents the x-coordinates of velocity grid points, or what array represents a pressure solution.)

We use XDMF format to provide such information to the post-processing software.
(Software such as ParaView and VisIt understand XDMF format.)
The present folder contains the XDMF file `postprocessing.xmf` that is specific to this example and can be open with ParaView and/or VisIt.

For those who are just interested in using PetIBM as a CFD solver, there's no need to write XDMF format.
The utility program, `petibm-createxdmf`, installed in the `bin` directory of the PetIBM installation directory, can be used to create XDMF files for the flow variables, as long as the simulations were performed with one of the solvers from the `applications` directory.
