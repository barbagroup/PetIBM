# Example of using PetIBM API: a simple Navier-Stokes solver

This folder contains a `main.cpp` file which implements a simple Navier-Stokes
solver with block LU decomposition type of projection method shown in Perot's paper
(1993). The aim here is to demonstration the basic usage of PetIBM API for those
who want to use PetIBM as a library or toolbox to build their own flow solvers.

In order to simplify the code in this example, the temporal integration is 
hard-coded with 1st order implicit and explicit Euler methods for diffusion and
convection terms, respectively. And some flow and numerical parameters are also
hard coded, such as viscosity and time-step size. So this example can focus on
the following concepts:

1. Reading and parsing a YAML text file that has simulation parameters in it.
2. Creating basic objects: mesh, boundary conditions (ghost points), solution,
   linear solvers.
3. Creating operators (matrices).

Following the comments in the `main.cpp`, readers should be able to understand
the basic usage of PetIBM API. For more details of the classes, functions, and
namespaces in PetIBM, please refer to 
[API documentation](https://barbagroup.github.io/PetIBM/modules.html).

This example omits many elements that a complete CFD solver should have, such as
restarting, outputting transient solutions, flexible temporal integration schemes,
etc. To see how to implement these elements, please refer to the source code of
the Navier-Stokes solver in the application folder under the PetIBM source folder.

## Build this example

We provide a shell script, `build.sh` for readers to build this example. It's
not our intention to demonstration how to compile and build software under Linux
environment, so we don't provide any Makefile or CMake files in this example.
But the shell script should be enough to show how to link against PetIBM library.

In a nutshell, when using PetIBM as a library, we also need to link against PETSc
library, AmgX library (if used), and other user-provided third-party libraries (
for example, if users provide their own yaml-cpp, instead of using the one 
provided by PetIBM configuration system).

In order to use `build.sh`, readers should first modify the variables in the
first section of the scrip to reflect the path of PETSc and PetIBM installation. 
The remaining sections in the script show how to grab the configurations from 
PETSc library and how to link PetIBM library. Then, simply run `sh ./build.sh` 
to build this example.

## Post-processing

The data file of the solution output is in HDF5. After building this
example and running it, the solver should create two HDF5, `grid.h5` and a 
`0002000.h5` under solution folder. They contain data of grid points and solutions
at the 2000th time-step.

An HDF5 file is just a file containing raw data, and most post-processing software 
don't know what the data in an HDF5 file means. So we need to tell post-processing 
software what data in as HDF5 file mean what. For example, what array in the file 
represents the x-coordinates of velocity grid points, or what array represents a 
pressure solution.

To do this, we utilize XDMF format to provide such information to post-processing
software. Software such as ParaView and VisIt understand XDMF format. In this
example, we also provide an XDMF file `postprocessing.xmf` as an example of
how to write an XDMF file. 

For those who are just interested in using PetIBM as a CFD solver, there's no 
need to write XDMF format. The utility provided in the application folder, 
`petibm-createxdmf` can create XDMF files for them as long as the simulations are 
done with the solvers in the application folder. 
