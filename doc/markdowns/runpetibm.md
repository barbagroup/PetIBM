# Run PetIBM

PetIBM is a library (a toolbox to solve the Navier-Stokes equations with an immersed-boundary method) that comes with several application codes.
Upon successful installation, the library (shared and/or static) and the application programs should be respectively located in the `lib` and `bin` folders of your installation directory.

Once PetIBM is installed, the libraries (shared and/or static) are located in the `lib` folder of your installation directory.
The present software package comes with 5 application codes that use the PetIBM library.
Upon successful installation, the binary executables for these applications are located in the `bin` folder of your installation directory.
For convenience, you can prepend your PATH environment variable with the `bin` directory to use the binary executables:

    export PATH=<petibm-installation-directory>/bin:$PATH

List of binary executables:
    * `petibm-navierstokes`
    * `petibm-ibpm`
    * `petibm-decoupledibpm`
    * `petibm-writemesh`
    * `petibm-vorticity`
    * `petibm-createxdmf`

These programs work for 2D **and** 3D configurations and can be run in serial or parallel (with `mpiexec` or `mpirun`).
The following sub-sections provide more details about each executable.

## Program `petibm-navierstokes`

The program solves the 2D or 3D Navier-Stokes equations using a projection method: an approximate block-LU decomposition (Perot, 1993).
The equations are fully discretized in space (second-order central differences) and in time.
At each time step, the program solves a system for the velocity field and then a Poisson system for the pressure to project the velocity onto the space of divergence-free fields.

To run the program:

    cd <simulation-directory>
    petibm-navierstokes

In the current working directory differs from the simulation directory, provide the path of the simulation directory:

    petibm-navierstokes -directory <simulation-directory>

If the YAML configuration file in not located in the simulation directory:

    petibm-navierstokes -config <config-path>

To run on two CPU processes:

    mpiexec -np 2 petibm-navierstokes


## Program `petibm-ibpm`

The program solves the 2D and 3D Navier-Stokes equations with an immersed-boundary method: the Immersed Boundary Projection Method (IBPM) proposed by Taira and Colonius (2007).
The method is based on the projection approach of Perot (1993) and is suitable to solve the flow around one or several, fixed or moving, rigid immersed boundaries.
Note that the current implementation of the method in PetIBM only works for fixed immersed boundaries.
At each time step, the program solves a system for the velocity field and then a modified Poisson system for the pressure field and Lagrangian forces.
Finally, the velocity field  is projected onto the space of divergence-free fields that satisfy the no-slip condition at the location of the Lagrangian boundary points.

To run a simulation with the IBPM:

    cd <simulation-directory>
    mpiexec -np n petibm-ibpm

You can also provide the path of the simulation directory with the command-line argument `-directory <path>` and/or the path of the YAML configuration file with `-config <path>`.


## Program `petibm-decoupledibpm`

The program solves the 2D and 3D Navier-Stokes equations with an immersed-boundary method: a decoupled version of the Immersed Boundary Projection Method proposed by Li et al. (2016).
The method extends from the work of Taira and Colonius (2007) and decouples the pressure field from the Lagrangian boundary forces.
Three linear systems are solved every time step: a system for the velocity field, a Poisson system for the pressure field, and a system for the Lagrangian forces.
The divergence-free and the no-slip constraints are imposed in a sequential fashion.

To run a simulation with a decoupled version of the IBPM:

    cd <simulation-directory>
    mpiexec -np n petibm-decoupledibpm

You can also provide the path of the simulation directory with the command-line argument `-directory <path>` and/or the path of the YAML configuration file with `-config <path>`.

## Program `petibm-writemesh`

This program is a simple (and optional) pre-processing utility that creates a structured Cartesian mesh based on the configuration provided in a given YAML file.
The gridline coordinates are then written into a HDF5 file (`grid.h5`) saved in the simulation directory.


## Program `petibm-vorticity`

This program is a post-processing utility to compute the vorticity vector field from the velocity vector field; it works for 2D and 3D configurations.

In 2D, `petibm-vorticity`:
- reads the grids for the x- and y-components of the velocity field (group names `u` and `v` in the file `grid.h5`),
- computes the grid of the vorticity scalar field,
- writes the vorticity grid in `grid.h5` (group name `wz`),
- for each saved time step:
    + reads the x- and y-components of the velocity field (group names `u` and `v` in the HDF5 solution file),
    + computes the vorticity scalar field,
    + writes the vorticity values in the HDF5 solution field (group name `wz`).

In 3D, `petibm-vorticity`:
- reads the grids for the x-, y- and z-components of the velocity field (group names `u`, `v`, and `w` in the file `grid.h5`),
- computes the grid for each component of the vorticity vector field,
- writes the vorticity grids in `grid.h5` (group name `wx`, `wy`, and `wz`),
- for each saved time step:
    + reads the x-, y-, and z-components of the velocity field (group names `u`, `v`, and `w` in the HDF5 solution file),
    + computes the vorticity vector field,
    + writes the vorticity values in the HDF5 solution field (group names `wx`, `wy`, and `wz`).

To compute the vorticity field:

    cd <simulation-directory>
    mpiexec -np n petibm-vorticity

You can also provide the path of the simulation directory with the command-line argument `-directory <path>` and/or the path of the YAML configuration file with `-config <path>`.


## Program `petibm-createxdmf`

This program is a post-processing utility which creates XDMF files for the solution variables which can be used to visualize and process HDF5 fields with [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/) (versions 2.12.3 and 2.13.1 have been tested) or [ParaView](https://www.paraview.org/) (versions 5.4.1 and 5.5.0 have been tested).
The utility works for 2D and 3D configurations.

To create XDMF files:

    petibm-createxdmf

It will create XDMF files for the pressure (`p.xmf`), the velocity components (`u.xmf` and `v.xmf` for 2D configurations;`u.xmf`, `v.xmf`, and `w.xmf` for 3D configurations), and the vorticity components (`wz.xmf` for 2D configurations; `wx.xmf`, `wy.xmf`, and `wz.xmf` for 3D configurations).


## Running PetIBM using NVIDIA AmgX

To solve one or several linear systems on CUDA-capable GPU devices, PetIBM calls the [NVIDIA AmgX](https://github.com/NVIDIA/AMGX) library.
Solving on GPU devices is supported in the programs `petibm-navierstokes`, `petibm-ibpm`, `petibm-decoupledibpm`.

For example to run a simulation using 2 MPI processes and 2 GPU devices, simply use:

    mpiexec -np 2 <petibm-executable>

If you choose to run on 4 MPI processes but only using 2 GPU devices:

    export CUDA_VISIBLE_DEVICES=0,1
    mpiexec -np 4 <petibm-executable>


## References

* Li, R. Y., Xie, C. M., Huang, W. X., & Xu, C. X. (2016). An efficient immersed boundary projection method for flow over complex/moving boundaries. Computers & Fluids, 140, 122-135.
* Perot, J. B. (1993). An analysis of the fractional step method. Journal of Computational Physics, 108(1), 51-58.
* Taira, K., & Colonius, T. (2007). The immersed boundary method: a projection approach. Journal of Computational Physics, 225(2), 2118-2137.
