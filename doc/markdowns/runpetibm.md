# Run PetIBM

Once PetIBM is installed, the executables should be located in the `bin` folder of your installation directory.
For convenience, you can prepend the `bin` directory to your PATH environment variable:

    export PATH=<petibm-installation-directory>/bin:$PATH

To run a pure Navier-Stokes simulation (no immersed boundary in the computational domain) in serial:

## Running PetIBM

    cd <simulation-directory>
    petibm-navierstokes

In the current working directory differs from the simulation directory, provide the path of the simulation directory:

    petibm-navierstokes -directory <simulation-directory>

If the YAML configuration file in not located in the simulation directory:

    petibm-navierstokes -config <config-path>

To run on two CPU processes:

    mpiexec -np 2 petibm-navierstokes

To run a simulation with the Immersed Boundary Projection Method (Taira and Colonius, 2007):

    mpiexec -np n petibm-tairacolonius

To run a simulation with a decoupled version of the Immersed Boundary Projection Method (Li et al., 2016):

    mpiexec -np n petibm-decoupledibpm


## Post-processing applications

We also provide application codes to compute the vorticity field and to create XDMF files (to visualize HDF5 fields with [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/)).

To compute the vorticity field at saved time steps:

    petibm-vorticity

It will compute the vorticity fields (`wz` for 2D configurations; `wx`, `wy`, and `wz` for 3D configurations) and write them in the time-step solution HDF5 files
The utility will also write the grid associated with each vorticity component in the file `grid.h5` (under the names `wx`, `wy`, and `wz`).

To create XDMF files:

    petibm-createxdmf

It will create XDMF files for the pressure (`p.xmf`), the velocity components (`u.xmf` and `v.xmf` for 2D configurations;`u.xmf`, `v.xmf`, and `w.xmf` for 3D configurations), and the vorticity components (`wz.xmf` for 2D configurations; `wx.xmf`, `wy.xmf`, and `wz.xmf` for 3D configurations).
These files can be open with VisIt to read, visualize, and process the field solutions.


## Running PetIBM using NVIDIA AmgX

To solve one or several linear systems on CUDA-capable GPU devices, PetIBM calls the [NVIDIA AmgX](https://github.com/NVIDIA/AMGX) library.

For example to run a simulation using 2 MPI processes and 2 GPU devices, simply use:

    mpiexec -np 2 <petibm-executable>

If you choose to run on 4 MPI processes but only using 2 GPU devices:

    export CUDA_VISIBLE_DEVICES=0,1
    mpiexec -np 2 <petibm-executable>
