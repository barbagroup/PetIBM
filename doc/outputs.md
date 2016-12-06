All PetIBM output files created during the course of the simulation are located in the simulation directory and its sub-folders. Here is a list of such files:

* `grid.txt`: ASCII file containing the points along a gridline of the structured Cartesian mesh in each direction. The first line stores the number of cells in each direction. Then, coordinates are stored in ascending order along each direction, starting with the x-direction, immediately followed by the y-direction, and finally the z-direction (for 3d problem).

* `forces.txt`: this files stores the hydrodynamic forces acting on the immersed body at each time-step. The first column contains the time values. Following columns contains the force component in the x-direction, the y-direction (and the z-direction for a 3D simulation). This file is only generated when a body is immersed in the computational domain.

* `iterationCounts.txt`: this file consists of three columns - the first columns contains the time-step index; the second and third columns contain the number of iterations at a time-step required by the iterative solvers to converge (second column for the velocity solver, third one for the Poisson solver).

* The sub-folder `grids` is generated **only** when HDF5 is chosen as output format (by adding the line `outputFormat: hdf5` to your input file `simulationParameters.yaml`). The folder contains files that store the locations in the computational domain of a cell-centered quantity (`cell-centered.h5`) and of the vector components of a staggered quantity (`staggered-x.h5`, `staggered-y.h5`, and `staggered-z.h5` for 3D runs).

* Every given time-step interval (see `nsave` in the input file `simulationParameters.yaml`), the numerical solution if saved in a sub-folder whose name is the time-step index. The content of these folders depends on the type of output requested (we support PETSc binary and HDF5 formats) and the variables you choose to save (velocity and/or flux components; the pressure field is always saved).

For example, if you choose the PETSc binary format and decide to save the flux components, each sub-folder will contain the following files: 
  - `qx.dat`, `qy.dat`, and `qz.dat` (for 3D runs): files that store the velocity flux components through each cell faces in the domain;
  - `phi.dat`: file that stores the values of the discrete pressure field;
  - `fTilde.dat` (if immersed boundary present in the flow): file that stores a vector of the rescaled body forces calculated at every boundary point;
  - a `.info` file for each quantity saved that is used by PETSc to read the files.

In case where you choose a HDF5 format and decide to save the velocity components, a time-step sub-folder will contain the following files: `phi.h5`, `ux.h5`, `uy.h5`, `uz.h5` (for 3D runs), `fTilde.h5` (if you have an immersed boundary).