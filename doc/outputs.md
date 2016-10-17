All PetIBM output files created during the course of the simulation are located in the simulation directory and its sub-folders. Here is a list of such files:

* `grid.txt`: ASCII file containing the points along a gridline of the structured Cartesian mesh in each direction. The first line stores the number of cells in each direction. Then, coordinates are stored in ascending order along each direction, starting with the x-direction, immediately followed by the y-direction, and finally the z-direction (for 3d problem).

* `forces.txt`: this files stores the hydrodynamic forces acting on the immersed body at each time-step. The first column contains the time values. Following columns contains the force component in the x-direction, the y-direction (and the z-direction for a 3D simulation). This file is only generated when a body is immersed in the computational domain.

* `iterationCounts.txt`: this file consists of three columns - the first columns contains the time-step index; the second and third columns contain the number of iterations at a time-step required by the iterative solvers to converge (second column for the velocity solver, third one for the Poisson solver).

* Every given time-step interval (see `nsave` in the input file `simulationParameters.yaml`), the numerical solution if saved in a sub-folder whose name is the time-step index. Each sub-folder contains the following files:
  * `qx.dat`, `qy.dat`, and `qz.dat`: PETSc binary files that store the velocity flux components through each cell faces in the domain.
  * `phi.dat`: PETSc binary file that stores the values of discrete pressure field.
  * `fTilde.dat`: this file is created only if an immersed boundary is present in the computational domain. It is a PETSc binary file that stores a vector of the rescaled body forces calculated at every boundary point.
  * Along with the above files, the corresponding `.info` files are also written, which are used by PETSc.
