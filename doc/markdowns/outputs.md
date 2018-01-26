# Output files

All PetIBM output files created during the course of the simulation are located in the simulation directory and its sub-folders. Here is a list of such files:

* `grid.h5`: HDF5 file containing the gridline stations in the x, y, and z directions (datasets `x`, `y`, and `z`) for each fluid variables (group `p` for the pressure and groups `u`, `v`, and `w` for the velocity components).

* `forces.txt`: ASCII file that contains the hydrodynamic forces acting on the immersed body at each time-step. The first column contains the time values. The next two columns (or three columns for 3D runs) contains the forces in the x and y directions (and in the z direction for 3D runs). If there is a second immersed boundary in the domain, the force columns will be append to the right. (Note that this file does not exist for pure Navier-Stokes simulations, i.e. when there is no immersed boundary in the computation domain.)

* `iterations.txt`: ASCII file reporting the number of iterations to converge and the residuals for each linear solver: velocity solver, Poisson solver, and forces solver (when using the decoupled version of the immersed-boundary projection method). The first column contains the time-step index; the second and third columns contains the number of iterations to converge and the residuals for the first linear solver (velocity); etc.

* The subfolder `solution` contains the numerical solution in HDF5 format at certain time steps. (The frequency of saving is prescribed in the YAML configuration file.) For example, the numerical solution after 100 time steps is saved in the file `solution/0000100.h5`. The velocity field, the pressure field, the boundary forces (when body present in the domain), and the convection and diffusion terms (for restarting) are saved in the same time-step solution file.
