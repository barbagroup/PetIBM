A simulation directory must contain some input files in order to run the program PetIBM.

* The file [[flowDescription.yaml]] is required to prescribe the characteristics of the fluid such as the Reynolds numbers and the boundary and initial conditions of the velocity.

* The file [[cartesianMesh.yaml]] is required to generate a structured Cartesian grid (uniform or stretched).

* The file [[simulationParameters.yaml]] contains the parameters related to the time-scheme. They include the number of time-steps, the time-increment, the saving time-step interval, the starting time-step, the name of the immersed-boundary method to use (if body present in the domain), as well as the temporal scheme for the convective and diffusive terms of the momentum equation.

The following files are optional:

* The file [[solversPetscOptions.info]] contains the parameters of the two iterative solvers. The solver for the intermediate velocity and the Poisson solver can be configured using the prefix `-velocity_` and `-poisson_`, respectively. If the file is not present, the solvers will be configured with default values.

* The file [[bodies.yaml]] lists the bodies that are present in the flow. It is not required when solving pure fluid flows, i.e. in the absence of an immersed boundary in the computational domain.
* Body coordinates files.

All the input files make use of [YAML](http://en.wikipedia.org/wiki/YAML) to describe their associated data so that they are human-readable. Each of these files are described in details in the pages their names link to.
