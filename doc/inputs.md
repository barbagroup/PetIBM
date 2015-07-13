A simulation directory must contain some input files in order to run the program PetIBM.

* The file [[flowDescription.yaml]] is required to prescribe the characteristics of the fluid, as well as the initial and boundary conditions of the flow.

* The file [[cartesianMesh.yaml]] is required to generate a Cartesian grid (uniform or stretched mesh).


* The file [[simulationParameters.yaml]] contains the parameters related to the computational aspects of the simulation. They include the number of time-steps, the time-increment, the saving time-step interval, whether the simulation must be restart from given initial velocity field, as well as the type of numerical schemes to discretize the equations.

The following files are optional:

* The file [[bodies.yaml]] lists the bodies that are present in the flow. It is not required when solving pure fluid flows, i.e. when no immersed body is present in the computational domain.
* Body coordinates files.

All the input files make use of [YAML](http://en.wikipedia.org/wiki/YAML) to describe their associated data so that they are human-readable. Each of these files are described in detail in the pages their names link to.
