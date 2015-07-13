The main input that must be provided to PetIBM is the case folder, which contains the files that describe the simulation.

    > mpiexec -n 4 bin/PetIBM3d -caseFolder path/to/folder

This folder must contain the following files:

* [[flowDescription.yaml]]: This file specifies the properties of the flow - the number of dimensions (2 or 3), the kinematic viscosity of the fluid, and the boundary conditions of the domain. 
* [[cartesianMesh.yaml]]: This file describes the mesh used to solve the flow. PetIBM only makes use of cartesian meshes (that can be non-uniform).
* [[simulationParameters.yaml]]: The parameters related to the computational aspects of the simulation are provided in this file. They include the number of time steps, the size of each time step, whether the simulation must be restart from previously saved data, and the type of numerical schemes to use for the solution.

The following files are optional:
* [[bodies.yaml]]: This file lists the bodies that are present in the flow. It is not required when solving pure fluid flows, i.e. when no immersed bodies are present.
* Other body files: These files are required when bodies are described by a set of points.

All the input files make use of [YAML](http://en.wikipedia.org/wiki/YAML) to describe their associated data so that they are human-readable. Each of these files are described in detail in the pages their names link to.