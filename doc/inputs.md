PetIBM reads input parameters from a configuration file typically called [[config.yaml]].

The simulation directory (where the numerical solution will be saved) should contain a YAML configuration file.

In the presence of an immersed boundary, the path of the file that contains the boundary points is prescribed in the configuration file (using a relative path).

The parameters for the linear solvers can be stored in text files and their relative paths should be written in the YAML configuration file.
