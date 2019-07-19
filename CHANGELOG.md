# PetIBM Change Log

---

## 0.4.2

---

### Added

* Add possibility to choose the order of the truncated Taylor series expansion for the inverse of the implicit operator. The order can be set in the YAML configuration file using the key `BN` under the YAML node `parameters`. (The default value is `BN: 1`.)

### Changed

* Use `-pc_factor_mat_solver_type` instead of `-pc_factor_mat_solver_package` in the configuration files for the forces solver. (This removes the deprecation warning when using PETSc-3.9+.)

### Fixed

* Probes: fix index for pressure field when using a volume probe (see PR [#145](https://github.com/barbagroup/PetIBM/pull/145)).

### Removed

---

## 0.4.1

---

### Added

### Changed

### Fixed

* Probes: monitor the correct sub-domain using indices with PETSc ordering, not natural ordering (see PR [#144](https://github.com/barbagroup/PetIBM/pull/144)).

### Removed

---

## 0.4.0

---

### Added

* Add possibility to use configuration arguments `--with-<package>-include=<path>` and `--with-<package>-lib=<path>` (instead of `--with-<package>-dir=<path>`) for yaml-cpp, gtest, AmgX, and AmgXWrapper. (`--with-<package>-dir=<path>` is still supported for those packages.)
* Add classes to monitor the solution in sub-regions or at specific points (using linear interpolation). YAML configuration should be provided in the node "probes". (See Markdown documentation `doc/markdowns/inputs.md` for details on how to use.)
* Install application header files upon make install call. The header files are installed in the include folder of the install directory. A user can now create new classes that inherits from an application class.
* Add a simple class `RigidKinematicsSolver` in the applications folder to handle cases with moving rigid bodies (with prescribed kinematics). The user should create a class that inherits from `RigidKinematcsSolver` and that implements the methods to update the coordinates and velocity of the Lagrangian points, `setCoordinatesBodies` and `setVelocityBodies`.) For example, if you want to compute the flow around a flapping wing, you just have to create a class (outside the PetIBM source directory) that inherits from the application class `RigidKinematicsSolver` and implement the case-specific methods to update the location and velocity of the wing.
* Add possibility to implement different kernels for the regularized delta function. (The 3-point delta function from Roma et al. (1999) and the 4-point delta function from Peskin (2002) are currently available in PetIBM.)

### Changed

* Upgrade to yaml-cpp-0.6.2 when downloading building yaml-cpp at PetIBM configuration time (configuration flag `--enable-yamlcpp`).
* Re-format code with clang-format. The clang-format style file is added to the repository.
* Re-write I/O functions to provide file paths that include the extension.
* Re-write application codes to simplify main functions (move PetIBM object instantiations into the init method of the application class).
* Change the scheme for the prediction of the forces at the beginning of a time step. Now using the forces from the previous time steps as prediction. (This is named "scheme 2" in Li et al., 2016.) Using this scheme avoids having to reset the vector of the Lagrangian momentum forcing every time step. (According to the authors there is little difference in the results between the different forcing schemes.)
* For each example, move the configuration files for the solvers into the sub-folder `config` of the simulation directory.
* Update the example for the 2D in-line oscillating cylinder to use the newly implemented application class `RigidKinematicsSolver`.

### Fixed


### Removed

* Remove Boost dependency; configuration does not check for Boost anymore.

---

## 0.3.1

---

JOSS revision.

### Added

* Add section in README and documentation on how to use the API.
* Add simple Navier-Stokes solver using the PetIBM API in the folder `examples/api_examples`.

### Changed

* Rename application code `tairacolonius` into `ibpm` (Immersed Boundary Projection Method); binary program now named `petibm-ibpm`.
* Make YAML converters and operator< private.
* Move example of oscillating cylinder to folder `examples/api_examples`.
* Update optional dependency `AmgXWrapper` with latest version (1.4).

### Fixed

* Remove MPI barrier in routine that creates a directory.
* Fixed memory leaks in solver programs and utility programs.

### Removed

---

## 0.3.0

---

WARNING: this is a major re-factorization of the code that is not backward compatible.

### Added

* PetIBM is now a library.
* Add `include` folder with all header files.
* Add `applications` folder with code applications that use the PetIBM library.
* Add some unit-tests.
* Add instructions on how to contribute (`CONTRIBUTING.md`).
* Add multi-body example (2 vertically aligned cylinders at Re=100).
* Add moving-body example (oscillating cylinder at Re=100).
* Require PETSc-3.8.

### Changed

* Code is now under the BSD 3-clause license.
* Move from a single PetIBM executable to a PetIBM library with a set of application codes (Navier-Stokes solver, immersed-boundary method solvers, and post-processing utilities to compute the vorticity field and generate XMDF files).
* Update user's documentation and API documentation. (Markdown files for user's documentation are now in `doc/markdowns`; API documentation is generated with Doxygen into `doc/hmtl`.)
* External dependencies such as yaml-cpp, Boost, Gtest, and AmgXWrapper are now optional; they can be installed separately or when configuring PetIBM.
* Application codes solve the systems for the velocity components, not anymore for the velocity fluxes.
* Update previous examples with new configuration file.
* Update Travis CI build to cache PETSc-3.8.2.

### Fixed


### Removed

* Remove all Python scripts from the `scripts` folder. (It contained old Python scripts that were not used anymore.)
* Cannot output the fields into PETSc binary format; only HDF5 is supported.
* Remove Wiki pages. (User's documentation and API documentation now available [here](https://barbagroup.github.io/PetIBM).)

---

## 0.2.0

---

### Added

* Compatibility with PETSc-3.7 (last tested 3.7.4; cannot use 3.5 and 3.6 anymore).
* A change log.
* Possibility to output the numerical solution in HDF5 format.
* Possibility to output the flux and/or velocity variables.
* Python script `createXMFFile.py` that generates a `.xmf` file to let VisIt know how to read the HDF5 files.
* [AmgXWrapper](https://github.com/barbagroup/AmgXWrapper) as an optional external package (version v1.0-beta).
* Examples to solve the 2D flow around a circular cylinder using AmgX.
* Possibility to have multiple immersed boundaries.
* Decoupled IBPM solver (as published in [Li et al., 2016](http://www.sciencedirect.com/science/article/pii/S0045793016302833))

### Changed

* Comply Python pre- and post-processing scripts to PEP8.
* Move `validation_data` files to `scripts/python/verification/data` and `scripts/python/validation/data`.
* GAMG parameters for the 2D lid-driven cavity flow to make it run.
* Update external package AmgXWrapper to version 1.0-beta2.

### Fixed
* Update boundary ghost points only once per time step, at the beginning of the time step.
* Bug-fix in the methods `generateBNQ()` of the class `TairaColoniusSolver` and `generateET` of the class `LiEtAlSolver`: calculate the correct widths of the computational domain (fixed index).

### Removed

* Python script `restartFromSolution.py` (not finished and not necessary).
* Non-zero initial guess as default for KSPs; the user should add `-ksp_initial_guess_nonzero true` to the command-line (or configuration file) to switch on a nonzero initial for a specific KSP.
