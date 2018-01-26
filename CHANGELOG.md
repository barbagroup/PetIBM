# PetIBM Change Log

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