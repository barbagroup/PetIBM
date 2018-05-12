# PetIBM Change Log


---

## 0.2.1

---

### Added

### Changed

* Update version of AmgXWrapper (from `1.0-beta2` to `1.4`).
* Update calls to AmgXWrapper routines (due to changes in API).
* AmgXWrapper-1.4 requires to use PETSc-3.8+. Therefore, when configuring with AmgX, use PETSc-3.8+.

### Fixed

* When checking for AmgX, find and keep only the first occurrence of `libamgxsh.so`.

### Removed

* When configuring with AmgX, do not check the version of OpenMPI anymore.
* Do not check for Doxygen.

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