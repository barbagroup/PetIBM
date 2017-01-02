# PetIBM Change Log

---

## Current development

---

### Added

* Compatibility with PETSc-3.7 (last tested 3.7.4; cannot use 3.5 and 3.6 anymore).
* A change log.
* Possibility to output the numerical solution in HDF5 format.
* Possibility to output the flux and/or velocity variables.
* Python script `createXMFFile.py` that generates a `.xmf` file to let VisIt know how to read the HDF5 files.
* [AmgXWrapper](https://github.com/barbagroup/AmgXWrapper) as an optional external package (version v1.0-beta).
* Examples to solve the 2D flow around a circular cylinder using AmgX.

### Changed

* Comply Python pre- and post-processing scripts to PEP8.
* Move `validation_data` files to `scripts/python/verification/data` and `scripts/python/validation/data`.
* GAMG parameters for the 2D lid-driven cavity flow to make it run.

### Fixed

### Removed

* Python script `restartFromSolution.py` (not finished and not necessary).