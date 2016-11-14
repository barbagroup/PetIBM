# PetIBM Change Log

---

## Current development

---

### Added

* Compatibility with PETSc-3.7 (last tested 3.7.4; cannot use 3.5 and 3.6 anymore).
* A change log.

### Changed

* Comply Python pre- and post-processing scripts to PEP8.
* Move `validation_data` files to `scripts/python/verification/data` and `scripts/python/validation/data`.
* GAMG parameters for the 2D lid-driven cavity flow to make it run.

### Fixed

### Removed

* Python script `restartFromSolution.py` (not finished and not necessary).