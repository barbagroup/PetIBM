# PetIBM - A PETSc-based Immersed Boundary Method code

[![Build Status](https://travis-ci.org/barbagroup/PetIBM.png?branch=develop)](https://travis-ci.org/barbagroup/PetIBM)

PetIBM implements immersed-boundary methods to solve the 2D and 3D incompressible Navier-Stokes on structured Cartesian grids using a projection approach.

Currently, two immersed boundary methods are implemented:

* Immersed Boundary Projection Method (IBPM; Taira and Colonius, 2007);
* decoupled version of the IBPM (Li et al., 2016).

PetIBM runs on distributed-memory architectures and relies on the [PETSc](http://www.mcs.anl.gov/petsc/) library for data structures and parallel routines.

PetIBM solves the different linear systems either on CPUs using PETSc KSP objects or on multiple CUDA-capable GPU devices using the NVIDIA [AmgX](https://github.com/NVIDIA/AMGX) library.
The data transfer between PETSc and AmgX is handled by the wrapper [AmgXWrapper](https://github.com/barbagroup/AmgXWrapper).

PetIBM runs only on Unix-based systems (no support on Windows) and was last tested on Ubuntu 16.04 and MacOS Sierra 10.12.6.

To build PetIBM, please refer to the [installation instructions](https://github.com/barbagroup/PetIBM/wiki/installation).

---

## Dependencies (last tested)

* GNU C++ compiler g++ (4.9.2, 5.4.0)
* [PETSc](https://www.mcs.anl.gov/petsc/) (3.8.1)

Optional to solve linear systems on GPUs:

* NVIDIA's CUDA library (8.0)
* [AmgX](https://github.com/NVIDIA/AMGX) (2.0)

Optional for pre- and post-processing scripts:

* Python (3.6)
* Numpy (1.12.1)
* H5py (2.7.0)
* Matplotlib (2.0.2)

Note: Python and libraries have been installed using Anaconda (4.4.0).

---

## Documentation

User's documentation is available on the [Wiki](https://github.com/barbagroup/PetIBM/wiki) pages of the PetIBM repository.

Doxygen API documentation is available [here](http://barbagroup.github.io/PetIBM).

---

## Papers published using PetIBM

* Mesnard, O., & Barba, L. A. (2017). _Reproducible and Replicable Computational Fluid Dynamics: It's Harder Than You Think_. Computing in Science & Engineering, 19(4), 44-55, https://doi.org/10.1109/MCSE.2017.3151254.

---

## Contact

Please e-mail [Olivier Mesnard](mailto:mesnardo@gwu.edu) or [Pi-Yueh Chuang](mailto:pychuang@email.gwu.edu) if you have any questions, suggestions or feedback.

To report bugs, please use the GitHub issue tracking system.
We are also open to pull-requests.

---

## References

* Taira, K., & Colonius, T. (2007). The immersed boundary method: a projection approach. Journal of Computational Physics, 225(2), 2118-2137.
* Li, R. Y., Xie, C. M., Huang, W. X., & Xu, C. X. (2016). An efficient immersed boundary projection method for flow over complex/moving boundaries. Computers & Fluids, 140, 122-135.