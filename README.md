# PetIBM - A PETSc-based Immersed Boundary Method code

[![Build Status](https://travis-ci.org/barbagroup/PetIBM.png?branch=develop)](https://travis-ci.org/barbagroup/PetIBM)

PetIBM solves the 2D and 3D incompressible Navier-Stokes equations using a projection method with an immersed-boundary method (IBM).
Currently, two IBMs are implemented:
* the immersed-boundary projection method (Taira and Colonius, 2007);
* and its decoupled version (Li et al., 2016).

PetIBM works on distributed-memory architectures using the [PETSc](http://www.mcs.anl.gov/petsc/) library.
We use iterative solvers for the different systems of the problem.
Currently, the system can be solved on CPUs using PETSc or on multiple GPUs using the Nvidia [AmgX](https://developer.nvidia.com/amgx) library and [AmgXWrapper](https://github.com/barbagroup/AmgXWrapper).

PetIBM runs only on Unix-based systems (no support on Windows) and was last tested on Ubuntu 16.04 and MacOS Sierra 10.12.6.

To build PetIBM, please refer to the [installation instructions](https://github.com/barbagroup/PetIBM/wiki/installation).

---

## Dependencies (last tested)

* GNU C++ compiler g++ (4.9.2, 5.4.0)
* [PETSc](https://www.mcs.anl.gov/petsc/) (3.8.1)

Optional to solve linear systems on GPUs:

* NVIDIA's CUDA compiler nvcc (8.0)
* [AmgX](https://github.com/NVIDIA/AMGX) (2.0)

Optional for pre- and post-processing scripts:

* Python (3.6)
* Numpy (1.12.1)
* H5py (2.7.0)
* Matplotlib (2.0.2)

Note: Python and libraries have been installed using Anaconda (4.4.0).

---

## Documentation

User's documentation is available on the [Wiki](https://github.com/barbagroup/PetIBM/wiki) pages of the cuIBM repository.

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