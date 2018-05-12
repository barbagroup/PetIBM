# PetIBM - A PETSc-based Immersed Boundary Method code

[![Build Status](https://travis-ci.org/barbagroup/PetIBM.png?branch=develop)](https://travis-ci.org/barbagroup/PetIBM)

PetIBM solves the 2D and 3D incompressible Navier-Stokes equations using a projection method with an immersed-boundary method (IBM).
Currently, two IBMs are implemented:
* the immersed-boundary projection method (Taira and Colonius, 2007);
* and its decoupled version (Li et al., 2016).

PetIBM works on distributed-memory architectures using the [PETSc](http://www.mcs.anl.gov/petsc/) library.
We use iterative solvers for the different systems of the problem.
Currently, the system can be solved on CPUs using PETSc or on multiple GPUs using the Nvidia [AmgX](https://developer.nvidia.com/amgx) library and [AmgXWrapper](https://github.com/barbagroup/AmgXWrapper).


PetIBM solves the incompressible Navier-Stokes equations in two and three dimensions using the immersed boundary projection method from [Taira and Colonius](http://colonius.caltech.edu/pdfs/TairaColonius2007.pdf) (2007) and is implemented using [PETSc](http://www.mcs.anl.gov/petsc/), the Portable, Extensible Toolkit for Scientific Computation.

---

## Dependencies (last tested)

* g++-4.9.2, g++-5.4.0
* [PETSc](https://www.mcs.anl.gov/petsc/) 3.8.1
* Python 2.7 (optional, for pre- and post-processing)

Note: Python and libraries have been installed using [`conda`](http://conda.pydata.org/docs/get-started.html) (4.3.1).

To build `PetIBM`, please refer to the [installation instructions](https://github.com/barbagroup/PetIBM/wiki/installation).

Users-documentation is available on the [Wiki page](https://github.com/barbagroup/PetIBM/wiki) of this repository.

---

## Contact

Please e-mail [Olivier Mesnard](mailto:mesnardo@gwu.edu) or [Pi-Yueh Chuang](mailto:pychuang@email.gwu.edu) if you have any questions, suggestions or feedback.

To report bugs, please use the GitHub issue tracking system.
We are also open to pull-requests.


_References:_
* Taira, K., & Colonius, T. (2007). The immersed boundary method: a projection approach. Journal of Computational Physics, 225(2), 2118-2137.
* Li, R. Y., Xie, C. M., Huang, W. X., & Xu, C. X. (2016). An efficient immersed boundary projection method for flow over complex/moving boundaries. Computers & Fluids, 140, 122-135.