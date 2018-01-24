# PetIBM - A PETSc-based Immersed Boundary Method code

[Build Status]: https://travis-ci.org/barbagroup/PetIBM.png?branch=develop
[![Build Status]](https://travis-ci.org/barbagroup/PetIBM)

PetIBM implements immersed-boundary methods to solve 2D and 3D incompressible
Navier-Stokes on stretched Cartesian grids using a projection approach.

Currently, two immersed boundary methods are implemented:

* Immersed Boundary Projection Method (IBPM; Taira and Colonius, 2007);
* decoupled version of the IBPM (Li et al., 2016).

With object-oriented design, the objects and classes in PetIBM can be re-used to 
develop other solvers easily, as long as the numerical methods used can fit into
Perot's framework (Perot, 1993; Chang et. al, 2002). 
See [Doxygen pages](http://barbagroup.github.io/PetIBM/modules.html) for API
manual.

PetIBM relies on the [PETSc](http://www.mcs.anl.gov/petsc/) library for data 
structures and parallel routines. 
Linear systems can be solved either on CPUs using PETSc KSP objects or on multiple 
CUDA-capable GPU devices using the NVIDIA [AmgX](https://github.com/NVIDIA/AMGX) 
library. 
Data transfers between PETSc and AmgX are handled by 
[AmgXWrapper](https://github.com/barbagroup/AmgXWrapper).

PetIBM runs only on Unix-based systems (no support on Windows) and was last 
tested on Ubuntu 16.04, MacOS Sierra 10.12.6, and Arch Linux. 
PetIBM was also tested on the following HPC systems: 
[GW ColonialOne](https://colonialone.gwu.edu/) and 
[Titan at ORNL](https://www.olcf.ornl.gov/titan/).

Please see [Documentation](#documentation) for more details.

---

## Features

PetIBM supports:
* multiple bodies,
* moving bodies with prescribed motions,
* 2D and 3D stretched mesh,
* distributed-memory architectures,
* multiple GPUs on a single node,
* GPU clusters, and
* parallel HDF5 I/O.

---

## Documentation

* [Markdown files for quick start](doc/markdowns)
    * [Dependencies and Installation](doc/markdowns/installation.md)
    * [Run PetIBM](doc/markdowns/runpetibm.md)
    * [Input files](doc/markdowns/inputs.md)
    * [Output files](doc/markdowns/outputs.md)
    * [2D Examples](doc/markdowns/examples2d.md)
    * [3D Examples](doc/markdowns/examples3d.md)
* [Online API manual](http://barbagroup.github.io/PetIBM)
* [Change Log](CHANGELOG.md)
* [Contributing](CONTRIBUTING.md)

Offline API manual can be generated with 
[Doxygen](http://www.stack.nl/~dimitri/doxygen/).

---

## Papers published using PetIBM

* Mesnard, O., & Barba, L. A. (2017). _Reproducible and Replicable Computational
  Fluid Dynamics: It's Harder Than You Think_. Computing in Science & Engineering,
  19(4), 44-55, https://doi.org/10.1109/MCSE.2017.3151254.

---

## Contact

Please e-mail [Olivier Mesnard](mailto:mesnardo@gwu.edu) or 
[Pi-Yueh Chuang](mailto:pychuang@gwu.edu) if you have any questions, suggestions
or feedback.

To report bugs, please use the GitHub issue tracking system.
We are also open to pull-requests.

---

## References

* Perot, J. B. (1993). *An analysis of the fractional step method*. Journal of 
  Computational Physics, 108(1), 51-58.
* Chang, W., Giraldo, F., & Perot, B. (2002). *Analysis of an exact fractional 
  step method*. Journal of Computational Physics, 180(1), 183-199.
* Taira, K., & Colonius, T. (2007). *The immersed boundary method: a projection 
  approach*. Journal of Computational Physics, 225(2), 2118-2137.
* Li, R. Y., Xie, C. M., Huang, W. X., & Xu, C. X. (2016). *An efficient immersed
  boundary projection method for flow over complex/moving boundaries*. 
  Computers & Fluids, 140, 122-135.
