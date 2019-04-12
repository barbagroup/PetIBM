# PetIBM - toolbox and applications of the immersed-boundary method on distributed-memory architectures

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/barbagroup/PetIBM/raw/master/LICENSE)
[![Travis](https://img.shields.io/travis/barbagroup/PetIBM/develop.svg?logo=travis)](https://travis-ci.org/barbagroup/PetIBM)
[![Docs](https://img.shields.io/badge/docs-0.4-brightgreen.svg)](https://barbagroup.github.io/PetIBM)
[![DOI](https://img.shields.io/badge/JOSS-10.21105%2Fjoss.00558-brightgreen.svg)](https://doi.org/10.21105/joss.00558)
[![CITE_BIB](https://img.shields.io/badge/Cite%20PetIBM-bibtex-blue.svg)](https://www.doi2bib.org/bib/10.21105/joss.00558)

PetIBM implements immersed-boundary methods to solve 2D and 3D incompressible Navier-Stokes on stretched Cartesian grids using a projection approach.

Currently, two immersed boundary methods are implemented:

* Immersed Boundary Projection Method (IBPM; Taira and Colonius, 2007);
* decoupled version of the IBPM (Li et al., 2016).

With object-oriented design, the objects and classes in PetIBM can be re-used to develop other solvers easily, as long as the numerical methods used can fit into Perot's framework (Perot, 1993; Chang et. al, 2002). 
See [Doxygen pages](https://barbagroup.github.io/PetIBM/modules.html) for API manual.

PetIBM relies on the [PETSc](http://www.mcs.anl.gov/petsc/) library for data structures and parallel routines. 
Linear systems can be solved either on CPUs using PETSc KSP objects or on multiple CUDA-capable GPU devices using the NVIDIA [AmgX](https://github.com/NVIDIA/AMGX) library.
Data transfers between PETSc and AmgX are handled by [AmgXWrapper](https://github.com/barbagroup/AmgXWrapper).

PetIBM runs only on Unix-based systems (no support on Windows) and was last tested on Ubuntu 16.04, MacOS Sierra 10.12.6, and Arch Linux.
PetIBM was also tested on the following HPC systems: [GW ColonialOne](https://colonialone.gwu.edu/) and [Titan at ORNL](https://www.olcf.ornl.gov/titan/).

Please see [Documentation](#documentation) for more details.

---

## Features

PetIBM supports:
* multiple immersed bodies,
* moving bodies with prescribed kinematics,
* 2D and 3D stretched Cartesian meshes,
* distributed-memory architectures,
* multiple GPUs on a single node,
* GPU clusters, and
* HDF5 I/O.


---

## Documentation

* [Quick Start](doc/markdowns)
    * [Dependencies and Installation](doc/markdowns/installation.md)
    * [Run PetIBM](doc/markdowns/runpetibm.md)
    * [Input files](doc/markdowns/inputs.md)
    * [Output files](doc/markdowns/outputs.md)
    * [2D Examples](doc/markdowns/examples2d.md)
    * [3D Examples](doc/markdowns/examples3d.md)
    * [Use PetIBM API](doc/markdowns/usepetibmapi.md)
* [Online API manual](https://barbagroup.github.io/PetIBM)
* [Change Log](CHANGELOG.md)
* [Contributing](CONTRIBUTING.md)

Offline API manual can be generated with [Doxygen](http://www.stack.nl/~dimitri/doxygen/).

---

## Papers published using PetIBM

* Mesnard, O., & Barba, L. A. (2017). _Reproducible and Replicable Computational Fluid Dynamics: It's Harder Than You Think_. Computing in Science & Engineering, 19(4), 44-55, https://doi.org/10.1109/MCSE.2017.3151254.


---

## Contact

Please e-mail [Olivier Mesnard](mailto:mesnardo@gwu.edu) or [Pi-Yueh Chuang](mailto:pychuang@gwu.edu) if you have any questions, suggestions, or feedback.

To report bugs, please use the GitHub issue tracking system.
We are also open to pull-requests.

---

## References

* Chang, W., Giraldo, F., & Perot, B. (2002). *Analysis of an exact fractional step method*. Journal of Computational Physics, 180(1), 183-199.
* Li, R. Y., Xie, C. M., Huang, W. X., & Xu, C. X. (2016). *An efficient immersed boundary projection method for flow over complex/moving boundaries*. Computers & Fluids, 140, 122-135.
* Perot, J. B. (1993). *An analysis of the fractional step method*. Journal of Computational Physics, 108(1), 51-58.
* Taira, K., & Colonius, T. (2007). *The immersed boundary method: a projection approach*. Journal of Computational Physics, 225(2), 2118-2137.

---

## How to cite PetIBM

If PetIBM contributes to a project that leads to a scientific publication, please cite the project.
You can use this citation or the BibTeX entry below.

### PetIBM - toolbox and applications of the immersed-boundary method on distributed-memory architectures

> Pi-Yueh Chuang, Olivier Mesnard, Anush Krishnan, Lorena A. Barba (2018). PetIBM: toolbox and applications of the immersed-boundary method on distributed-memory architectures. _Journal of Open Source Software_, **3**(25), 558, [doi:10.21105/joss.00558](https://doi.org/10.21105/joss.00558) 

```console
@article{chuang2018petibm,
  doi = {10.21105/joss.00558},
  url = {https://doi.org/10.21105/joss.00558},
  year = {2018},
  month = {may},
  publisher = {The Open Journal},
  volume = {3},
  number = {25},
  pages = {558},
  author = {Pi-Yueh Chuang and Olivier Mesnard and Anush Krishnan and Lorena A. Barba},
  title = {{PetIBM}: toolbox and applications of the immersed-boundary method on distributed-memory architectures},
  journal = {The Journal of Open Source Software}
}
```
