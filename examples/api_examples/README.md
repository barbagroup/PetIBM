# API examples

This folder contains examples on how to use the PetIBM API to build flow solvers.
The examples provided here are just for demonstration purpose.
While they may still be used to conduct flow simulations, their functionalities are limited.

Solvers with more complete features can be found in the `applications` folder under the PetIBM directory.
(Binary executables for those solvers can be found in the `bin` folder of the PetIBM installation directory.)

Contents of the present folder:
* `liddrivencavity2d`: a basic Navier-Stokes solver using a projection method (Perot, 1993);
* `oscillatingcylinder2dRe100_GPU`: a Navier-Stokes solver using a decoupled immersed-boundary projection method (Li et al., 2016) to compute the 2D flow around an inline-oscillating cylinder (Poisson system solved on GPU devices).


__References:__

* Li, R. Y., Xie, C. M., Huang, W. X., & Xu, C. X. (2016). An efficient immersed boundary projection method for flow over complex/moving boundaries. Computers & Fluids, 140, 122-135.
* Perot, J. B. (1993). An analysis of the fractional step method. Journal of Computational Physics, 108(1), 51-58.