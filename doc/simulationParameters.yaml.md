The file `simulationParameters.yaml` is required to provide information to the different solvers.

Here, we provide an example of such a file followed by an explanation of its content.

Note: text is case-sensitive.


## Example

The simulation parameters do not depend on the dimensionality of the problem, and so the same template can be used for both 2d and 3d simulations:

    - dt: 0.01
      startStep: 0
      nt: 300
      nsave: 50
      nrestart: 300
      ibm: TAIRA_COLONIUS
      convection: ADAMS_BASHFORTH_2
      diffusion: CRANK_NICOLSON
      outputFormat: binary
      outputFlux: true
      outputVelocity: true
      vSolveType: CPU
      pSolveType: CPU


## File options

* `dt`: (mandatory) the time-increment. Can be any number greater than zero (as long as it satisfies the stability criterion for the numerical scheme used).
* `startStep`: (optional, default: `0`) starting time-step. If different than the default value `0`, then the program reads the numerical solution of the stating time-step given. This parameters should be used to restart a simulation.
* `nt`: (mandatory) number of time-steps to execute.
* `nsave`: (mandatory) time-step interval at which Eulerian and Lagrangian quantities are saved into files.
* `nrestart`: (optional, default: `nt`) time-step interval at which the convective terms from the previous time-step are saved. This is useful to restart properly when using Adams-Bashforth second-order time-scheme for the convective terms.
* `ibm`: (optional) specifies the immersed boundary method used in the simulation. Currently, there are two immersed boundary methods implemented in PetIBM: `TAIRA_COLONIUS` and `LI_ET_AL`. `TAIRA_COLONIUS` is an immersed-boundary projection method where the pressure field and the Lagrangian forces are coupled together and a modified Poisson system is solved at each time step. `LI_ET_AL` is a decoupled version of the immersed-boundary projection method where the no-slip constraint and the divergence-free constraint are solved sequentially at each time step. If no immersed boundary are present in the computational domain, once should remove this line.
* `convection`: (optional, default: `EULER_EXPLICIT`) specifies the time-scheme to use for the convective terms of the momentum equation. In PetIBM, the convective terms can be temporally discretized using an explicit Euler method (`EULER_EXPLICIT`, default value) or a second-order Adams-Bashforth scheme (`ADAMS_BASHFORTH_2`).
* `diffusion`: (optional, default: `EULER_IMPLICIT`) specifies the time-scheme to use for the diffusive terms of the momentum equation. In PetIBM, the diffusive terms can be  treated explicitly (`EULER_EXPLICIT`), implicitly (`EULER_IMPLICIT`, default), or using a second-order Crank-Nicolson scheme (`CRANK_NICOLSON`).
* `outputFormat`: (optional, default: `binary`) specifies the format of the output files in which the numerical solution is stored. Right now, two formats are supported: `binary` and `hdf5`.
* `outputFlux`: (optional, default: `true`) writes the flux variable into files when set to `true`.
* `outputVelocity`: (optional, default: `false`) writes the velocity variable into files when set to `true`.
* `vSolveType`: (optional, default: `CPU`) to define which hardware to use to solve the iterative velocity system. Note: the `GPU` option is available only if PetIBM has been built with AmgXWrapper (see installation instructions for details about how to configure and build PetIBM with AmgXWrapper).
* `pSolveType`: (optional, default: `CPU`) to define which hardware to use to solve the iterative Poisson system. Note: the `GPU` option is available only if PetIBM has been built with AmgXWrapper (see installation instructions for details about how to configure and build PetIBM with AmgXWrapper).
* `decoupling`: (optional) information about the decoupling procedure when using the decoupled immersed-boundary projection method (`LI_ET_AL`). The YAML node contains the following parameters:
      - `algorithm`: (optional, default: `1`) index of the algorithm to apply. `1` satisfies the no-slip constraint first, then the divergence-free one. `3` satisfies the divergence-free constraint, then the no-slip one.
      - `forceEstimator`: (optional, default: `2`) index of the scheme to use to estimate the momentum forcing at the beginning of the time step. `1` sets the forcing to zero; `2` uses the forcing from the previous time step; `3` solves a system for the Lagrangian forces where the right-hand side is computed from the velocity fluxes at the previous time step.
      - `maxIters`: (optional, default: `1`) maximum number of iterations for the sub-iterative process of applying the constraints.
      - `atol`: (optional, default: `1.0E-05`) absolute tolerance criterion to stop the sub-iterative process. When the L2-norm of the Lagrangian forces variation vector is smaller than the provided absolute tolerance, the iterative process stops.
      - `rtol`: (optional, default: `1.0E-05`) relative tolerance criterion to stop the sub-iterative process. When the L2-norm of the Lagrangian forces variation vector is smaller than the provided relative tolerance times the L2-norm of the total Lagrangian forces vector, the iterative process stops.
      - `printStats`: (optional, default: false) when `true` (and when `maxIters` is greater than `1`), prints information of the sub-iterative process.
