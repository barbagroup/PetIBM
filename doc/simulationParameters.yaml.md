The file `simulationParameters.yaml` is required to provide information to the different solvers.

Here, we provide an example of such a file followed by an explanation of its content.

Note: text is case-sensitive.


## Example

The simulation parameters do not depend on the dimensionality of the problem, and so the same template can be used for both 2d and 3d simulations:

    - dt: 0.01
      startStep: 0
      nt: 300
      nsave: 50
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
* `ibm`: (optional) specifies the immersed boundary method used in the simulation. Currently, the only immersed boundary method available in PetIBM is `TAIRA_COLONIUS`. If no immersed boundary are present in the computational domain, once should remove this line.
* `convection`: (optional, default: `EULER_EXPLICIT`) specifies the time-scheme to use for the convective terms of the momentum equation. In PetIBM, the convective terms can be temporally discretized using an explicit Euler method (`EULER_EXPLICIT`, default value) or a second-order Adams-Bashforth scheme (`ADAMS_BASHFORTH_2`).
* `diffusion`: (optional, default: `EULER_IMPLICIT`) specifies the time-scheme to use for the diffusive terms of the momentum equation. In PetIBM, the diffusive terms can be  treated explicitly (`EULER_EXPLICIT`), implicitly (`EULER_IMPLICIT`, default), or using a second-order Crank-Nicolson scheme (`CRANK_NICOLSON`).
* `outputFormat`: (optional, default: `binary`) specifies the format of the output files in which the numerical solution is stored. Right now, two formats are supported: `binary` and `hdf5`.
* `outputFlux`: (optional, default: `true`) writes the flux variable into files when set to `true`.
* `outputVelocity`: (optional, default: `false`) writes the velocity variable into files when set to `true`.
* `vSolveType`: (optional, default: `CPU`) to define which hardware to use to solve the iterative velocity system. Note: the `GPU` option is available only if PetIBM has been built with AmgXWrapper (see installation instructions for details about how to configure and build PetIBM with AmgXWrapper).
* `pSolveType`: (optional, default: `CPU`) to define which hardware to use to solve the iterative Poisson system. Note: the `GPU` option is available only if PetIBM has been built with AmgXWrapper (see installation instructions for details about how to configure and build PetIBM with AmgXWrapper).
