The file `simulationParameters.yaml` is required to provide information to the different solvers.

Here, we provide an example of such a file followed by an explanation of its content.

Note: text is case-sensitive.


## Example

The simulation parameters do not depend on the dimensionality of the problem, and so the same template can be used for both 2d and 3d simulations:

    - type: simulation
      dt: 0.01
      startStep: 0
      nt: 300
      nsave: 50
      restartFromSolution: false
      timeScheme: [ADAMS_BASHFORTH_2, CRANK_NICOLSON]
      ibmScheme: TAIRA_COLONIUS


## File options

* **dt**: (mandatory) the time-increment. Can be any number greater than zero (as long as it satisfies the stability criterion for the numerical scheme used).
* **startStep**: (optional, default: `0`) starting time-step. If different than the default value `0.0`, then the program reads the numerical solution of the stating time-step given. This parameters should be used to restart a simulation.
* **nt**: (mandatory) number of time-steps to execute.
* **nsave**: (mandatory) time-step interval at which Eulerian and Lagrangian quantities are saved into files.
* **restartFromSolution**: (optional) Boolean flag to tell PetIBM to start the simulation from time-step `0` with initial conditions for the velocity field different than those given on the `flowDescription.yaml` file. The initial velocity field is read from the file `q.dat` located in the sub-folder `0000000` of the simulation directory.
* **timeScheme**: (optional, default: `[EULER_EXPLICIT, EULER_IMPLICIT]`) specifies the time-stepping schemes for the convection and the diffusion terms, respectively. The convection term can take the options `EULER_EXPLICIT` and `ADAMS_BASHFORTH_2`, and the diffusion terms can take the options `EULER_EXPLICIT`, `EULER_IMPLICIT` and `CRANK_NICOLSON`.
* **ibmScheme**: (optional, default: `NAVIER_STOKES`) specifies the immersed boundary method used in the simulation. Currently, the only immersed boundary method available in PetIBM is `TAIRA_COLONIUS`. When no immersed body is present in the flow, pass the option `NAVIER_STOKES`.
