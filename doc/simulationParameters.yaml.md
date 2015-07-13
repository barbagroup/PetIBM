## Example

The simulation parameters do not depend on the dimensionality of the flow, and so the same template can be used for both 2-D and 3-D flows:

    - type: simulation
      dt: 0.01
      nt: 300
      nsave: 50
      restart: false
      startStep: 0
      timeScheme: [ADAMS_BASHFORTH_2, CRANK_NICOLSON]
      ibmScheme: TAIRA_COLONIUS

## Options

* **dt**: The size of the time step. Can be any number greater than zero (as long as it satisfies the stability criterion for the numerical scheme used)
* **nt**: Number of time steps for which the simulation will run.
* **nsave**: The interval at which the flow data is saved. The velocity fluxes, pressure, and body forces in the entire domain are saved, whenever **nsave** divides the current time step number exactly.
* **restart**: A flag to tell the solver whether it needs to restart the simulation using previously saved data, or whether it should start it from the beginning. The simulation will restart from a previous save point if this option is set to `yes` or `true`, and will start from the beginning when the option is `no` or `false`.
  * When a simulation is restarted, any of the other simulation parameters can be changed. Only the initial conditions will be read from the specified save point.
* **startStep**: When **restart** is set to `false`, this option is ignored and the default value is `0`. If restart is `true`, then the data from the time step `startStep` is used as the initial data to restart the simulation.
* **timeScheme**: Specify the time-stepping schemes for the convection and the diffusion terms, respectively. The convection term can take the options `EULER_EXPLICIT` and `ADAMS_BASHFORTH_2`, and the diffusion terms can take the options `EULER_EXPLICIT`, `EULER_IMPLICIT` and `CRANK_NICOLSON`.
* **ibmScheme**: Specify the immersed boundary method used in the simulation. Currently, the only immersed boundary method available in PetIBM is `TAIRA_COLONIUS`. When no immersed body is present in the flow, pass the option `NAVIER_STOKES`.