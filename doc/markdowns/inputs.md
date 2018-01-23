# Input files

PetIBM reads input parameters from a YAML configuration file (typically called `config.yaml`) that is located in the simulation directory (i.e. directory where the numerical solution will be written).

If the current working directory is not the simulation directory, the user should provide the command-line argument `-directory <simulation-directory>`.

If the configuration file is located in another directory, the user should provide the path of the file using the command-line argument `-config <file-path>`.
(The path can be absolute or relative to the current working directory.)

---

## Basic YAML configuration file

A basic YAML configuration file for a 2D run will look like:

```yaml
mesh:
  - direction: x
    start: 0.0
    subDomains:
      - end: 1.0
        cells: 32
        stretchRatio: 1.0
  - direction: y
    start: 0.0
    subDomains:
      - end: 1.0
        cells: 32
        stretchRatio: 1.0

flow:
    nu: 0.01
    initialVelocity: [0.0, 0.0]
    boundaryConditions:
      - location: xMinus
        u: [DIRICHLET, 0.0]
        v: [DIRICHLET, 0.0]
      - location: xPlus
        u: [DIRICHLET, 0.0]
        v: [DIRICHLET, 0.0]
      - location: yMinus
        u: [DIRICHLET, 0.0]
        v: [DIRICHLET, 0.0]
      - location: yPlus
        u: [DIRICHLET, 1.0]
        v: [DIRICHLET, 0.0]

parameters:
    dt: 0.01
    startStep: 0
    nt: 1000
    nsave: 1000
    nrestart: 1000
    convection: ADAMS_BASHFORTH_2
    diffusion: CRANK_NICOLSON
    velocitySolver:
      type: CPU
      config: solversPetscOptions.info
    poissonSolver:
      type: CPU
      config: solversPetscOptions.info

bodies:
  - type: points
    file: circle.body
```

The configuration contains three mandatory nodes: `mesh`, `flow`, `parameters`.
The node `bodies` is required when running a simulation with at least one boundary immersed in the computational domain.

---

## YAML node `mesh`

`mesh` controls the parameters of the structured Cartesian grid to be created.
We need to provide sufficient information to generate the gridline coordinates in the `direction` `x`, `y`, and `z`.
(For 2D runs, the `z` direction should be omitted.)

A gridline is defined by the location of a starting point (`start`), which corresponds to the minimum value along the gridline (for example, location of the left boundary in the x-direction).

A gridline is divided into consecutive sub-domains (`subDomains`), each one being defined with the coordinate of a ending point (`end`), a number of cells (`cells`), and optionally a stretching ratio (`stretchRatio`).
The stretching ratio represents the ratio between the width of two consecutive cells.
(A ratio greater than `1.0` means that the cell width increases along the gridline in the positive direction.)

Based on those information, PetIBM computes the cell widths and the user should make sure that the cells at the interface between two sub-domains are approximately of equal size.

Example:

The following node will generate a 2D structured Cartesian grid with 1704 cells in the x and y directions in the domain [-15, 15]x[-15, 15] where the mesh is uniform in the sub-domain [-0.52, 1.48]x[-2, 2] and stretched to the external boundaries with a constant ratio of 1.01.

```yaml
mesh:
  - direction: x
    start: -15.0
    subDomains:
      - cells: 363
        end: -0.52
        stretchRatio: 0.9900990099009901
      - cells: 1000
        end: 3.48
        stretchRatio: 1.0
      - cells: 341
        end: 15.0
        stretchRatio: 1.01
  - direction: y
    start: -15.0
    subDomains:
      - cells: 352
        end: -2.0
        stretchRatio: 0.9900990099009901
      - cells: 1000
        end: 2.0
        stretchRatio: 1.0
      - cells: 352
        end: 15.0
        stretchRatio: 1.01
```

---

## YAML node `flow`

`flow` prescribed the characteristics of the fluid as well as the initial and boundary conditions of the flow velocity.

`nu` is the kinematic viscosity of the fluid.

`initialVelocity` describes the initial uniform velocity vector field.

`boundaryConditions` lists the type and value of a velocity component (`u`, `v`, or `w`) for all boundaries.
The boundary locations are: `xMinus` and `xPlus` for the left and right, `yMinus` and `yPlus` for the bottom and top, and `zMinus` and `zPlus` for the front and back boundaries.
(`zMinus` and `zPLus` should be omitted for 2D runs.)
PetIBM implements different type of boundary conditions:

- Dirichlet (`DIRICHLET`); the value corresponds the value of the velocity component at the boundary;
- Neumman (`NEUMANN`); the value represents the value of the derivative normal to the boundary of the velocity component;
- convective (`CONVECTIVE`); the value is the speed at which the velocity component is convected;
- periodic (`PERIODIC`); the value has no meaning and can be set to whatever number.

The following node corresponds to flow characteristics for a 3D cavity flow (initially at rest) at Reynolds number 500 (based on the kinematic viscosity, the length of the cavity, and the speed of the lid-driven wall) with a lid-driven wall at the top boundary (moving with speed 1 in the x direction) and with periodic boundary conditions in the z direction.

```yaml
flow:
    nu: 0.002
    initialVelocity: [0.0, 0.0, 0.0]
    boundaryConditions:
      - location: xMinus
        u: [DIRICHLET, 0.0]
        v: [DIRICHLET, 0.0]
        w: [DIRICHLET, 0.0]
      - location: xPlus
        u: [DIRICHLET, 0.0]
        v: [DIRICHLET, 0.0]
        w: [DIRICHLET, 0.0]
      - location: yMinus
        u: [DIRICHLET, 0.0]
        v: [DIRICHLET, 0.0]
        w: [DIRICHLET, 0.0]
      - location: yPlus
        u: [DIRICHLET, 1.0]
        v: [DIRICHLET, 0.0]
        w: [DIRICHLET, 0.0]
      - location: zMinus
        u: [PERIODIC, 0.0]
        v: [PERIODIC, 0.0]
        w: [PERIODIC, 0.0]
      - location: zPlus
        u: [PERIODIC, 0.0]
        v: [PERIODIC, 0.0]
        w: [PERIODIC, 0.0]
```

---

## YAML node `parameters`

This node gather various `parameters` regarding the advancement of a simulation.

- `dt`: the time-step size.
- `startStep`: index of the stating time step (default is zero).
- `nt`: number of time steps to compute.
- `nsave`: frequency (in number of time steps) of saving for the numerical solution.
- `nrestart`: frequency (in number of time steps) of saving for the convective and diffusive terms; those terms will required upon restart of a run at a time step different from 0.
- `convection`: time scheme for the convective terms; choices are the default explicit Euler method (`EULER_EXPLICIT`) or an explicit second-order Adams-Bashforth scheme (`ADAMS_BASHFORTH_2`).
- `diffusion`: time scheme for the diffusive terms; choices are the default implicit Euler method (`EULER_IMPLICIT`), an explicit Euler method (`EULER_EXPLICIT`), or a second-order Crank-Nicolson scheme (`CRANK_NICOLSON`).
- `velocitySolver`, `poissonSolver`, and `forcesSolver` (for the decoupled version of the immersed-boundary projection method) each references the type of hardware used to solve the linear system (either `CPU` or `GPU`) and the path (relative to the YAML configuration file) of the file containing the parameters for the linear solver.

In the following example, PetIBM will run 1000 time steps (from time step 0) with a time increment of 0.01, saving the numerical solution (velocity vector field, pressure scalar field, and Lagrangian boundary forces) every 100 time steps and saving the convective and diffusive terms every 200 time steps.

```yaml
parameters:
    dt: 0.01
    startStep: 0
    nt: 1000
    nsave: 100
    nrestart: 200
    convection: ADAMS_BASHFORTH_2
    diffusion: CRANK_NICOLSON
    velocitySolver:
      type: CPU
      config: solversPetscOptions.info
    poissonSolver:
      type: GPU
      config: solversAmgXOptions.info
```

In the example above, the velocity  system will be solved on CPU with PETSc KSP and the parameters of the linear solver can be prescribed in the file `solversPetscOptions.info`.
(The PETSc prefix for the arguments is `-velocity_`, `-poisson_`, and `-forces_` for the velocity system, the Poisson system, and the Lagrangian forces system, respectively.)

The Poisson system will be solved on GPU devices with the [NVIDIA AmgX library](https://github.com/NVIDIA/AMGX) and the parameters of the linear solver are prescribed in the file `solversAmgXOptions.info`.

---

## YAML node `bodies`

`bodies` provides information about the boundaries immersed in the computational domain.

For example, the following node stipulates that there are two immersed boundaries defined by `points` coordinates that can be read in the ASCII files `boundary1.txt` and `boundary2.txt`, respectively.
(The path of the files can be either absolute or relative to the simulation directory.)

```yaml
bodies:
  - type: points
    file: boundary1.txt
  - type: points
    file: boundary2.txt
```

The file `boundary1.txt` contains the coordinates of the first boundary with the first line storing the number of points.
The other lines contain the coordinates of points describing the boundary (2 columns for 2D simulations, 3 columns for 3D simulations).

Example of a 2D boundary with 4 points:
```
4
0.0    0.0
1.0    0.0
1.0    1.0
0.0    1.0
```
