Examples of the file are first provided, followed by an explanation of the contents. All text is case-sensitive.

## 2D Example
The following file could be used to simulate external flows over bodies:

    - dimensions: 2
      nu: 0.025
      initialVelocity: [1.0, 0.0]
      initialPerturbation: [0.0, 0.0]
      boundaryConditions:
        - location: xMinus
          u: [DIRICHLET, 1.0]
          v: [DIRICHLET, 0.0]
        - location: xPlus
          u: [CONVECTIVE, 1.0]
          v: [CONVECTIVE, 0.0]
        - location: yMinus
          u: [DIRICHLET, 1.0]
          v: [DIRICHLET, 0.0]
        - location: yPlus
          u: [DIRICHLET, 1.0]
          v: [DIRICHLET, 0.0]

## File options
* **dimensions**: Use `2` or `3` to specify whether the flow is 2-D or 3-D, respectively.
* **nu**: The kinematic viscosity of the fluid. Can be any positive number.
* **initialVelocity**: This is the initial velocity of the fluid in the flow field. This value is set throughout the domain. For 2-D flows, two components are required within square brackets, separated by a comma.
* **initialPerturbation**: Adds a perturbation field to the initial flow field which is proportional to the values specified in this option. This is useful when an instability in the flow needs to be triggered. For a uniform initial flow field, either do not include this option in the file, or set all values to zero.
* **boundaryConditions**: Specify the boundary conditions on each edge of the domain. These are listed in subsections, which need to be indented using two spaces.
  * Each boundary edge is represented by its relative **location** on the cartesian axes - `xMinus`, `xPlus`, `yMinus` and `yPlus`.
  * For each component of velocity (**u**, **v**) on the boundary, you can specify the type of boundary condition and an associated value. There are four types currently available: `DIRICHLET`, `NEUMANN`, `CONVECTIVE` and `PERIODIC`.
  * The numbers provided with the boundary conditions are interpreted as follows:
    * `DIRICHLET`: The value of that component of velocity at the boundary.
    * `NEUMANN`: The value of the normal derivative of that component at the boundary.
    * `CONVECTIVE`: The speed at which the fluid is convected out of the boundary. The values for the components normal to the boundary are ignored.
    * `PERIODIC`: All values are ignored.
  * Specific types of boundary conditions can be created by using a suitable combination of the above. For example, a slip boundary condition on the `yMinus` boundary would use homogeneous Neumann for **u** and zero Dirichlet for **v**.
* In the above example, `xMinus` is the inlet, `yMinus` and `yPlus` are the outer boundaries, and the velocity is fixed at (1,0) on all of them. `xPlus` is the outlet and the fluid is convected out with velocity (1,0). The initial velocity of the field is (1,0) with no perturbations, and the kinematic viscosity of the fluid is 0.025.

## 3D Example

    - dimensions: 3
      nu: 0.025
      initialVelocity: [1.0, 0.0, 0.0]
      initialPerturbation: [0.0, 0.0, 0.0]
      boundaryConditions:
        - location: xMinus
          u: [DIRICHLET, 1.0]
          v: [DIRICHLET, 0.0]
          w: [DIRICHLET, 0.0]
        - location: xPlus
          u: [CONVECTIVE, 1.0]
          v: [CONVECTIVE, 0.0]
          w: [CONVECTIVE, 0.0]
        - location: yMinus
          u: [DIRICHLET, 1.0]
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

The differences in a 3-D input file are the following:
* The number of **dimensions** is `3`.
* Three components need to be specified for **initialVelocity** and **initialPerturbation**.
* Two additional boundary locations: `zMinus` and `zPlus`.
* Boundary conditions are also required for the component of velocity **w**.