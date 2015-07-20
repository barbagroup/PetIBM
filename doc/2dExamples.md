## Lid-driven cavity flow at Re=100

The input files for this example are located in the directory `$PETIBM_DIR/examples/2d/lidDrivenCavity/Re100`.

The fluid, with kinematic viscosity `0.01` is enclosed in a two-dimensional box that spans `[0.0, 1.0]` in both directions. All solid walls are fixed except the top one that moves in the x-direction at velocity speed `1.0`. The Reynolds number of the flow, based on the box-length and the top-wall speed, is `100`. The fluid is initially at rest and the top-wall impulsively starts moving.

The computational domains discretized using a uniform `32x32` Cartesian grid. We choose a time-increment of `0.02`.

To run the case in parallel using `4` MPI processes, change directory to your PetIBM build directory, then:

    > cd examples
    > make examples
    > make cavity2dRe100Parallel

The program runs for `1000` time-steps, and saves the numerical solution at the final time-step (at which we consider the flow has reach a steady-state).

We may now post-process the results:

    > python $PETIBM_DIR/scripts/python plotFields2d.py --case 2d/lidDrivenCavity/Re100

The command-line above, generates the contours of the velocity components, the pressure field and the vorticity field after `1000` time-steps. The figures are saved as `.png` files in the sub-folder `images` of the simulation directory.

Here are the figures generated:

[[images/lidDrivenCavity2dRe100/pressure0001000.png]]

[[images/lidDrivenCavity2dRe100/uVelocity0001000.png]]

[[images/lidDrivenCavity2dRe100/vVelocity0001000.png]]

[[images/lidDrivenCavity2dRe100/vorticity0001000.png]]


## flow over a circular cylinder at Re=40

A circular cylinder of diameter `1.0` is placed at the center of a two-dimensional domain spanning `[-15,15]x[-15,15]`. The initial velocity of the fluid in the domain is `(1, 0)`. Dirichlet boundary conditions for the velocity are on all edges (velocity set to `(1.0, 0.0)`), except at the outlet where a the fluid is convective outside the domain in the x-direction at speed `1.0`. The Reynolds number, based on the freestream speed and the diameter of the circular cylinder, is `40`.

The computational domain is discretized using an stretched grid with `186x186` cells. The uniform region consists of the box `[-0.6, 0.6]x[-0.6, 0.6]` that contains `48x48` cells, giving a grid-resolution of `0.025` in both directions. Then, a stretching-ratio of `1.05` is applied in all directions from the edges of the uniform region to the domain limits.

We run the simulation for `300` time-steps with time-increment of `0.01`, saving the numerical solution every `50` iterations.

To run this example using `2` MPI processes, navigate to your PetIBM build directory, then:

    > cd examples
    > make examples
    > make cylinder2dRe40

You may now post-process the results:

    > python $PETIBM_DIR/scripts/python/plotFields2d.py --case 2d/cylinder/Re40 \
                  -- bottom-left -2.0 -2.0 \
                  --top-right 4.0 2.0 \
                  --no-velocity \
                  --pressure-range -1.0 0.5 101 \
                  --vorticity-range -1.0 1.0 101

The command-line above generates the contours within the box `[-2.0, 4.0]x[-2.0, 2.0]` of the pressure field, the vorticity field (between `-1.0` and `1.0` using `101` contours), excluding the velocity fields, at each time-step saved.

Here are the figures generated:

[[images/cylinder2dRe40/pressure0000300.png]]

[[images/cylinder2dRe40/vorticity0000300.png]]
