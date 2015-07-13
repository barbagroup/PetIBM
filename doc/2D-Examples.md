## Lid-driven cavity at Re=100

The fluid is in a square domain spanning [0,1]x[0,1], surrounded by solid walls on all sides. The top wall moves tangentially with velocity 1.0, and the fluid has kinematic viscosity 0.01. Using the side of the square as the characteristic length and the speed of the top wall as the characteristic velocity, the flow has Reynolds number 100.

In the current example, the fluid region is discretized using a uniform 32x32 grid, the fluid starts from rest, and the simulation is advanced in time in steps of 0.02 for 1000 time steps. The flow reaches steady state at the end of the simulation.

To run the case in parallel, first switch to the PetIBM root directory and call the corresponding rule in the `makefile`:

    > cd $PETIBM_DIR
    > make cavityRe100Parallel

The code is run using 4 MPI processes. Information related to the simulation is printed on the screen, along with the time steps when the flow data is saved.

After the simulation is complete, post-process the results:

    > python scripts/python/plot.py -folder cases/2d/lidDrivenCavity/Re100

The x-velocity, y-velocity, vorticity and pressure are plotted and saved in the folder `cases/2d/lidDrivenCavity/Re100/output`. The following figures show the velocity contours at the end of the simulation:

<img src="https://cloud.githubusercontent.com/assets/6268735/3593211/3731e81e-0c80-11e4-8c60-03e1a39658d5.png" alt="x-velocity: Lid-driven cavity at Re=100" height="400px"/>

<img src="https://cloud.githubusercontent.com/assets/6268735/3593213/3bef0bac-0c80-11e4-9d1a-9bdce81f214c.png" alt="y-velocity: Lid-driven cavity at Re=100" height="400px"/>

## Flow over a circular cyinder at Re=40

A circular cylinder of diameter 1.0 is placed at the center of a domain spanning [-15,15]x[-15,15]. The initial velocity of the fluid in the domain is (1, 0). The boundary conditions on the velocity are as follows: Dirichlet conditions at the inlet and the top and bottom boundaries, with value equal to the freestream velocity (1,0), and a convective boundary condition at the outlet, with the fluid transported out of the boundary at freestream velocity.

Run the simulation:

    > cd $PETIBM_DIR
    > make cylinderRe40

The flow is calculated for a total time of 3.0, in steps of 0.01. 2 MPI processes are used. The data is saved every 50 time steps.

Postprocess the results:

    > python scripts/python/plot.py -folder cases/2d/cylinder/Re40 -ymin -2 -ymax 2 -xmin -1.5 -xmax 3

We specify the extent of the plotting region through the command line options (the default behavior is to plot the entire domain). The output images are stored in the folder `cases/2d/cylinder/Re40/output`. The following are the vorticity and pressure plots at time step 300:

<img src="https://cloud.githubusercontent.com/assets/6268735/3594444/d7226534-0c9d-11e4-9733-9d801110f736.png" alt="Vorticity: Flow over a cylinder at Re=40" height="400px"/>

<img src="https://cloud.githubusercontent.com/assets/6268735/3594445/da5c5a5c-0c9d-11e4-8642-78edd6c04113.png" alt="Pressure: Flow over a cylinder at Re=40" height="400px"/>