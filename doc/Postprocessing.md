Three postprocessing scripts are available in the folder `${PETIBM_DIR}/scripts/python`: `plot.py`, `plot2d.py` and `plot3d.py`.

## `plot.py`

This script can be used to plot the x and y components of velocity, the pressure and the vorticity at every saved time step of a 2-D PetIBM simulation. It uses the `python` module `matplotlib` to generate the plots.

Provide the case folder when running the script:

	> ./scripts/python/plot.py -folder cases/2d/lidDrivenCavity/Re100

To specify the extents of the plot region, use the following options:

	> ./scripts/python/plot.py -folder cases/2d/cylinder/Re40 -xmin -2 -xmax 4 -ymin -3 -ymax 3

For a complete list of available command line options, run the script with the `-h` or `--help` flag:
    
	> ./scripts/python/plot.py -h

## `plot2d.py`

This script generates VTK files of the velocity field in a specified rectangular region of the flow. These files are generated for every save point, or can be generated for a subset of save points by using the options `-startStep` and `-nsave`.

	> ./scripts/python/plot2d.py -folder cases/2d/cylinder/Re40 -startStep 100 -nsave 100

On a staggered Cartesian grid, each velocity component is specified on a cell face perpendicular to it. In a VTK file, all three components of the velocity must be specified at the same point. What the above script does is average the components on opposite faces of a cell, and print the values at the locations of the cell centers.

## `plot3d.py`

This script generates VTK files of the velocity field in a specified cuboidal region of the flow using the same procedure detailed above. List all the available command line options using the `-h` or `--help` flags:

	> ./scripts/python/plot3d.py --help

## Visualizing the VTK files

The files output by `plot2d.py` or `plot3d.py` that contain the velocity fields can be visualized using [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/) or [Paraview](http://www.paraview.org/). These software also have options to view derived quantities such as vorticity or streamlines, and allow you to define custom flow quantities such as the Q-criterion.