The file `bodies.yaml` is required to provide information about the boundaries immersed in the computational domain.

Note: text is case-sensitive.

Warning: currently, PetIBM can only handle a single stationary immersed body.

### Providing the coordinates

Here is an example of what the file `bodies.yaml` should look like:

    - type: points
      pointsFile: boundaryLocations.body

When the user chooses the `type` `points`, they also need to provide the path (relative to the simulation directory) of the file containing the body points.
In the example above, the file `boundaryLocations.body` contains the body points and should be located in the simulation directory.
The first line of the file must specify the number of points in the file, and the remaining lines must list the coordinates of the points (x-coordinates in the first column, y-coordinates in the second, and z-coordinates in the third one for a 3D simulation). An example for such a file is:

    4
    0.0    0.0
    1.0    0.0
    1.0    1.0
    0.0    1.0

The solution behaves well and converges quicker when the distance between successive points on the boundary is approximately equal to the size of the Eulerian mesh  in the vicinity of the boundary.
