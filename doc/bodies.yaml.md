The file `bodies.yaml` is required to provide information about the bodies immersed in the computational domain.

Note: text is case-sensitive.

There are different ways to specify bodies in the flow. Currently, PetIBM handles only one body per simulation.

The body description begins with the **type** of body, and the subsequent options are read accordingly. The following are the different types of bodies that can be input to PetIBM.


## 2d bodies

### Circle
A circular cylinder in the flow can be specified in the following way:

    - type: circle
      circleOptions: [0.0, 0.0, 0.5, 63]

When the type of body is **circle**, **circleOptions** with four parameters needs to be provided. The first two represent the coordinates of the center of the circle, the third represents the radius, and the fourth represents the number of points on the boundary (these will be uniformly distributed).

### Line Segment

    - type: lineSegment
      segmentOptions: [0.0, 0.0, 1.0, 1.0, 20]

A line segment is drawn when the **type** is `lineSegment`. **segmentOptions** takes five arguments: the first two are the coordinates of the start point of the segment, the next two specify the end point, and the final options gives the number of points the segment will be divided into. The segment is divided that many equal portions, and the points are created at the midpoints of the sections.

### Providing the coordinates

    - type: points
      pointsFile: bodyCoordinates.body

When the user chooses the **type** `points`, they also need to provide a file that contains a list of points that describe the immersed boundary. In the above example, this is the file `bodyCoordinates.body`. This file must be present in the case folder of the simulation. The first line of the file must specify the number of points in the file, and the remaining lines must list the coordinates of the points (x-coordinates in the first column, y-coordinates in the second). An example for such a file is:

    4
    0.0    0.0
    1.0    0.0
    1.0    1.0
    0.0    1.0

The solution behaves well and converges quicker when the distance between successive points on the boundary is approximately equal to the size of them mesh near the body.


## 3d bodies

### Providing the ooordinates

    - type: points
      pointsFile: bodyCoordinates.body

This is the same as for 2d bodies, but here, three coordinates need to be provided for every point (the third column represents the z-coordinate of each points). The first line of the file gives the number of points listed in the file. An example for such a file is:

    8
    0.0    0.0    0.0
    1.0    0.0    0.0
    1.0    1.0    0.0
    0.0    1.0    0.0
    0.0    0.0    1.0
    1.0    0.0    1.0
    1.0    1.0    1.0
    0.0    1.0    1.0
