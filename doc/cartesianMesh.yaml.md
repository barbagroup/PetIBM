## 3-D Example

    - direction: x
      start: -15.0
      subDomains:
        - end: -0.6
          cells: 69
          stretchRatio: 0.952380952
        - end: 0.6
          cells: 48
          stretchRatio: 1.0
        - end: 15.0
          cells: 69
          stretchRatio: 1.05

    - direction: y
      start: -15.0
      subDomains:
        - end: -0.6
          cells: 69
          stretchRatio: 0.952380952
        - end: 0.6
          cells: 48
          stretchRatio: 1.0
        - end: 15.0
          cells: 69
          stretchRatio: 1.05

    - direction: z
      start: -1.0
      subDomains:
        - end: -1.0
          cells: 80
          stretchRatio: 1.0

## File Options

* **direction**: `x`, `y` or `z` 

The nodes along each direction of the mesh are specified separately. Each edge of the domain is split into smaller subdomains, and the cells within these subdomains either constitute a uniform grid or an exponential grid. For each of these subdomains, we specify the number of cells and the stretching ratio between consecutive cells. For a uniform section, the stretching ratio is 1.

Consider the edge along the x-axis in the above example file:

* **start**: This is the coordinate of the starting node of the mesh, which is the lowest coordinate along an edge in a particular direction.
* **subDomains**: We now specify the properties of each section of the edge. In the above example we have three such sections. The options for each are:
  * **end**: The start node of the first section is obtained from **start**, and the end node is specified by **end**. For subsequent sections, the **end** node of the previous section is the start node for the current one.
  * **cells**: Number of cells that the current section must be divided into.
  * **stretchRatio**: The ratio between the widths of consecutive cells. The value is greater than one in a region where the cell widths increase along the positive direction, less than one if the cell widths decrease along the positive direction, and equal to 1 for a uniform grid.  

The program automatically computes the widths of the cells when the above information is provided. Care must be taken by the user to ensure that the cells at the interface between two subDomains are approximately of equal size.

To illustrate: In the above grid, the extent of the domain along the x-direction is from -15.0 to +15.0. The mesh is described in three parts. The first part is the region [-15.0, -0.6], which consists of 69 cells that decrease in size in geometric progression, with the ratio 0.952. The second region consists of a uniform grid stretching from -0.6 to +0.6, with 48 cells. And the third region is [0.6, 15.0], with 69 cells increasing in size with a stretching ratio of 1.05. The grid along the y-direction is the same. And the grid along the z-direction is uniform and stretches from -1 to 1, consisting of 80 cells of width 0.025.

For a 2-D flow problem, the portion with **direction**: `z` must be omitted.