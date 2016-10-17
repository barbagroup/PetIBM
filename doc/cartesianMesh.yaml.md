The file `cartesianMesh.yaml` is required to generate a Cartesian grid (uniform or stretched mesh).

Here, we provide an example of such a file followed by an explanation of its content.

Note: text is case-sensitive.


## 3d example

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


## File options

To generate a Cartesian grid, one needs to provide sufficient information to generate the coordinate of the coordinates along a grid-line in a given `direction` (could be `x` or `y` -- or `z` for a 3D problem).

A grid-line is defined by a starting-node (`start`, which is the lowest coordinate on the grid-line) and is divided into consecutive sub-domains (keyword `subDomains`), each one defined by: an ending-node (`end`), a number of cells (`cells`) and possibly a stretching-ratio (`stetchRatio`, optional with default value `1.0`). When specified, the stretching-ration corresponds to the ratio between the widths of two consecutive cells. A ratio greater than `1.0` means that the cell-width increases along the grid-line in the positive direction.

The program automatically computes the widths of the cells when the above information is provided. Care must be taken by the user to ensure that the cells at the interface between two sub-domains are approximately of equal size.

To illustrate: In the above grid, the extent of the domain along the x-direction is from `-15.0` to `+15.0`. The x-direction is divided into three parts. The first part is the region from `-15.0` to `-0.6`, which consists of `69` cells that decrease in size in geometric progression, with the ratio `0.952`. The second region consists of a uniform grid from `-0.6` to `+0.6`, with `48` cells. And the third region is `+0.6`, `+15.0`, with `69` cells increasing in size with a stretching ratio of `1.05`. The grid along the y-direction is the same. And the grid along the z-direction is uniform and stretches from -1 to 1, consisting of 80 cells of width 0.025.

For a 2D flow problem, the block `- direction: z` should be omitted.
