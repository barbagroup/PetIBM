"""
Create a circle.
"""

import pathlib
import math
import numpy


root_dir = pathlib.Path(__file__).absolute().parents[1]

# Circle's parameters.
R = 0.5  # radius
xc, yc = 0.0, 0.0  # center's coordinates
ds = 0.025  # distance between two consecutive points

# Create coordinates of the circle.
n = math.ceil(2 * numpy.pi * R / ds)  # number of divisions
theta = numpy.linspace(0.0, 2.0 * numpy.pi, num=n + 1)[:-1]
x, y = xc + R * numpy.cos(theta), yc + R * numpy.sin(theta)

# Write coordinates into file.
filepath = root_dir / 'circle.body'
with open(filepath, 'w') as outfile:
    outfile.write('{}\n'.format(n))
with open(filepath, 'ab') as outfile:
    numpy.savetxt(outfile, numpy.c_[x, y])
