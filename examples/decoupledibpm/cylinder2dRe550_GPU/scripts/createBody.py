"""
Create a circle.
"""

import os
import numpy


# Circle's parameters.
R = 0.5  # radius
xc, yc = 0.0, 0.0  # center's coordinates
ds = 0.01  # distance between two consecutive points

script_dir = os.path.dirname(os.path.realpath(__file__))
simu_dir = os.sep.join(script_dir.split(os.sep)[:-1])

# Generate coordinates of the circle.
n = int(2.0 * numpy.pi * R / ds)
theta = numpy.linspace(0.0, 2.0 * numpy.pi, num=n, endpoint=False)
x, y = xc + R * numpy.cos(theta), yc + R * numpy.sin(theta)

# Write coordinates into file.
filepath = os.path.join(simu_dir, 'circle.body')
with open(filepath, 'w') as outfile:
  outfile.write('{}\n'.format(n))
with open(filepath, 'ab') as outfile:
  numpy.savetxt(outfile, numpy.c_[x, y])
