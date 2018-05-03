"""
Create a circle.
"""

import os
import math
import numpy


# Circle's parameters.
R = 0.5  # radius
xc, yc = 0.0, 0.0  # center's coordinates
ds = 8.0 / 512  # mesh spacing

script_dir = os.path.dirname(os.path.realpath(__file__))
simu_dir = os.path.dirname(script_dir)

# Generate coordinates of the circle.
n = math.ceil(2.0 * numpy.pi * R / ds)
theta = numpy.linspace(0.0, 2.0 * numpy.pi, num=n)[:-1]
x, y = xc + R * numpy.cos(theta), yc + R * numpy.sin(theta)

# Write coordinates into file.
filepath = os.path.join(simu_dir, 'circle.body')
with open(filepath, 'w') as outfile:
  outfile.write('{}\n'.format(x.size))
with open(filepath, 'ab') as outfile:
  numpy.savetxt(outfile, numpy.c_[x, y])
