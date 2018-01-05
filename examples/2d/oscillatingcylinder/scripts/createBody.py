"""
Create a circle.
"""

import os
import math
import numpy


R = 0.5
xc, yc = 0.0, 0.0
ds = 8.0 / 512
n = math.ceil(2.0 * numpy.pi * R / ds)

script_dir = os.path.dirname(os.path.realpath(__file__))
simu_dir = os.sep.join(script_dir.split(os.sep)[:-1])

theta = numpy.linspace(0.0, 2.0 * numpy.pi, num=n)[:-1]
x, y = xc + R * numpy.cos(theta), yc + R * numpy.sin(theta)

filepath = os.path.join(simu_dir, 'circle.body')
with open(filepath, 'w') as outfile:
  outfile.write('{}\n'.format(x.size))
with open(filepath, 'ab') as outfile:
  numpy.savetxt(outfile, numpy.c_[x, y])
