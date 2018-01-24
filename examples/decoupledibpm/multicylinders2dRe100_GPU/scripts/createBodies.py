"""
Create two circles and write coordinates into files.
"""

import os
import numpy


def create_circle(R=0.5, center=(0.0, 0.0), ds=0.01):
  xc, yc = center
  n = int(2.0 * numpy.pi * R / ds)
  theta = numpy.linspace(0.0, 2.0 * numpy.pi, num=n, endpoint=False)
  x, y = xc + R * numpy.cos(theta), yc + R * numpy.sin(theta)
  return x, y


def write_coordinates(coords, filepath):
  x, y = coords
  n = x.size
  with open(filepath, 'w') as outfile:
    outfile.write('{}\n'.format(n))
  with open(filepath, 'ab') as outfile:
    numpy.savetxt(outfile, numpy.c_[x, y])


script_dir = os.path.dirname(os.path.realpath(__file__))
simu_dir = os.sep.join(script_dir.split(os.sep)[:-1])

circle1 = create_circle(R=0.5, center=(0.0, -2.5), ds=0.02)
filepath = os.path.join(simu_dir, 'circle1.body')
write_coordinates(circle1, filepath)

circle2 = create_circle(R=0.5, center=(0.0, 2.5), ds=0.02)
filepath = os.path.join(simu_dir, 'circle2.body')
write_coordinates(circle2, filepath)
