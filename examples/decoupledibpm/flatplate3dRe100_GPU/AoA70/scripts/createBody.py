"""
Create a flat plate of length 1.0 with aspect ratio 2.0 and a 70-degree
inclination.
The plate is discretize with spacing 0.04 in the x-y plane and with spacing
0.04 along the z-direction.
"""

import os
import numpy


L = 1.0
AR = 2.0
xc, yc, zc = 0.0, 0.0, 0.0
aoa = 70  # angle of inclination in degrees
ds = 0.04

script_dir = os.path.dirname(os.path.realpath(__file__))
simu_dir = os.sep.join(script_dir.split(os.sep)[:-1])

n = int(L / ds)
section = numpy.linspace(xc - L / 2, xc + L / 2, n)

x = xc + numpy.cos(numpy.radians(-aoa)) * section
y = yc + numpy.sin(numpy.radians(-aoa)) * section

nz = int(L * AR / ds)
z = numpy.linspace(zc - L * AR / 2, zc + L * AR / 2, nz)

filepath = os.path.join(simu_dir, 'flatplateAoA{}.body'.format(aoa))
with open(filepath, 'w') as outfile:
  outfile.write('{}\n'.format(n * nz))
for zi in z:
  with open(filepath, 'ab') as outfile:
    numpy.savetxt(outfile, numpy.c_[x, y, zi * numpy.ones(n)])
