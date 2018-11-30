"""
Create a flat plate of length 1.0 with aspect ratio 2.0 and a 80-degree
inclination.
The plate is discretized with spacing 0.04 in the x-y plane and with spacing
0.04 along the z-direction.
"""

import math
import pathlib
import numpy


# Flat-plate's parameters.
L = 1.0  # chord length
AR = 2.0  # aspect ratio
xc, yc, zc = 0.0, 0.0, 0.0  # center's coordinates
aoa = 80.0  # angle of inclination in degrees
ds = 0.04  # mesh spacing

simu_dir = pathlib.Path(__file__).absolute().parents[1]

# Generate coordinates of the flat plate.
n = math.ceil(L / ds)
s = numpy.linspace(xc - L / 2, xc + L / 2, num=n + 1)

x = xc + numpy.cos(numpy.radians(-aoa)) * s
y = yc + numpy.sin(numpy.radians(-aoa)) * s

nz = math.ceil(L * AR / ds)
z = numpy.linspace(zc - L * AR / 2, zc + L * AR / 2, num=nz + 1)

# Write coordinates into file.
filepath = simu_dir / 'flatplateAoA{}.body'.format(aoa)
with open(filepath, 'w') as outfile:
    outfile.write('{}\n'.format(x.size * z.size))
for zi in z:
    with open(filepath, 'ab') as outfile:
        numpy.savetxt(outfile, numpy.c_[x, y, zi * numpy.ones(x.size)])
