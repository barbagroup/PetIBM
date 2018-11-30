"""
Plots the instantaneous drag coefficient between 0 and 3 time-units of flow
simulation and compares with numerical results from
Koumoutsakos and Leonard (1995).

_References:_
* Koumoutsakos, P., & Leonard, A. (1995).
  High-resolution simulations of the flow around an impulsively started
  cylinder using vortex methods.
  Journal of Fluid Mechanics, 296, 1-38.
"""

import os
import pathlib
import numpy
import collections
from matplotlib import pyplot


simu_dir = pathlib.Path(__file__).absolute().parents[1]
root_dir = os.environ.get('PETIBM_EXAMPLES')
if not root_dir:
    root_dir = simu_dir.parents[1]

data = collections.OrderedDict({})

# Reads forces from file.
label = 'PetIBM'
filepath = simu_dir / 'forces-0.txt'
with open(filepath, 'r') as infile:
    t, fx = numpy.loadtxt(infile, dtype=numpy.float64,
                          unpack=True, usecols=(0, 1))
data[label] = {'t': t, 'cd': 2 * fx}
data[label]['kwargs'] = {}

# Reads drag coefficient of Koumoutsakos and Leonard (1995) for Re=3000.
label = 'Koumoutsakos and Leonard (1995)'
filename = 'koumoutsakos_leonard_1995_cylinder_dragCoefficientRe3000.dat'
filepath = root_dir / 'data' / filename
with open(filepath, 'r') as infile:
    t, cd = numpy.loadtxt(infile, dtype=numpy.float64, unpack=True)
data[label] = {'t': 0.5 * t, 'cd': cd}
data[label]['kwargs'] = {'linewidth': 0, 'marker': 'o',
                         'markerfacecolor': 'none', 'markeredgecolor': 'black'}

pyplot.rc('font', family='serif', size=16)

# Plots the instantaneous drag coefficients.
fig, ax = pyplot.subplots(figsize=(8.0, 6.0))
ax.grid()
ax.set_xlabel('Non-dimensional time')
ax.set_ylabel('Drag coefficient')
for label, subdata in data.items():
    ax.plot(subdata['t'], subdata['cd'], label=label, **subdata['kwargs'])
ax.axis((0.0, 3.0, 0.0, 2.0))
ax.legend()

pyplot.show()

# Save figure.
fig_dir = simu_dir / 'figures'
fig_dir.mkdir(parents=True, exist_ok=True)
filepath = fig_dir / 'dragCoefficient.png'
fig.savefig(str(filepath), dpi=300)
