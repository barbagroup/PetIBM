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
import numpy
from matplotlib import pyplot


if not os.environ.get('PETIBM_EXAMPLES'):
  raise KeyError('Environment variable PETIBM_EXAMPLES is not set; '
                 'Set PETIBM_EXAMPLES as the root directory of the examples.')

script_dir = os.path.dirname(os.path.realpath(__file__))
simu_dir = os.sep.join(script_dir.split(os.sep)[:-1])
root_dir = os.environ['PETIBM_EXAMPLES']

# Reads forces from file.
filepath = os.path.join(simu_dir, 'forces.txt')
with open(filepath, 'r') as infile:
  data = numpy.loadtxt(infile, dtype=numpy.float64, unpack=True)
times, cd = data[0], 2.0 * data[1]

# Reads drag coefficient of Koumoutsakos and Leonard (1995) for Re=3000.
filename = 'koumoutsakos_leonard_1995_cylinder_dragCoefficientRe3000.dat'
filepath = os.path.join(root_dir, 'data', filename)
with open(filepath, 'r') as infile:
  times_kl, cd_kl = numpy.loadtxt(infile, dtype=numpy.float64, unpack=True)
times_kl *= 0.5

# Plots the instantaneous drag coefficients.
pyplot.style.use('seaborn-dark')
kwargs_data = {'label': 'PetIBM',
               'color': '#336699',
               'linestyle': '-',
               'linewidth': 3,
               'zorder': 10}
kwargs_kl1995 = {'label': 'Koumoutsakos and Leonard (1995)',
                 'color': '#993333',
                 'linewidth': 0,
                 'markeredgewidth': 2,
                 'markeredgecolor': '#993333',
                 'markerfacecolor': 'none',
                 'marker': 'o',
                 'markersize': 4,
                 'zorder': 10}
fig, ax = pyplot.subplots(figsize=(6, 6))
ax.grid(True, zorder=0)
ax.set_xlabel('non-dimensional time', fontsize=16)
ax.set_ylabel('drag coefficient', fontsize=16)
ax.plot(times, cd, **kwargs_data)
ax.plot(times_kl, cd_kl, **kwargs_kl1995)
ax.axis([0.0, 3.0, 0.0, 2.0])
ax.legend(prop={'size': 16})

pyplot.show()

# Saves figure.
figures_dir = os.path.join(simu_dir, 'figures')
if not os.path.isdir(figures_dir):
  os.makedirs(figures_dir)
filepath = os.path.join(figures_dir, 'dragCoefficient.png')
fig.savefig(filepath)
