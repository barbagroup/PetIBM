"""
Plots the instantaneous drag coefficient up to 4 periods of flow simulation
for the case of an inline oscillating cylinder.
"""

import os
import numpy
from matplotlib import pyplot


Um = 1.0  # Maximum translational velocity of the body
f = 0.2  # frequency of oscillation
T = 1.0 / f  # period
w = 2.0 * numpy.pi * f  # angular frequency
Am = Um / w  # amplitude of oscillation

script_dir = os.path.dirname(os.path.realpath(__file__))
simu_dir = os.sep.join(script_dir.split(os.sep)[:-1])

# Reads forces from file.
filepath = os.path.join(simu_dir, 'forces.txt')
with open(filepath, 'r') as infile:
  times, fx, _ = numpy.loadtxt(infile, dtype=numpy.float64, unpack=True)

# Computes acceleration in the x-direction of the body.
acc_x = Am * w**2 * numpy.sin(w * times)

# Computes the drag coefficient.
cd = 2.0 * (fx + acc_x)

# Plots the instantaneous drag coefficient.
pyplot.style.use('seaborn-dark')
kwargs_data = {'label': 'PetIBM',
               'color': '#336699',
               'linestyle': '-',
               'linewidth': 2,
               'zorder': 10}
fig, ax = pyplot.subplots(figsize=(8, 4))
ax.grid(True, zorder=0)
ax.set_xlabel('$t / T$', fontsize=16)
ax.set_ylabel('$C_D$', fontsize=16)
ax.plot(times / T, cd, **kwargs_data)
ax.axis([0.0, 4.0, -6.0, 8.0])
ax.legend(loc='upper right', prop={'size': 16})

pyplot.show()

# Saves figure.
figures_dir = os.path.join(simu_dir, 'figures')
if not os.path.isdir(figures_dir):
  os.makedirs(figures_dir)
filepath = os.path.join(figures_dir, 'dragCoefficient.png')
fig.savefig(filepath)
