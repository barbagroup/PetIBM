"""
Computes, plots, and saves the 2D vorticity field from a PetIBM simulation
after 3000 time steps (3 non-dimensional time-units).
"""

import os
import h5py
import numpy
from matplotlib import pyplot


script_dir = os.path.dirname(os.path.realpath(__file__))
simu_dir = os.path.dirname(script_dir)

# Read vorticity field from file.
name = 'wz'
filepath = os.path.join(simu_dir, 'grid.h5')
f = h5py.File(filepath, 'r')
x, y = f[name]['x'][:], f[name]['y'][:]
X, Y = numpy.meshgrid(x, y)
timestep = 3000
filepath = os.path.join(simu_dir, 'solution', '{:0>7}.h5'.format(timestep))
f = h5py.File(filepath, 'r')
wz = f[name][:]

# Read boundary coordinates from file.
filepath = os.path.join(simu_dir, 'circle.body')
with open(filepath, 'r') as infile:
  xb, yb = numpy.loadtxt(infile, dtype=numpy.float64, unpack=True, skiprows=1)

# Plot the filled contour of the vorticity.
pyplot.style.use('seaborn-dark')
fig, ax = pyplot.subplots(figsize=(6.0, 6.0))
ax.grid()
ax.set_xlabel('x', fontname='DejaVu Serif', fontsize=16)
ax.set_ylabel('y', fontname='DejaVu Serif', fontsize=16)
levels = numpy.linspace(-56.0, 56.0, 28)
ax.contour(X, Y, wz,
           levels=levels, colors='black')
ax.plot(xb, yb, color='red')
ax.set_xlim(-0.6, 1.6)
ax.set_ylim(-0.8, 0.8)
ax.set_aspect('equal')
for method in ['get_xticklabels', 'get_yticklabels']:
  for label in getattr(ax, method)():
    label.set_fontname('DejaVu Serif')
    label.set_fontsize(14)
fig.tight_layout()
pyplot.show()

# Saves figure.
figures_dir = os.path.join(simu_dir, 'figures')
if not os.path.isdir(figures_dir):
  os.makedirs(figures_dir)
filepath = os.path.join(figures_dir, 'vorticity{:0>7}.png'.format(timestep))
fig.savefig(filepath)
