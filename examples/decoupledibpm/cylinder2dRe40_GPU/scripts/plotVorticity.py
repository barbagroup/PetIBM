"""
Computes, plots, and saves the 2D vorticity field from a PetIBM simulation
after 2000 time steps (20 non-dimensional time-units).
"""

import pathlib
import h5py
import numpy
from matplotlib import pyplot


simu_dir = pathlib.Path(__file__).absolute().parents[1]

# Read vorticity field and its grid from files.
name = 'wz'
filepath = simu_dir / 'grid.h5'
f = h5py.File(filepath, 'r')
x, y = f[name]['x'][:], f[name]['y'][:]
X, Y = numpy.meshgrid(x, y)
timestep = 2000
filepath = simu_dir / 'solution' / '{:0>7}.h5'.format(timestep)
f = h5py.File(filepath, 'r')
wz = f[name][:]

# Read body coordinates from file.
filepath = simu_dir / 'circle.body'
with open(filepath, 'r') as infile:
    xb, yb = numpy.loadtxt(infile, dtype=numpy.float64,
                           unpack=True, skiprows=1)

pyplot.rc('font', family='serif', size=16)

# Plot the filled contour of the vorticity.
fig, ax = pyplot.subplots(figsize=(6.0, 6.0))
ax.grid()
ax.set_xlabel('x')
ax.set_ylabel('y')
levels = numpy.linspace(-3.0, 3.0, 16)
ax.contour(X, Y, wz, levels=levels, colors='black')
ax.plot(xb, yb, color='red')
ax.set_xlim(-1.0, 4.0)
ax.set_ylim(-2.0, 2.0)
ax.set_aspect('equal')
fig.tight_layout()

pyplot.show()

# Save figure.
fig_dir = simu_dir / 'figures'
fig_dir.mkdir(parents=True, exist_ok=True)
filepath = fig_dir / 'wz{:0>7}.png'.format(timestep)
fig.savefig(str(filepath), dpi=300)
