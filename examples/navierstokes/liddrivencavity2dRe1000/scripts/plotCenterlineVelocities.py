"""
Plots the velocities along the centerlines of the 2D cavity at Reynolds number
1000 and compares with the numerical data reported in Ghia et al. (1982).

_References:_
* Ghia, U. K. N. G., Ghia, K. N., & Shin, C. T. (1982).
  High-Re solutions for incompressible flow using the Navier-Stokes equations
  and a multigrid method.
  Journal of computational physics, 48(3), 387-411.
"""

import os
import pathlib
import numpy
import h5py
from matplotlib import pyplot


# User's parameters
Re = 1000.0  # Reynolds number
time_step = 10000  # Time step at which to read the solution
# End of user's parameters

simu_dir = pathlib.Path(__file__).absolute().parents[1]
root_dir = os.environ.get('PETIBM_EXAMPLES')
if not root_dir:
    root_dir = simu_dir.parents[1]


def get_gridline_velocity(x_target, u, x, axis=0):
    i = numpy.where(x < x_target)[0][-1]
    x_a, x_b = x[i], x[i + 1]
    if axis == 0:
        u_a, u_b = u[:, i], u[:, i + 1]
    elif axis == 1:
        u_a, u_b = u[i], u[i + 1]
    return (u_a * (x_b - x_target) + u_b * (x_target - x_a)) / (x_b - x_a)


def read_data_ghia_et_al_1982(filepath, Re):
    with open(filepath, 'r') as infile:
        data = numpy.loadtxt(infile, dtype=numpy.float64, unpack=True)
    re2col = {100.0: (1, 7), 1000.0: (2, 8), 3200.0: (3, 9), 5000.0: (4, 10),
              10000.0: (5, 11)}
    return {'vertical': {'y': data[0], 'u': data[re2col[Re][0]]},
            'horizontal': {'x': data[6], 'v': data[re2col[Re][1]]}}


def read_field_hdf5(name, fieldpath, gridpath):
    field = {}
    f = h5py.File(gridpath, 'r')
    field['x'], field['y'] = f[name]['x'][:], f[name]['y'][:]
    f = h5py.File(fieldpath, 'r')
    field['values'] = f[name][:]
    return field


# Reads data from Ghia et al. (1982).
data_dir = root_dir / 'data'
filepath = data_dir / 'ghia_et_al_1982_lid_driven_cavity.dat'
ghia = read_data_ghia_et_al_1982(filepath, Re)

# Reads gridlines and velocity fields.
gridpath = simu_dir / 'grid.h5'
filepath = simu_dir / 'solution' / '{:0>7}.h5'.format(time_step)
u = read_field_hdf5('u', filepath, gridpath)
v = read_field_hdf5('v', filepath, gridpath)

# Computes x-velocity along vertical gridline at mid-cavity.
x_target = 0.5
u['vertical'] = get_gridline_velocity(x_target, u['values'], u['x'], axis=0)
# Computes y-velocity along horizontal gridline at mid-cavity.
y_target = 0.5
v['horizontal'] = get_gridline_velocity(y_target, v['values'], v['y'], axis=1)

pyplot.rc('font', family='serif', size=16)

# Plots the centerline velocities.
simu_kwargs = {'label': 'PetIBM',
               'color': '#336699', 'linestyle': '-', 'linewidth': 3,
               'zorder': 10}
ghia_kwargs = {'label': 'Ghia et al. (1982)',
               'color': '#993333', 'linewidth': 0,
               'markeredgewidth': 2, 'markeredgecolor': '#993333',
               'markerfacecolor': 'none',
               'marker': 'o', 'markersize': 8,
               'zorder': 10}
fig, ax = pyplot.subplots(nrows=2, figsize=(8.0, 8.0))
fig.suptitle('Re = {}'.format(int(Re)))
ax[0].grid()
ax[0].set_xlabel('y')
ax[0].set_ylabel('u (x={})'.format(x_target))
ax[0].plot(u['y'], u['vertical'], **simu_kwargs)
ax[0].plot(ghia['vertical']['y'], ghia['vertical']['u'], **ghia_kwargs)
ax[0].axis((0.0, 1.0, -0.75, 1.25))
ax[0].legend(loc='upper left')
ax[1].grid()
ax[1].set_xlabel('x')
ax[1].set_ylabel('v (y={})'.format(y_target))
ax[1].plot(v['x'], v['horizontal'], **simu_kwargs)
ax[1].plot(ghia['horizontal']['x'], ghia['horizontal']['v'], **ghia_kwargs)
ax[1].axis((0.0, 1.0, -0.75, 1.25))
ax[1].legend(loc='upper left')

pyplot.show()

# Save figure.
fig_dir = simu_dir / 'figures'
fig_dir.mkdir(parents=True, exist_ok=True)
filepath = fig_dir / 'centerlineVelocities{:0>7}.png'.format(time_step)
fig.savefig(str(filepath), dpi=300)
