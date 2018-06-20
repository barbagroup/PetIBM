"""
Computes the observed order of convergence for the velocity components and the
pressure using the solution on 4 consistently refined grids.
"""

import os
import numpy
import h5py
import pprint


def read_fields_from_hdf5(filepath, gridpath, names=[]):
  fields = {}
  f = h5py.File(filepath, 'r')
  fg = h5py.File(gridpath, 'r')
  for name in names:
    x, y = fg[name]['x'][:], fg[name]['y'][:]
    values = f[name][:]
    fields[name] = {'values': values, 'grid': {'x': x, 'y': y}}
  return fields


def restrict_field(field, grid, atol=1.0E-12):
  def intersection(a, b, atol=atol):
    return numpy.any(numpy.abs(a - b[:, numpy.newaxis]) <= atol, axis=0)
  values, x, y = field['values'], field['grid']['x'], field['grid']['y']
  mask_x = intersection(x, grid['x'], atol=atol)
  mask_y = intersection(y, grid['y'], atol=atol)
  return {'grid': {'x': x[mask_x], 'y': y[mask_y]},
          'values': numpy.array([values[j][mask_x]
                                 for j in range(y.size) if mask_y[j]])}


script_dir = os.path.dirname(os.path.realpath(__file__))
root_dir = os.path.dirname(script_dir)

timestep = 500
ncells = [30, 90, 270, 810]
folders = [str(n) for n in ncells]
ratio = 3.0

fields = {}
field_names = ['u', 'v', 'p']
for folder in folders:
  directory = os.path.join(root_dir, folder)
  filepath = os.path.join(directory, 'solution', '{:0>7}.h5'.format(timestep))
  gridpath = os.path.join(directory, 'grid.h5')
  fields[folder] = read_fields_from_hdf5(filepath, gridpath, names=field_names)
  # Restrict fields onto coarse grid.
  for name in field_names:
    fields[folder][name] = restrict_field(fields[folder][name],
                                          fields[folders[0]][name]['grid'])


# Computes observed orders of convergence.
alpha = {'first': {}, 'last': {}}
# Using the first three grids.
coarse = fields[folders[0]]
medium = fields[folders[1]]
fine = fields[folders[2]]
for name in field_names:
  alpha['first'][name] = (numpy.log(numpy.linalg.norm(medium[name]['values'] -
                                                      coarse[name]['values'],
                                                      ord=None) /
                                    numpy.linalg.norm(fine[name]['values'] -
                                                      medium[name]['values'],
                                                      ord=None)) /
                          numpy.log(ratio))
# Using the last three grids.
coarse = fields[folders[1]]
medium = fields[folders[2]]
fine = fields[folders[3]]
for name in field_names:
  alpha['last'][name] = (numpy.log(numpy.linalg.norm(medium[name]['values'] -
                                                     coarse[name]['values'],
                                                     ord=None) /
                                   numpy.linalg.norm(fine[name]['values'] -
                                                     medium[name]['values'],
                                                     ord=None)) /
                         numpy.log(ratio))

pprint.pprint(alpha)
