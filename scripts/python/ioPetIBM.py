#!/usr/bin/env python

# file: ioPetIBM.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Collection of IO functions for PetIBM.


import os
import sys

import numpy
from matplotlib import pyplot
sys.path.append(os.path.join(os.environ['PETSC_DIR'], 'bin', 'pythonscripts'))
import PetscBinaryIO


def read_grid_old(case_directory):
  """Reads the coordinates from the file grid.txt
  and computes the cell-widths in each direction.

  Arguments
  ---------
  case_directory -- directory of the simulation
  """
  print('Reading the mesh ...')
  grid_path = '{}/grid.txt'.format(case_directory)
  with open(grid_path, 'r') as infile:
    nCells = numpy.array([int(n) for n in infile.readline().strip().split()])
    nodes = numpy.loadtxt(infile, dtype=float)
  coords = numpy.array(numpy.split(nodes, numpy.cumsum(nCells[:-1]+1)))
  cell_widths = coords[:, 1:] - coords[:, :-1]
  return coords, cell_widths

def read_grid(case_directory,
              bottom_left=[float('-inf'), float('-inf'), float('-inf')],
              top_right=[float('inf'), float('inf'), float('inf')]):
  """Reads the coordinates from the file grid.txt
  and computes the cell-widths in each direction.

  Arguments
  ---------
  case_directory -- directory of the simulation
  bottom_left -- bottom-left corner coordinates of the view 
                 (default -inf, -inf, -inf)
  top_right -- top-right corner coordinates of the view 
               (default inf, inf, inf)
  """
  print('Reading the mesh ...')
  grid_path = '{}/grid.txt'.format(case_directory)
  with open(grid_path, 'r') as infile:
    nCells = numpy.array([int(n) for n in infile.readline().strip().split()])
    nodes = numpy.loadtxt(infile, dtype=float)
  coords = numpy.array(numpy.split(nodes, numpy.cumsum(nCells[:-1]+1)))
  coords = numpy.ma.masked_array([numpy.ma.masked_outside(coords[i], 
                                                          bottom_left[i], 
                                                          top_right[i], 
                                                          copy=False) 
                                  for i in xrange(len(nCells))])
  cell_widths = coords[:, 1:] - coords[:, :-1]
  return coords, cell_widths


def read_velocity_old(case_directory, time_step, widths, periodic=[]):
  """Reads the velocity field at a given time-step.

  Arguments
  ---------
  case_directory -- directory of the simulation
  time_step -- time-step at which the field will be read
  widths -- cell-widths in each direction
  periodic -- list of directions with periodic boundary conditions (default [])
  """
  print('Reading the velocity field at time-step {} ...'.format(time_step))
  dx, dy, dz = widths[0], widths[1], (None if len(widths) < 3 else widths[2])
  # number of velocity nodes in each direction (depends on type of bc)
  nx, ny = dx.size, dy.size
  nxu, nyv = (nx if 'x' in periodic else nx-1), (ny if 'y' in periodic else ny-1)
  if dz != None:
    nz = dz.size
    nzw = (nz if 'z' in periodic else nz-1)
  # folder with numerical solution
  time_step_directory = '{}/{:0>7}'.format(case_directory, time_step)
  # x-flux
  flux_path = '{}/qx.dat'.format(time_step_directory)
  qx = PetscBinaryIO.PetscBinaryIO().readBinaryFile(flux_path)[0]
  # y-flux
  flux_path = '{}/qy.dat'.format(time_step_directory)
  qy = PetscBinaryIO.PetscBinaryIO().readBinaryFile(flux_path)[0]
  if dz == None:
    # x-velocity at x-velocity nodes
    u = qx.reshape((ny, nxu))/numpy.outer(dy, numpy.ones(nxu))
    # y-velocity at y-velocity nodes
    v = qy.reshape((nyv, nx))/numpy.outer(numpy.ones(nyv), dx)
    return u, v
  else:
    # z-flux
    flux_path = '{}/qz.dat'.format(time_step_directory)
    qz = PetscBinaryIO.PetscBinaryIO().readBinaryFile(flux_path)[0]
    u = qx.reshape((nz, ny, nxu))/ reduce(numpy.multiply, 
                                          numpy.ix_(dz, dy, numpy.ones(nxu)))
    v = qy.reshape((nz, nyv, nx))/ reduce(numpy.multiply, 
                                          numpy.ix_(dz, numpy.ones(nyv), dx))
    w = qz.reshape((nzw, ny, nx))/ reduce(numpy.multiply, 
                                          numpy.ix_(numpy.ones(nzw), dy, dx))
    return u, v, w


def read_velocity(case_directory, time_step, widths, periodic=[]):
  """Reads the velocity field at a given time-step.

  Arguments
  ---------
  case_directory -- directory of the simulation
  time_step -- time-step at which the field will be read
  widths -- cell-widths in each direction
  periodic -- list of directions with periodic boundary conditions (default [])
  """
  print('Reading the velocity field at time-step {} ...'.format(time_step))
  dx, dy, dz = widths[0], widths[1], (None if len(widths) < 3 else widths[2])
  # number of velocity nodes in each direction (depends on type of bc)
  nx, ny = dx.size, dy.size
  nxu, nyv = (nx if 'x' in periodic else nx-1), (ny if 'y' in periodic else ny-1)
  if dz != None:
    nz = dz.size
    nzw = (nz if 'z' in periodic else nz-1)
  # folder with numerical solution
  time_step_directory = '{}/{:0>7}'.format(case_directory, time_step)
  # x-flux
  flux_path = '{}/qx.dat'.format(time_step_directory)
  qx = PetscBinaryIO.PetscBinaryIO().readBinaryFile(flux_path)[0]
  # y-flux
  flux_path = '{}/qy.dat'.format(time_step_directory)
  qy = PetscBinaryIO.PetscBinaryIO().readBinaryFile(flux_path)[0]
  if dz == None:
    # x-velocity at x-velocity nodes
    u = qx.reshape((ny, nxu))/numpy.outer(dy, numpy.ones(nxu))
    # y-velocity at y-velocity nodes
    v = qy.reshape((nyv, nx))/numpy.outer(numpy.ones(nyv), dx)
    return u, v
  else:
    # z-flux
    flux_path = '{}/qz.dat'.format(time_step_directory)
    qz = PetscBinaryIO.PetscBinaryIO().readBinaryFile(flux_path)[0]
    u = qx.reshape((nz, ny, nxu))/ reduce(numpy.multiply, 
                                          numpy.ix_(dz, dy, numpy.ones(nxu)))
    v = qy.reshape((nz, nyv, nx))/ reduce(numpy.multiply, 
                                          numpy.ix_(dz, numpy.ones(nyv), dx))
    w = qz.reshape((nzw, ny, nx))/ reduce(numpy.multiply, 
                                          numpy.ix_(numpy.ones(nzw), dy, dx))
    return u, v, w


def read_pressure(case_directory, time_step, n):
  """Reads the pressure fields from file given the time-step.

  Arguments
  ---------
  case_directory -- directory of the simulation
  time_step -- time-step at which the field will be read
  n -- list of number of cells in each direction
  """
  print('Reading the pressure field at time-step {} ...'.format(time_step))
  # folder with numerical solution
  time_step_directory = '{}/{:0>7}'.format(case_directory, time_step)
  # pressure
  pressure_path = '{}/phi.dat'.format(time_step_directory)
  p = PetscBinaryIO.PetscBinaryIO().readBinaryFile(pressure_path)[0]
  if len(n) > 2:
    return p.reshape((n[2], n[1], nx[0]))
  else:
    return p.reshape((n[1], n[0]))


if __name__ == '__main__':
  pass