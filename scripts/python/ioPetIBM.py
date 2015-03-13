#!/usr/bin/env python

# file: ioPetIBM.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Collection of IO functions for PetIBM.


import os
import sys

import numpy
sys.path.append(os.path.join(os.environ['PETSC_DIR'], 'bin', 'pythonscripts'))
import PetscBinaryIO


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
  print('Read the mesh grid ...')
  grid_path = '{}/grid.txt'.format(case_directory)
  with open(grid_path, 'r') as infile:
    nCells = numpy.array([int(n) for n in infile.readline().strip().split()])
    nodes = numpy.loadtxt(infile, dtype=float)
  coords = numpy.array(numpy.split(nodes, numpy.cumsum(nCells[:-1]+1)))
  return numpy.ma.masked_array([numpy.ma.masked_outside(coords[i], 
                                                          bottom_left[i], 
                                                          top_right[i], 
                                                          copy=False) 
                                  for i in xrange(len(nCells))])


def read_velocity(case_directory, time_step, coords, periodic=[]):
  """Reads the velocity field at a given time-step.

  Arguments
  ---------
  case_directory -- directory of the simulation
  time_step -- time-step at which the field will be read
  coords -- coordinates in each direction
  periodic -- list of directions with periodic boundary conditions (default [])
  """
  print('Read the velocity field at time-step {} ...'.format(time_step))
  dim3 = (True if len(coords) == 3 else False)
  x, y, z = coords[0], coords[1], (None if not dim3 else coords[2])
  # compute cell-widths
  dx, dy, dz = x[1:]-x[:-1], y[1:]-y[:-1], (None if not dim3 else z[1:]-z[:-1])
  # number of of cells
  nx, ny, nz = dx.size, dy.size, (None if not dim3 else dz.size)
  # number of nodes depends on the boundary condition type
  nxu = (nx if 'x' in periodic else nx-1)
  nyv = (ny if 'y' in periodic else ny-1)
  nzw = (nz if 'z' in periodic or not nz else nz-1)
  # folder with numerical solution
  time_step_directory = '{}/{:0>7}'.format(case_directory, time_step)
  # read x-flux
  flux_path = '{}/qx.dat'.format(time_step_directory)
  qx = PetscBinaryIO.PetscBinaryIO().readBinaryFile(flux_path)[0]
  # read y-flux
  flux_path = '{}/qy.dat'.format(time_step_directory)
  qy = PetscBinaryIO.PetscBinaryIO().readBinaryFile(flux_path)[0]
  # get velocity nodes coordinates
  xu, yu = x[1:-1], 0.5*(y[:-1]+y[1:])
  xv, yv = 0.5*(x[:-1]+x[1:]), y[1:-1]
  if dim3:
    # get third-dimension coordinate of x-velocity nodes
    zu = 0.5*(z[:-1]+z[1:])
    # get mask to apply
    mask_x, xu = ~numpy.ma.getmaskarray(xu), xu.compressed()
    mask_y, yu = ~numpy.ma.getmaskarray(yu), yu.compressed()
    mask_z, zu = ~numpy.ma.getmaskarray(zu), zu.compressed()
    mask = reduce(numpy.multiply, numpy.ix_(mask_z, mask_y, mask_x))
    # compute x-velocity field
    u = qx.reshape((nz, ny, nxu))/ reduce(numpy.multiply, 
                                          numpy.ix_(dz, dy, numpy.ones(nxu)))
    u = u[mask].reshape((zu.size, yu.size, xu.size))
    # get third-dimension coordinate of y-velocity nodes
    zv = 0.5*(z[:-1]+z[1:])
    # get mask to apply
    mask_x, xv = ~numpy.ma.getmaskarray(xv), xv.compressed()
    mask_y, yv = ~numpy.ma.getmaskarray(yv), yv.compressed()
    mask_z, zv = ~numpy.ma.getmaskarray(zv), zv.compressed()
    mask = reduce(numpy.multiply, numpy.ix_(mask_z, mask_y, mask_x))
    # compute y-velocity field
    v = qy.reshape((nz, nyv, nx))/ reduce(numpy.multiply, 
                                          numpy.ix_(dz, numpy.ones(nyv), dx))
    v = v[mask].reshape((zv.size, yv.size, xv.size))
    # read z-flux
    flux_path = '{}/qz.dat'.format(time_step_directory)
    qz = PetscBinaryIO.PetscBinaryIO().readBinaryFile(flux_path)[0]
    # get coordinates of z-velocity nodes
    xw, yw, zw = 0.5*(x[:-1]+x[1:]), 0.5*(y[:-1]+y[1:]), z[1:-1]
    # get mask to apply
    mask_x, xw = ~numpy.ma.getmaskarray(xw), xw.compressed()
    mask_y, yw = ~numpy.ma.getmaskarray(yw), yw.compressed()
    mask_z, zw = ~numpy.ma.getmaskarray(zw), zw.compressed()
    mask = reduce(numpy.multiply, numpy.ix_(mask_z, mask_y, mask_x))
    # compute z-velocity field
    w = qz.reshape((nzw, ny, nx))/ reduce(numpy.multiply, 
                                          numpy.ix_(numpy.ones(nzw), dy, dx))
    w = w[mask].reshape((zw.size, yw.size, xw.size))
    return ({'x': xu, 'y': yu, 'z': zu, 'values': u}, 
            {'x': xv, 'y': yv, 'z': zv, 'values': v}, 
            {'x': xw, 'y': yw, 'z': zw, 'values': w})
  else:
    # get mask to apply
    mask_x, xu = ~numpy.ma.getmaskarray(xu), xu.compressed()
    mask_y, yu = ~numpy.ma.getmaskarray(yu), yu.compressed()
    mask = numpy.outer(mask_y, mask_x)
    # compute x-velocity field
    u = qx.reshape((ny, nxu))/numpy.outer(dy, numpy.ones(nxu))
    u = u[mask].reshape((yu.size, xu.size))
    # get mask to apply
    mask_x, xv = ~numpy.ma.getmaskarray(xv), xv.compressed()
    mask_y, yv = ~numpy.ma.getmaskarray(yv), yv.compressed()
    mask = numpy.outer(mask_y, mask_x)
    # compute y-velocity field
    v = qy.reshape((nyv, nx))/numpy.outer(numpy.ones(nyv), dx)
    v = v[mask].reshape((yv.size, xv.size))
    return ({'x': xu, 'y': yu, 'values': u}, 
            {'x': xv, 'y': yv, 'values': v})


def read_pressure(case_directory, time_step, coords):
  """Reads the pressure fields from file given the time-step.

  Arguments
  ---------
  case_directory -- directory of the simulation
  time_step -- time-step at which the field will be read
  n -- list of number of cells in each direction
  """
  print('Read the pressure field at time-step {} ...'.format(time_step))
  dim3 = (True if len(coords) == 3 else False)
  x, y, z = coords[0], coords[1], (None if not dim3 else coords[2])
  # folder with numerical solution
  time_step_directory = '{}/{:0>7}'.format(case_directory, time_step)
  # pressure
  pressure_path = '{}/phi.dat'.format(time_step_directory)
  p = PetscBinaryIO.PetscBinaryIO().readBinaryFile(pressure_path)[0]
  # get pressure nodes coordinates
  xp, yp = 0.5*(x[:-1]+x[1:]), 0.5*(y[:-1]+y[1:])
  nx, ny = xp.size, yp.size
  # get mask to apply
  mask_x, xp = ~numpy.ma.getmaskarray(xp), xp.compressed()
  mask_y, yp = ~numpy.ma.getmaskarray(yp), yp.compressed()
  if dim3:
    # get third-dimension coordinates of pressure nodes
    zp = 0.5*(z[:-1]+z[1:])
    nz = zp.size
    # get mask to apply
    mask_z, zp = ~numpy.ma.getmaskarray(zp), zp.compressed()
    mask = reduce(numpy.multiply, numpy.ix_(mask_z, mask_y, mask_x))
    # compute pressure field
    p = p.reshape((nz, ny, nx))[mask].reshape((zp.size, yp.size, xp.size))
    return {'x': xp, 'y':yp, 'z': zp, 'values': p}
  else:
    # get mask to apply
    mask = numpy.outer(mask_y, mask_x)
    # compute pressure field
    p = p.reshape((ny, nx))[mask].reshape((yp.size, xp.size))
    return {'x': xp, 'y': yp, 'values': p}


if __name__ == '__main__':
  pass