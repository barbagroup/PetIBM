#!/usr/bin/env python

# file: ioPetIBM.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Collection of IO functions for PetIBM.


import os
import sys

import numpy
sys.path.append(os.path.join(os.environ['PETSC_DIR'], 'bin', 'pythonscripts'))
import PetscBinaryIO


def get_time_steps(case_directory, time_steps_range):
  """Returns a list of the time-steps to post-process.

  Arguments
  ---------
  case_directory -- directory of the simulation
  time_steps_range -- initial, final and stride of the time-steps
  """
  if any(time_steps_range):
    return range(time_steps_range[0],
                 time_steps_range[1]+1,
                 time_steps_range[2])
  else:
    return sorted(int(folder) for folder in os.listdir(case_directory)
                              if folder[0] == '0')


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
    m_x, xu = numpy.where(~numpy.ma.getmaskarray(xu))[0], xu.compressed()
    m_y, yu = numpy.where(~numpy.ma.getmaskarray(yu))[0], yu.compressed()
    m_z, zu = numpy.where(~numpy.ma.getmaskarray(zu))[0], zu.compressed()
    # compute x-velocity field
    u = qx.reshape((nz, ny, nxu))/ reduce(numpy.multiply, 
                                          numpy.ix_(dz, dy, numpy.ones(nxu)))
    u = u[m_z[0]:m_z[-1]+1, m_y[0]:m_y[-1]+1, m_x[0]:m_x[-1]+1]
    # get third-dimension coordinate of y-velocity nodes
    zv = 0.5*(z[:-1]+z[1:])
    # get mask to apply
    m_x, xv = numpy.where(~numpy.ma.getmaskarray(xv))[0], xv.compressed()
    m_y, yv = numpy.where(~numpy.ma.getmaskarray(yv))[0], yv.compressed()
    m_z, zv = numpy.where(~numpy.ma.getmaskarray(zv))[0], zv.compressed()
    # compute y-velocity field
    v = qy.reshape((nz, nyv, nx))/ reduce(numpy.multiply, 
                                          numpy.ix_(dz, numpy.ones(nyv), dx))
    v = v[m_z[0]:m_z[-1]+1, m_y[0]:m_y[-1]+1, m_x[0]:m_x[-1]+1]
    # read z-flux
    flux_path = '{}/qz.dat'.format(time_step_directory)
    qz = PetscBinaryIO.PetscBinaryIO().readBinaryFile(flux_path)[0]
    # get coordinates of z-velocity nodes
    xw, yw, zw = 0.5*(x[:-1]+x[1:]), 0.5*(y[:-1]+y[1:]), z[1:-1]
    # get mask to apply
    m_x, xw = numpy.where(~numpy.ma.getmaskarray(xw))[0], xw.compressed()
    m_y, yw = numpy.where(~numpy.ma.getmaskarray(yw))[0], yw.compressed()
    m_z, zw = numpy.where(~numpy.ma.getmaskarray(zw))[0], zw.compressed()
    # compute z-velocity field
    w = qz.reshape((nzw, ny, nx))/ reduce(numpy.multiply, 
                                          numpy.ix_(numpy.ones(nzw), dy, dx))
    w = w[m_z[0]:m_z[-1]+1, m_y[0]:m_y[-1]+1, m_x[0]:m_x[-1]+1]
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


def write_velocity_vtk2d(u, v, case_directory, time_step):
  """Writes the velocity fields in a .vtk file.

  Arguments
  ---------
  u, v -- velocity fields to write
  case_directory -- directory of the simulation
  time_step -- time-step to write
  """
  print('Write .vtk file for 2d velocity fields ...')
  # create directory where .vtk file will be saved
  vtk_directory = '{}/vtk_files/velocity'.format(case_directory)
  if not os.path.isdir(vtk_directory):
    print('Make directory: {}'.format(vtk_directory))
    os.makedirs(vtk_directory)
  vtk_file_path = '{}/velocity{:0>7}.vtk'.format(vtk_directory, time_step)
  # get velocity fields
  u, v, x, y = u['values'], v['values'], u['x'], u['y']
  # write .vtk file
  with open(vtk_file_path, 'w') as outfile:
    outfile.write('# vtk DataFile Version 3.0\n')
    outfile.write('contains 2d velocity fields\n')
    outfile.write('ASCII\n')
    outfile.write('DATASET RECTILINEAR_GRID\n')
    outfile.write('DIMENSIONS {} {} 1\n'.format(x.size, y.size))
    outfile.write('X_COORDINATES {} double\n'.format(x.size))
    numpy.savetxt(outfile, x, fmt='%f')
    outfile.write('Y_COORDINATES {} double\n'.format(y.size))
    numpy.savetxt(outfile, y, fmt='%f')
    outfile.write('Z_COORDINATES 1 double\n0.0\n')
    outfile.write('POINT_DATA {}\n'.format(x.size*y.size))
    outfile.write('\nVECTORS velocity double\n')
    numpy.savetxt(outfile, 
                  numpy.c_[u.flatten(), v.flatten(), numpy.zeros(u.size)],
                  fmt='%.6f', delimiter='\t')


def write_velocity_vtk3d(u, v, w, case_directory, time_step):
  """Writes the velocity fields in a .vtk file.

  Arguments
  ---------
  u, v, w -- velocity fields to write
  case_directory -- directory of the simulation
  time_step -- time-step to write
  """
  print('Write .vtk file for 3d velocity fields ...')
  # create directory where .vtk file will be saved
  vtk_directory = '{}/vtk_files/velocity'.format(case_directory)
  if not os.path.isdir(vtk_directory):
    print('Make directory: {}'.format(vtk_directory))
    os.makedirs(vtk_directory)
  vtk_file_path = '{}/velocity{:0>7}.vtk'.format(vtk_directory, time_step)
  # get velocity fields
  u, v, w, x, y, z = u['values'], v['values'], w['values'], u['x'], u['y'], u['z']
  # write .vtk file
  with open(vtk_file_path, 'w') as outfile:
    outfile.write('# vtk DataFile Version 3.0\n')
    outfile.write('contains 3d velocity fields\n')
    outfile.write('ASCII\n')
    outfile.write('DATASET RECTILINEAR_GRID\n')
    outfile.write('DIMENSIONS {} {} {}\n'.format(x.size, y.size, z.size))
    outfile.write('X_COORDINATES {} double\n'.format(x.size))
    numpy.savetxt(outfile, x, fmt='%f')
    outfile.write('Y_COORDINATES {} double\n'.format(y.size))
    numpy.savetxt(outfile, y, fmt='%f')
    outfile.write('Z_COORDINATES {} double\n'.format(z.size))
    numpy.savetxt(outfile, z, fmt='%f')
    outfile.write('POINT_DATA {}\n'.format(x.size*y.size*z.size))
    outfile.write('\nVECTORS velocity double\n')
    numpy.savetxt(outfile, 
                  numpy.c_[u.flatten(), v.flatten(), w.flatten()],
                  fmt='%.6f', delimiter='\t')


def write_pressure_vtk2d(p, case_directory, time_step):
  """Writes the pressure field in a .vtk file.

  Arguments
  ---------
  p -- pressure field
  case_directory -- directory of the simulation
  time_step -- time-step to write
  """
  print('Write .vtk file for 2d pressure field ...')
  # create directory where .vtk file will be saved
  vtk_directory = '{}/vtk_files/pressure'.format(case_directory)
  if not os.path.isdir(vtk_directory):
    print('Make directory: {}'.format(vtk_directory))
    os.makedirs(vtk_directory)
  vtk_file_path = '{}/pressure{:0>7}.vtk'.format(vtk_directory, time_step)
  # get pressure field
  p, x, y = p['values'], p['x'], p['y']
  with open(vtk_file_path, 'w') as outfile:
    outfile.write('# vtk DataFile Version 3.0\n')
    outfile.write('contains 2d pressure field\n')
    outfile.write('ASCII\n')
    outfile.write('DATASET RECTILINEAR_GRID\n')
    outfile.write('DIMENSIONS {} {} 1\n'.format(x.size, y.size))
    outfile.write('X_COORDINATES {} double\n'.format(x.size))
    numpy.savetxt(outfile, x, fmt='%f')
    outfile.write('Y_COORDINATES {} double\n'.format(y.size))
    numpy.savetxt(outfile, y, fmt='%f')
    outfile.write('Z_COORDINATES 1 double\n0.0\n')
    outfile.write('POINT_DATA {}\n'.format(x.size*y.size))
    outfile.write('\nSCALARS pressure double 1\nLOOKUP_TABLE default\n')
    numpy.savetxt(outfile, p.flatten(), fmt='%.6f', delimiter='\t')


def write_pressure_vtk3d(p, case_directory, time_step):
  """Writes the pressure field in a .vtk file.

  Arguments
  ---------
  p -- pressure field
  case_directory -- directory of the simulation
  time_step -- time-step to write
  """
  print('Write .vtk file for 3d pressure field ...')
  # create directory where .vtk file will be saved
  vtk_directory = '{}/vtk_files/pressure'.format(case_directory)
  if not os.path.isdir(vtk_directory):
    print('Make directory: {}'.format(vtk_directory))
    os.makedirs(vtk_directory)
  vtk_file_path = '{}/pressure{:0>7}.vtk'.format(vtk_directory, time_step)
  # get pressure field
  p, x, y, z = p['values'], p['x'], p['y'], p['z']
  with open(vtk_file_path, 'w') as outfile:
    outfile.write('# vtk DataFile Version 3.0\n')
    outfile.write('contains 3d pressure field\n')
    outfile.write('ASCII\n')
    outfile.write('DATASET RECTILINEAR_GRID\n')
    outfile.write('DIMENSIONS {} {} {}\n'.format(x.size, y.size, z.size))
    outfile.write('X_COORDINATES {} double\n'.format(x.size))
    numpy.savetxt(outfile, x, fmt='%f')
    outfile.write('Y_COORDINATES {} double\n'.format(y.size))
    numpy.savetxt(outfile, y, fmt='%f')
    outfile.write('Z_COORDINATES {} double\n'.format(z.size))
    outfile.write('POINT_DATA {}\n'.format(x.size*y.size*z.size))
    outfile.write('\nSCALARS pressure double 1\nLOOKUP_TABLE default\n')
    numpy.savetxt(outfile, p.flatten(), fmt='%.6f', delimiter='\t')


if __name__ == '__main__':
  pass