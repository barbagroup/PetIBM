#!/usr/bin/env python

# file: generate3DVTKFile.py
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: Converts PETSc output to VTK format for 3D cases.


import os
import sys
import argparse

import numpy
sys.path.append(os.path.join(os.environ['PETSC_DIR'], 'bin', 'pythonscripts'))
import PetscBinaryIO


def read_inputs():
  """Parses the command-line."""
  # create parser
  parser = argparse.ArgumentParser(description='Converts PETSc output to VTK '
                                               'format for 2D case',
                        formatter_class= argparse.ArgumentDefaultsHelpFormatter)
  # fill parser with arguments
  parser.add_argument('--case', dest='case_directory', type=str, 
                      default=os.getcwd(), help='directory of the simulation')
  parser.add_argument('--bottom-left', '-bl', dest='bottom_left', type=float,
                      nargs='+', 
                      default=[float('-inf'), float('-inf'), float('-inf')],
                      help='coordinates of the bottom-left corner')
  parser.add_argument('--top-right', '-tr', dest='top_right', type=float,
                      nargs='+', 
                      default=[float('inf'), float('inf'), float('inf')],
                      help='coordinates of the top-right corner')
  parser.add_argument('--time-steps', '-t', dest='time_steps', type=float,
                      nargs='+', default=[None, None, None],
                      help='time-steps to convert (start, end, increment)')
  parser.add_argument('--stride', '-s', dest='stride', type=int, default=1,
                      help='stride at which vector are written')
  parser.add_argument('--periodic', '-p', dest='periodic', type=str, nargs='+',
                      default=[], help='direction(s) (x and/or y and/or z) '
                                       'with periodic boundary conditions')
  return parser.parse_args()

def main():
  """Converts PETSc output to VTK format for 2D case."""
  # parse the command-line
  args = read_inputs()
  print ('[case-directory] %s' % args.case_directory)

  with open('%s/grid.txt' % args.case_directory, 'r') as infile:
    nx, ny, nz = [int(v) for v in infile.readline().strip().split()]
    grid = numpy.loadtxt(infile, dtype=float)

  # compute cell-vertex coordinates
  x, y, z = grid[:nx+1], grid[nx+1:nx+1+ny+1], grid[nx+1+ny+1:]
  # compute cell-widths
  dx, dy, dz = x[1:]-x[:-1], y[1:]-y[:-1], z[1:]-z[:-1]

  # number of velocity nodes in each direction (depends on type of bc)
  nxu, nyu, nzu = (nx if 'x' in args.periodic else nx-1), ny, nz
  nxv, nyv, nzv = nx, (ny if 'y' in args.periodic else ny-1), nz
  nxw, nyw, nzw = nx, ny, (nz if 'z' in args.periodic else nz-1)

  # calculate cell-center coordinates
  x = 0.5*(x[1:nxu]+x[2:nxu+1])
  y = 0.5*(y[1:nyv]+y[2:nyv+1])
  z = 0.5*(z[1:nzw]+z[2:nzw+1])

  # get masks to account for boundary-limits and stride
  mask_x = numpy.where(numpy.logical_and(x >= args.bottom_left[0], 
                                         x <= args.top_right[0]))[0][::args.stride]
  mask_y = numpy.where(numpy.logical_and(y >= args.bottom_left[1], 
                                         y <= args.top_right[1]))[0][::args.stride]
  mask_z = numpy.where(numpy.logical_and(z >= args.bottom_left[2], 
                                         z <= args.top_right[2]))[0][::args.stride]
  
  # apply mask on cell-center coordinates
  x = x[mask_x]
  y = y[mask_y]
  z = z[mask_z]

  # get time-steps to write .vtk files
  if any(args.time_steps):
    time_steps = range(args.time_steps[0], 
                       args.time_steps[1]+1, 
                       args.time_steps[2])
  else:
    time_steps = sorted(int(folder) for folder in os.listdir(args.case_directory)
                                    if folder[0] == '0')

  # create directory where .vtk files will be saved
  vtk_directory = '%s/vtk_files' % args.case_directory
  print ('[vtk-directory] %s' % vtk_directory)
  if not os.path.isdir(vtk_directory):
    os.makedirs(vtk_directory)

  for time_step in time_steps:
    # read fluxes in x-direction
    qx = PetscBinaryIO.PetscBinaryIO().readBinaryFile('%s/%07d/qx.dat' 
                                                      % (args.case_directory,
                                                         time_step))[0]
    qx = qx.reshape((nzu, nyu, nxu))
    # compute u-velocity at cell-centers
    u = ( 0.5 * (qx[1:nzw, 1:nyv, :nxu-1] + qx[1:nzw, 1:nyv, 1:nxu])
              / reduce(numpy.multiply, 
                       numpy.ix_(*[dz[1:nzw], dy[1:nyv], numpy.ones(nxu-1)])) )
    # apply mask for boundary limits and stride
    u = numpy.array([[u[k][j][mask_x] for j in mask_y] for k in mask_z])
    # read fluxes in y-direction
    qy = PetscBinaryIO.PetscBinaryIO().readBinaryFile('%s/%07d/qy.dat' 
                                                      % (args.case_directory,
                                                         time_step))[0]
    qy = qy.reshape((nzv, nyv, nxv))
    # compute v-velocity at cell-centers
    v = ( 0.5 * (qy[1:nzw, :nyv-1, 1:nxu] + qy[1:nzw, 1:nyv, 1:nxu])
              / reduce(numpy.multiply, 
                       numpy.ix_(*[dz[1:nzw], numpy.ones(nyv-1), dx[1:nxu]])) )
    # apply mask for boundary limits and stride
    v = numpy.array([[v[k][j][mask_x] for j in mask_y] for k in mask_z])
    # read fluxes in z-direction
    qz = PetscBinaryIO.PetscBinaryIO().readBinaryFile('%s/%07d/qz.dat' 
                                                      % (args.case_directory,
                                                         time_step))[0]
    qz = qz.reshape((nzw, nyw, nxw))
    # compute w-velocity at cell-centers
    w = ( 0.5 * (qy[:nzw-1, 1:nyv, 1:nxu] + qy[1:nzw, 1:nyv, 1:nxu])
              / reduce(numpy.multiply, 
                       numpy.ix_(*[numpy.ones(nzw-1), dy[1:nyv], dx[1:nxu]])) )
    # apply mask for boundary limits and stride
    w = numpy.array([[w[k][j][mask_x] for j in mask_y] for k in mask_z])

    print ('writing vtk file at time-step %d ...' % time_step)
    vtk_file_path = '%s/velocity%07d.vtk' % (vtk_directory, time_step)
    with open(vtk_file_path, 'w') as outfile:
      outfile.write('# vtk DataFile Version 3.0\n')
      outfile.write('Header\n')
      outfile.write('ASCII\n')
      outfile.write('DATASET RECTILINEAR_GRID\n')
      outfile.write('DIMENSIONS %d %d %d\n' % (x.size, y.size, z.size))
      outfile.write('X_COORDINATES %d double\n' % x.size)
      numpy.savetxt(outfile, x, fmt='%f')
      outfile.write('Y_COORDINATES %d double\n' % y.size)
      numpy.savetxt(outfile, y, fmt='%f')
      outfile.write('Z_COORDINATES %d double\n' % z.size)
      numpy.savetxt(outfile, z, fmt='%f')
      outfile.write('POINT_DATA %d\n' % (x.size*y.size*z.size))
      outfile.write('VECTORS velocity double\n')
      numpy.savetxt(outfile, 
                    numpy.c_[u.flatten(), v.flatten(), w.flatten()],
                    fmt='%.6f', delimiter='\t')


if __name__ == "__main__":
  main()