#!/usr/bin/env python

# file: generateVTKFiles2d.py
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: Converts PETSc output to VTK format for 2D case.


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
  parser.add_argument('--variables', '-v', dest='variables', type=str, 
                      nargs='+', default=['velocity', 'pressure'],
                      help='list of variables to generate (velocity, pressure)')
  parser.add_argument('--bottom-left', '-bl', dest='bottom_left', type=float,
                      nargs='+', default=[float('-inf'), float('-inf')],
                      help='coordinates of the bottom-left corner')
  parser.add_argument('--top-right', '-tr', dest='top_right', type=float,
                      nargs='+', default=[float('inf'), float('inf')],
                      help='coordinates of the top-right corner')
  parser.add_argument('--time-steps', '-t', dest='time_steps', type=float,
                      nargs='+', default=[None, None, None],
                      help='time-steps to convert (start, end, increment)')
  parser.add_argument('--stride', '-s', dest='stride', type=int, default=1,
                      help='stride at which vector are written')
  parser.add_argument('--periodic', '-p', dest='periodic', type=str, nargs='+',
                      default=[], help='direction(s) (x and/or y) with '
                                       'periodic boundary conditions')
  return parser.parse_args()

def main():
  """Converts PETSc output to VTK format for 2D case."""
  # parse the command-line
  args = read_inputs()
  print ('[case-directory] %s' % args.case_directory)
  print ('[variables] %s' % args.variables)

  with open('%s/grid.txt' % args.case_directory, 'r') as infile:
    nx, ny = [int(v) for v in infile.readline().strip().split()]
    grid = numpy.loadtxt(infile, dtype=float)
  
  # compute cell-vertex coordinates
  x, y = grid[:nx+1], grid[nx+1:]
  # compute cell-widths
  dx, dy = x[1:]-x[:-1], y[1:]-y[:-1]

  # number of velocity nodes in each direction (depends on type of bc)
  nxu, nyu = (nx if 'x' in args.periodic else nx-1), ny
  nxv, nyv = nx, (ny if 'y' in args.periodic else ny-1)

  # calculate cell-center coordinates
  x = 0.5*(x[1:nxu]+x[2:nxu+1])
  y = 0.5*(y[1:nyv]+y[2:nyv+1])

  # get masks to account for boundary-limits and stride
  mask_x = numpy.where(numpy.logical_and(x >= args.bottom_left[0], 
                                         x <= args.top_right[0]))[0][::args.stride]
  mask_y = numpy.where(numpy.logical_and(y >= args.bottom_left[1], 
                                         y <= args.top_right[1]))[0][::args.stride]
  
  # apply mask on cell-center coordinates
  x = x[mask_x]
  y = y[mask_y]

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
    if 'velocity' in args.variables:
      # read fluxes in x-direction
      qx = PetscBinaryIO.PetscBinaryIO().readBinaryFile('%s/%07d/qx.dat' 
                                                        % (args.case_directory,
                                                           time_step))[0]
      qx = qx.reshape((nyu, nxu))
      # compute u-velocity at cell-centers
      u = ( 0.5 * (qx[1:nyv, :nxu-1] + qx[1:nyv, 1:nxu])
                / numpy.outer(dy[1:nyv], numpy.ones(nxu-1)) )
      # apply mask for boundary limits and stride
      u = numpy.array([u[j][mask_x] for j in mask_y])
      # read fluxes in y-direction
      qy = PetscBinaryIO.PetscBinaryIO().readBinaryFile('%s/%07d/qy.dat' 
                                                        % (args.case_directory,
                                                           time_step))[0]
      qy = qy.reshape((nyv, nxv))
      # compute v-velocity at cell-centers
      v = ( 0.5 * (qy[:nyv-1, 1:nxu] + qy[1:nyv, 1:nxu])
                / numpy.outer(numpy.ones(nyv-1), dx[1:nxu]) )
      # apply mask for boundary limits and stride
      v = numpy.array([v[j][mask_x] for j in mask_y])
    if 'pressure' in args.variables:
      # read pressure field
      p = PetscBinaryIO.PetscBinaryIO().readBinaryFile('%s/%07d/phi.dat'
                                                       % (args.case_directory,
                                                          time_step))[0]
      p = p.reshape((ny, nx))
      # apply mask for boundary limits and stride
      p = numpy.array([p[j][mask_x] for j in mask_y])

    print ('writing .vtk file at time-step %d ...' % time_step)
    vtk_file_path = '%s/fields%07d.vtk' % (vtk_directory, time_step)
    with open(vtk_file_path, 'w') as outfile:
      outfile.write('# vtk DataFile Version 3.0\n')
      outfile.write('contains velocity and pressure fields\n')
      outfile.write('ASCII\n')
      outfile.write('DATASET RECTILINEAR_GRID\n')
      outfile.write('DIMENSIONS %d %d 1\n' % (x.size, y.size))
      outfile.write('X_COORDINATES %d double\n' % x.size)
      numpy.savetxt(outfile, x, fmt='%f')
      outfile.write('Y_COORDINATES %d double\n' % y.size)
      numpy.savetxt(outfile, y, fmt='%f')
      outfile.write('Z_COORDINATES 1 double\n0.0\n')
      outfile.write('POINT_DATA %d\n' % (x.size*y.size))
      if 'velocity' in args.variables:
        outfile.write('\nVECTORS velocity double\n')
        numpy.savetxt(outfile, 
                      numpy.c_[u.flatten(), v.flatten(), numpy.zeros(u.size)],
                      fmt='%.6f', delimiter='\t')
      if 'pressure' in args.variables:
        outfile.write('\nSCALARS pressure double 1\nLOOKUP_TABLE default\n')
        numpy.savetxt(outfile, p.flatten(), fmt='%.6f', delimiter='\t')


if __name__ == "__main__":
  main()