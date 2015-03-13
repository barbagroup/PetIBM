#!/usr/bin/env python

# file: generateVTKFiles2d.py
# author: Olivier Mesnard (mesnardo@gwu.edu), Anush Krishnan (anush@bu.edu)
# description: Converts PETSc output to VTK format for 2D case.


import os
import argparse

import numpy

import ioPetIBM


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
  parameters = read_inputs()
  print('[case-directory] {}'.format(args.case_directory))
  print('[variables] {}'.format(parameters.variables))

  # get time-steps to write .vtk files
  if any(args.time_steps):
    time_steps = range(parameters.time_steps[0], 
                       parameters.time_steps[1]+1, 
                       parameters.time_steps[2])
  else:
    time_steps = sorted(int(folder) for folder 
                                    in os.listdir(parameters.case_directory)
                                    if folder[0] == '0')

  # create directory where .vtk files will be saved
  vtk_directory = '%s/vtk_files' % args.case_directory
  print('[vtk-directory] {}'.format(vtk_directory))
  if not os.path.isdir(vtk_directory):
    os.makedirs(vtk_directory)

  # read mesh grid
  x, y = ioPetIBM.read_grid(parameters.case_directory)
  x_centers, y_centers = 0.5*(x[:-1]+x[1:]), 0.5*(y[:-1]+y[1:])

  for time_step in time_steps:
    if 'velocity' in parameters.variables:
      u, v = ioPetIBM.read_velocity(parameters.case_directory, time_step, [x, y],
                                    periodic=parameters.periodic)
    if 'pressure' in parameters.variables:
      p = ioPetIBM.read_pressure(parameters.case_directory, time_step, [x, y])

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