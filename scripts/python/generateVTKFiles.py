#!/usr/bin/env python

# file: generateVTKFiles.py
# author: Olivier Mesnard (mesnardo@gwu.edu), Anush Krishnan (anush@bu.edu)
# description: Converts PETSc output to .vtk format.


import os
import argparse

import numpy

import ioPetIBM


def read_inputs():
  """Parses the command-line."""
  # create parser
  parser = argparse.ArgumentParser(description='Converts PETSc output to VTK '
                                               'format for 3D case',
                        formatter_class= argparse.ArgumentDefaultsHelpFormatter)
  # fill parser with arguments
  parser.add_argument('--case', dest='case_directory', type=str, 
                      default=os.getcwd(), help='directory of the simulation')
  parser.add_argument('--variables', '-v', dest='variables', type=str, 
                      nargs='+', default=['velocity', 'pressure'],
                      help='list of variables to generate (velocity, pressure)')
  parser.add_argument('--bottom-left', '-bl', dest='bottom_left', type=float,
                      nargs='+', 
                      default=[float('-inf'), float('-inf'), float('-inf')],
                      help='coordinates of the bottom-left corner')
  parser.add_argument('--top-right', '-tr', dest='top_right', type=float,
                      nargs='+', 
                      default=[float('inf'), float('inf'), float('inf')],
                      help='coordinates of the top-right corner')
  parser.add_argument('--time-steps', '-t', dest='time_steps', type=int,
                      nargs='+', default=[None, None, None],
                      help='time-steps to convert (start, end, increment)')
  parser.add_argument('--stride', '-s', dest='stride', type=int, default=1,
                      help='stride at which vector are written')
  parser.add_argument('--periodic', '-p', dest='periodic', type=str, nargs='+',
                      default=[], help='direction(s) (x and/or y and/or z) '
                                       'with periodic boundary conditions')
  # parse command-line
  return parser.parse_args()


def interpolate_cell_centers(velocity):
  """Interpolates the velocity field at the cell-centers.

  Arguments
  ---------
  velocity -- velocity field on a staggered grid
  """
  dim3 = (True if len(velocity) == 3 else False)
  x_centers, y_centers = velocity[1]['x'][1:-1], velocity[0]['y'][1:-1]
  u, v = velocity[0]['values'], velocity[1]['values']
  if dim3:
    z_centers = velocity[0]['z'][1:-1]
    w = velocity[2]['values']
    u = 0.5*(u[1:-1, 1:-1, :-1] + u[1:-1, 1:-1, 1:])
    v = 0.5*(v[1:-1, :-1, 1:-1] + v[1:-1:, 1:, 1:-1])
    w = 0.5*(w[:-1, 1:-1, 1:-1] + w[1:, 1:-1, 1:-1])
    # tests
    assert (z_centers.size, y_centers.size, x_centers.size) == u.shape
    assert (z_centers.size, y_centers.size, x_centers.size) == v.shape
    assert (z_centers.size, y_centers.size, x_centers.size) == w.shape
    return [{'x': x_centers, 'y': y_centers, 'z': z_centers, 'values': u}, 
            {'x': x_centers, 'y': y_centers, 'z': z_centers, 'values': v}, 
            {'x': x_centers, 'y': y_centers, 'z': z_centers, 'values': w}]
  else:
    u = 0.5*(u[1:-1, :-1] + u[1:-1, 1:])
    v = 0.5*(v[:-1, 1:-1] + v[1:, 1:-1])
    # tests
    assert (y_centers.size, x_centers.size) == u.shape
    assert (y_centers.size, x_centers.size) == v.shape
    return [{'x': x_centers, 'y': y_centers, 'values': u},
            {'x': x_centers, 'y': y_centers, 'values': v}]


def main():
  """Converts PETSc output to .vtk format."""
  # parse command-line
  parameters = read_inputs()
  print ('[case-directory] %s' % parameters.case_directory)
  print ('[variables] %s' % parameters.variables)

  # list of time-steps to post-process
  time_steps = ioPetIBM.get_time_steps(parameters.case_directory, 
                                       parameters.time_steps)

  # read mesh grid
  coordinates = ioPetIBM.read_grid(parameters.case_directory)

  for time_step in time_steps:
    if 'velocity' in parameters.variables:
      velocity = ioPetIBM.read_velocity(parameters.case_directory, time_step, 
                                        coordinates, 
                                        periodic=parameters.periodic)
      # need to get velocity at cell-centers, not staggered arrangement
      velocity = interpolate_cell_centers(velocity)
      ioPetIBM.write_vtk(velocity, parameters.case_directory, time_step, 
                         name='velocity',
                         view=[parameters.bottom_left, parameters.top_right],
                         stride=parameters.stride)
    if 'pressure' in parameters.variables:
      pressure = ioPetIBM.read_pressure(parameters.case_directory, time_step, 
                                        coordinates)
      ioPetIBM.write_vtk(pressure, parameters.case_directory, time_step,
                         name='pressure',
                         view=[parameters.bottom_left, parameters.top_right],
                         stride=parameters.stride)

  print('\n[{}] DONE'.format(os.path.basename(__file__)))


if __name__ == "__main__":
  main()