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


def interpolate_cell_centers(u, v):
  """Interpolates the velocity fields at the cell-centers.

  Arguments
  ---------
  u, v -- x- and y- velocity fields on a staggered grid
  x, y -- coordinates of the mesh grid
  """
  u, xu, yu = u['values'], u['x'], u['y']
  v, xv, yv = v['values'], v['x'], v['y']
  mask_x = numpy.where(numpy.logical_and(xv > xu[0], xv < xu[-1]))[0]
  mask_y = numpy.where(numpy.logical_and(yu > yv[0], yu < yv[-1]))[0]
  x_centers, y_centers = xv[mask_x], yu[mask_y]
  u, v = 0.5*(u[mask_y, :-1]+u[mask_y, 1:]), 0.5*(v[:-1, mask_x]+v[1:, mask_x])
  return ( {'x': x_centers, 'y': y_centers, 'values': u}, 
           {'x': x_centers, 'y': y_centers, 'values': v} )


def main():
  """Converts PETSc output to VTK format for 2D case."""
  # parse command-line
  parameters = read_inputs()
  print('[case-directory] {}'.format(parameters.case_directory))
  print('[variables] {}'.format(parameters.variables))

  # list of time-steps to post-process
  time_steps = ioPetIBM.get_time_steps(parameters.case_directory, 
                                       parameters.time_steps)

  # read mesh grid
  x, y = ioPetIBM.read_grid(parameters.case_directory, 
                            bottom_left=parameters.bottom_left,
                            top_right=parameters.top_right)

  for time_step in time_steps:
    if 'velocity' in parameters.variables:
      u, v = ioPetIBM.read_velocity(parameters.case_directory, time_step, [x, y],
                                    periodic=parameters.periodic)
      # need to get values at cell-centers, not staggered arrangement
      u, v = interpolate_cell_centers(u, v)
      ioPetIBM.write_velocity_vtk2d(u, v, parameters.case_directory, time_step)
    if 'pressure' in parameters.variables:
      p = ioPetIBM.read_pressure(parameters.case_directory, time_step, [x, y])
      ioPetIBM.write_pressure_vtk2d(p, parameters.case_directory, time_step)

  print('\n[{}] DONE'.format(os.path.basename(__file__)))


if __name__ == "__main__":
  main()