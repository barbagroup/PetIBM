#!/usr/bin/env python

# file: generateVTKFiles3d.py
# author: Olivier Mesnard (mesnardo@gwu.edu), Anush Krishnan (anush@bu.edu)
# description: Converts PETSc output to VTK format for 3D cases.


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
  parser.add_argument('--time-steps', '-t', dest='time_steps', type=float,
                      nargs='+', default=[None, None, None],
                      help='time-steps to convert (start, end, increment)')
  parser.add_argument('--stride', '-s', dest='stride', type=int, default=1,
                      help='stride at which vector are written')
  parser.add_argument('--periodic', '-p', dest='periodic', type=str, nargs='+',
                      default=[], help='direction(s) (x and/or y and/or z) '
                                       'with periodic boundary conditions')
  # parse command-line
  return parser.parse_args()


def interpolate_cell_centers(u, v, w, x, y, z):
  """Interpolates the velocity fields at the cell-centers.

  Arguments
  ---------
  u, v, w -- x-, y- and z- velocity fields on a staggered grid
  x, y, z -- coordinates of the mesh grid
  """
  u, xu, yu, zu = u['values'], u['x'], u['y'], u['z']
  v, xv, yv, zv = v['values'], v['x'], v['y'], v['z']
  w, xw, yw, zw = w['values'], w['x'], w['y'], w['z']
  mask_x = numpy.where(numpy.logical_and(xv > xu[0], xv < xu[-1]))[0]
  mask_y = numpy.where(numpy.logical_and(yu > yv[0], yu < yv[-1]))[0]
  mask_z = numpy.where(numpy.logical_and(zu > zw[0], zu < zw[-1]))[0]
  x_centers, y_centers, z_centers = xv[mask_x], yu[mask_y], zu[mask_z]
  u = 0.5*(u[mask_z[0]:mask_z[-1]+1, mask_y[0]:mask_y[-1]+1, :-1] + 
           u[mask_z[0]:mask_z[-1]+1, mask_y[0]:mask_y[-1]+1, 1:])
  v = 0.5*(v[mask_z[0]:mask_z[-1]+1, :-1, mask_x[0]:mask_x[-1]+1] + 
           v[mask_z[0]:mask_z[-1]+1, 1:, mask_x[0]:mask_x[-1]+1])
  w = 0.5*(w[:-1, mask_y[0]:mask_y[-1]+1, mask_x[0]:mask_x[-1]+1] + 
           w[1:, mask_y[0]:mask_y[-1]+1, mask_x[0]:mask_x[-1]+1])
  return ( {'x': x_centers, 'y': y_centers, 'z': z_centers, 'values': u}, 
           {'x': x_centers, 'y': y_centers, 'z': z_centers, 'values': v}, 
           {'x': x_centers, 'y': y_centers, 'z': z_centers, 'values': w} )


def main():
  """Converts PETSc output to VTK format for 3D case."""
  # parse command-line
  parameters = read_inputs()
  print ('[case-directory] %s' % parameters.case_directory)
  print ('[variables] %s' % parameters.variables)

  # list of time-steps to post-process
  time_steps = ioPetIBM.get_time_steps(parameters.case_directory, 
                                       parameters.time_steps)

  # read mesh grid
  x, y, z = ioPetIBM.read_grid(parameters.case_directory, 
                               bottom_left=parameters.bottom_left,
                               top_right=parameters.top_right)

  for time_step in time_steps:
    if 'velocity' in parameters.variables:
      u, v, w = ioPetIBM.read_velocity(parameters.case_directory, time_step, [x, y, z],
                                       periodic=parameters.periodic)
      # need to get values at cell-centers, not staggered arrangement
      u, v, w = interpolate_cell_centers(u, v, w, x, y, z)
      return
      ioPetIBM.write_velocity_vtk3d(u, v, w, parameters.case_directory, time_step)
    if 'pressure' in parameters.variables:
      p = ioPetIBM.read_pressure(parameters.case_directory, time_step, [x, y, z])
      ioPetIBM.write_pressure_vtk3d(p, parameters.case_directory, time_step)

  print('\n[{}] DONE'.format(os.path.basename(__file__)))


if __name__ == "__main__":
  main()