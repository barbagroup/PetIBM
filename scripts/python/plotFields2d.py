#!/usr/bin/env python

# file: plotFields2d.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Plots the 2D fields from PetIBM numerical solution.


import os
import argparse

import numpy
from matplotlib import pyplot, cm

import ioPetIBM


def read_inputs():
  """Parses the command-line."""
  # create parser
  parser = argparse.ArgumentParser(description='Plots the 2D vorticity, '
                                               'pressure and velocity fields',
                        formatter_class= argparse.ArgumentDefaultsHelpFormatter)
  # fill parser with arguments
  parser.add_argument('--case', dest='case_directory', type=str, 
                      default=os.getcwd(), help='directory of the simulation')
  parser.add_argument('--vorticity-range', '-wr', dest='vorticity_range', 
                      type=float, nargs='+', default=[-5.0, 5.0, 11],
                      help='vorticity range (min, max, number of levels)')
  parser.add_argument('--u-range', '-ur', dest='u_range', 
                      type=float, nargs='+', default=[-1.0, 1.0, 11],
                      help='u-velocity range (min, max, number of levels)')
  parser.add_argument('--v-range', '-vr', dest='v_range', 
                      type=float, nargs='+', default=[-1.0, 1.0, 11],
                      help='v-velocity range (min, max, number of levels)')
  parser.add_argument('--pressure-range', '-pr', dest='pressure_range', 
                      type=float, nargs='+', default=[-1.0, 1.0, 11],
                      help='pressure range (min, max, number of levels)')
  parser.add_argument('--bottom-left', '-bl', dest='bottom_left', type=float,
                      nargs='+', default=[float('-inf'), float('-inf')],
                      help='coordinates of the bottom-left corner of the view')
  parser.add_argument('--top-right', '-tr', dest='top_right', type=float,
                      nargs='+', default=[float('inf'), float('inf')],
                      help='coordinates of the top-right corner of the view')
  parser.add_argument('--time-steps', '-t', dest='time_steps', type=int,
                      nargs='+', default=[None, None, None],
                      help='time-steps to plot (initial, final, increment)')
  parser.add_argument('--periodic', dest='periodic', type=str, nargs='+',
                      default=[], help='direction(s) (x and/or y) with '
                                       'periodic boundary conditions')
  parser.add_argument('--no-velocity', dest='velocity', action='store_false',
                      help='does not plot the velocity fields')
  parser.add_argument('--no-pressure', dest='pressure', action='store_false',
                      help='does not plot the pressure field')
  parser.add_argument('--no-vorticity', dest='vorticity', action='store_false',
                      help='does not plot the vorticity field')
  parser.set_defaults(velocity=True, pressure=True, vorticity=True)
  # parse command-line
  return parser.parse_args()


def vorticity(u, v):
  """Computes the vorticity field for a two-dimensional simulation.

  Arguments
  ---------
  u, v -- dictionaries with x- and y- velocity fields
  """
  print('\tCompute the vorticity field ...')
  u, xu, yu = u['values'], u['x'], u['y']
  v, xv, yv = v['values'], v['x'], v['y']
  mask_x = numpy.where(numpy.logical_and(xu > xv[0], xu < xv[-1]))[0]
  mask_y = numpy.where(numpy.logical_and(yv > yu[0], yv < yu[-1]))[0]
  xu, yv = xu[mask_x], yv[mask_y]
  # compute vorticity node coordinates
  xw, yw = 0.5*(xv[:-1]+xv[1:]), 0.5*(yu[:-1]+yu[1:])
  # compute vorticity field
  w = ( (v[mask_y, 1:]-v[mask_y, :-1])
        /numpy.outer(numpy.ones(yw.size), xv[1:]-xv[:-1]) 
      - (u[1:, mask_x]-u[:-1, mask_x])
        /numpy.outer(yu[1:]-yu[:-1], numpy.ones(xw.size)) )
  # tests
  assert (yw.size, xw.size) == w.shape
  return {'x': xw, 'y': yw, 'values': w}


def plot_contour(variable, variable_name, variable_range, image_path, 
                 view=[[None, None], [None, None]]):
  """Plots and saves the variable field.

  Arguments
  ---------
  variable -- dictionary with the nodal coordinates and values of the variable
  variable_name -- name of the variable
  variable_range -- contour values to plot
  image_path -- path of the image to save
  """
  print('\tPlot the {} contour ...'.format(variable_name))
  fig, ax = pyplot.subplots()
  pyplot.xlabel('$x$')
  pyplot.ylabel('$y$')
  levels = numpy.linspace(variable_range[0], variable_range[1], variable_range[2])
  X, Y = numpy.meshgrid(variable['x'], variable['y'])
  color_map = {'pressure': cm.jet, 'vorticity': cm.RdBu_r,
               'u-velocity': cm.RdBu_r, 'v-velocity': cm.RdBu_r}
  cont = ax.contourf(X, Y, variable['values'], 
                     levels=levels, extend='both', 
                     cmap=color_map[variable_name])
  cont_bar = fig.colorbar(cont, label=variable_name, fraction=0.046, pad=0.04)
  x_start = max(view[0][0], variable['x'].min())
  x_end = min(view[1][0], variable['x'].max())
  y_start = max(view[0][1], variable['y'].min())
  y_end = min(view[1][1], variable['y'].max())
  ax.axis([x_start, x_end, y_start, y_end])
  ax.set_aspect('equal')
  pyplot.savefig(image_path)
  pyplot.close()


def main():
  """Plots the the velocity, pressure and vorticity fields at saved time-steps
  for a two-dimensional simulation.
  """
  # parse command-line
  parameters = read_inputs()
  print('[case directory] {}'.format(parameters.case_directory))

  time_steps = ioPetIBM.get_time_steps(parameters.case_directory, 
                                       parameters.time_steps)
 
  # create directory where images will be saved
  images_directory = '{}/images'.format(parameters.case_directory)
  print('[images directory] {}'.format(images_directory))
  if not os.path.isdir(images_directory):
    os.makedirs(images_directory)

  # read the grid nodes
  coords = ioPetIBM.read_grid(parameters.case_directory)

  # load default style of matplotlib figures
  pyplot.style.use('{}/style_PetIBM.mplstyle'.format(os.path.dirname(__file__)))

  for time_step in time_steps:
    if parameters.velocity or parameters.vorticity:
      # get velocity fields
      u, v = ioPetIBM.read_velocity(parameters.case_directory, time_step, coords,
                                    periodic=parameters.periodic)
      if parameters.velocity:
        # plot u-velocity field
        image_path = '{}/uVelocity{:0>7}.png'.format(images_directory, time_step)
        plot_contour(u, 'u-velocity', parameters.u_range, image_path,
                     view=[parameters.bottom_left, parameters.top_right]) 
        # plot v-velocity field
        image_path = '{}/vVelocity{:0>7}.png'.format(images_directory, time_step)
        plot_contour(v, 'v-velocity', parameters.v_range, image_path,
                     view=[parameters.bottom_left, parameters.top_right])
      if parameters.vorticity:
        # compute vorticity field
        w = vorticity(u, v)
        # plot vorticity field
        image_path = '{}/vorticity{:0>7}.png'.format(images_directory, time_step)
        plot_contour(w, 'vorticity', parameters.vorticity_range, image_path,
                     view=[parameters.bottom_left, parameters.top_right])
    if parameters.pressure:
      # get pressure field
      p = ioPetIBM.read_pressure(parameters.case_directory, time_step, coords)
      # plot pressure field
      image_path = '{}/pressure{:0>7}.png'.format(images_directory, time_step)
      plot_contour(p, 'pressure', parameters.pressure_range, image_path,
                   view=[parameters.bottom_left, parameters.top_right])

  print('\n[{}] DONE'.format(os.path.basename(__file__)))


if __name__ == '__main__':
  main()