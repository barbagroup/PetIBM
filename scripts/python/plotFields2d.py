#!/usr/bin/env python

# file: plotFields2d.py
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: Plots the 2D vorticity, pressure and velocity fields.


import sys
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
  return parser.parse_args()


def vorticity(u, v, dx, dy):
  """Computes the vorticity field for a two-dimensional simulation.

  Arguments
  ---------
  u, v -- x- and y- velocity fields
  dx, dy -- cell-widths in the x- and y- direction
  """
  return ( (v[:,1:]-v[:,:-1]) 
           / numpy.outer(numpy.ones(dy.size-1), 0.5*(dx[:-1]+dx[1:])) 
           - (u[1:,:]-u[:-1,:]) 
           / numpy.outer(0.5*(dy[:-1]+dy[1:]), numpy.ones(dx.size-1)) )


def plot_field(x, y, variable, variable_range, variable_name,
              bottom_left, top_right, image_path):
  """Plots and saves the variable field.

  Arguments
  ---------
  x, y -- node coordinates in the x- and y- directions
  variable -- field to plot
  variable_range -- contour values to plot
  variable_name -- name of the variable
  bottom_left, top_right -- limits of the plot
  image_path -- path of the image to save
  """
  pyplot.style.use('{}/style_PetIBM.mplstyle'.format(os.path.dirname(__file__)))
  fig, ax = pyplot.subplots()
  pyplot.xlabel('$x$')
  pyplot.ylabel('$y$')
  levels = numpy.linspace(variable_range[0], 
                          variable_range[1], 
                          variable_range[2])
  X, Y = numpy.meshgrid(x, y)
  color_map = {'pressure': cm.jet, 'vorticity': cm.RdBu_r,
               'u-velocity': cm.RdBu_r, 'v-velocity': cm.RdBu_r}
  print X.shape
  print Y.shape
  print variable.shape
  cont = ax.contourf(X, Y, variable, 
                     levels=levels, extend='both', 
                     cmap=color_map[variable_name])
  cont_bar = fig.colorbar(cont, label=variable_name, fraction=0.046, pad=0.04)
  ax.axis([bottom_left[0], top_right[0], bottom_left[1], top_right[1]])
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

  # get time-steps to plot
  if any(parameters.time_steps):
    time_steps = range(parameters.time_steps[0], 
                       parameters.time_steps[1]+1, 
                       parameters.time_steps[2])
  else:
    time_steps = sorted(int(folder) 
                        for folder in os.listdir(parameters.case_directory)
                        if folder[0] == '0')
 
  # create directory where images will be saved
  images_directory = '{}/images'.format(parameters.case_directory)
  print('[images directory] {}'.format(images_directory))
  if not os.path.isdir(images_directory):
    os.makedirs(images_directory)

  # read the grid nodes
  [x, y], [dx, dy] = ioPetIBM.read_grid(parameters.case_directory,
                                        bottom_left=parameters.bottom_left,
                                        top_right=parameters.top_right)
  nx, ny = dx.size, dy.size

  print nx, ny
  print parameters.bottom_left
  print parameters.top_right
  print x, y
  print dx, dy
  return

  for time_step in time_steps:
    # velocity fields
    u, v = ioPetIBM.read_velocity(parameters.case_directory, time_step, [dx, dy],
                                  periodic=parameters.periodic)
    # u-velocity plot
    nxu = (nx if 'x' in parameters.periodic else nx-1)
    image_path = '{}/uVelocity{:0>7}.png'.format(images_directory, time_step)
    plot_field(x[1:nxu+1], 0.5*(y[:-1]+y[1:]), u, 
               parameters.u_range, 'u-velocity', 
               parameters.bottom_left, parameters.top_right, image_path)
    # v-velocity plot
    nyv = (ny if 'y' in parameters.periodic else ny-1)
    image_path = '{}/vVelocity{:0>7}.png'.format(images_directory, time_step)
    plot_field(0.5*(x[:-1]+x[1:]), y[1:nyv+1], v, 
               parameters.v_range, 'v-velocity', 
               parameters.bottom_left, parameters.top_right, image_path)
    # vorticity field
    w = vorticity(u, v, dx, dy)
    # vorticity plot
    image_path = '{}/vorticity{:0>7}.png'.format(images_directory, time_step)
    plot_field(x[1:-1], y[1:-1], w, 
               parameters.vorticity_range, 'vorticity', 
               parameters.bottom_left, parameters.top_right, image_path)
    # pressure field
    p = ioPetIBM.read_pressure(parameters.case_directory, time_step, [nx, ny])
    # pressure plot
    image_path = '{}/pressure{:0>7}.png'.format(images_directory, time_step)
    plot_field(0.5*(x[:-1]+x[1:]), 0.5*(y[:-1]+y[1:]), p, 
               parameters.pressure_range, 'pressure', 
               parameters.bottom_left, parameters.top_right, image_path)

    print('\n[{}] DONE'.format(os.path.basename(__file__)))


if __name__ == '__main__':
  main()