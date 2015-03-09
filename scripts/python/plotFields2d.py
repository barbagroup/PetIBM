#!/usr/bin/env python

# file: plotFields2d.py
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: Plots the 2D vorticity, pressure and velocity fields.


import sys
import os
import argparse

import numpy
from matplotlib import pyplot
sys.path.append(os.path.join(os.environ['PETSC_DIR'], 'bin', 'pythonscripts'))
import PetscBinaryIO


def read_inputs():
  """Parses the command-line."""
  # create parser
  parser = argparse.ArgumentParser(description='Plots the 2D vorticity, '
                                               'pressure and velocity fields',
                        formatter_class= argparse.ArgumentDefaultsHelpFormatter)
  # fill parser with arguments
  parser.add_argument('--case', dest='case_directory', type=str, 
                      default=os.getcwd(), help='directory of the simulation')
  parser.add_argument('--vorticity-limits', '-w', dest='vorticity_limits', 
                      type=float, nargs='+', default=[-5.0, 5.0, 1.0],
                      help='range of vorticity iso-surfaces (min, max, stride)')
  parser.add_argument('--velocity-limits', '-u', dest='velocity_limits', 
                      type=float, nargs='+', default=[-1.0, 1.0, 0.1],
                      help='range of velocity iso-surfaces (min, max, stride)')
  parser.add_argument('--pressure-limits', '-p', dest='pressure_limits', 
                      type=float, nargs='+', default=[-1.0, 1.0, 0.1],
                      help='range of pressure iso-surfaces (min, max, stride)')
  parser.add_argument('--bottom-left', '-bl', dest='bottom_left', type=float,
                      nargs='+', default=[float('-inf'), float('-inf')],
                      help='coordinates of the bottom-left corner')
  parser.add_argument('--top-right', '-tr', dest='top_right', type=float,
                      nargs='+', default=[float('inf'), float('inf')],
                      help='coordinates of the top-right corner')
  parser.add_argument('--time-steps', '-t', dest='time_steps', type=int,
                      nargs='+', default=[None, None, None],
                      help='time-steps to plot (initial, final, increment)')
  parser.add_argument('--periodic', dest='periodic', type=str, nargs='+',
                      default=[], help='direction(s) (x and/or y) with '
                                       'periodic boundary conditions')
  return parser.parse_args()

def plotField(X, Y, variable, variable_range, variable_name,
              bottom_left, top_right, image_path):
  """Plots and saves the variable field.

  Arguments
  ---------
  X, Y -- mesh grid
  variable -- field to plot
  variable_range -- contour values to plot
  variable_name -- name of the variable
  bottom_left, top_right -- limits of the plot
  image_path -- path of the image to save
  """
  pyplot.figure()
  pyplot.xlabel(r'$x$', fontsize=18)
  pyplot.ylabel(r'$y$', fontsize=18)
  levels = numpy.arange(variable_range[0], 
                        variable_range[1]+variable_range[2], 
                        variable_range[2])
  if variable_name == 'vorticity':
    # remove iso-surface where vorticity is ~zero
    levels = numpy.delete(levels, numpy.where(numpy.abs(levels) < 1.0E-06)[0])
  cont = pyplot.contour(X, Y, variable, levels=levels)
  cont_bar = pyplot.colorbar(cont)
  cont_bar.set_label(variable_name)
  pyplot.axis([bottom_left[0], top_right[0], bottom_left[1], top_right[1]])
  pyplot.savefig(image_path)
  pyplot.clf()
  pyplot.close()

def main():
  """Plots the 2D vorticity fields at certain time-steps."""
  # parse the command-line
  args = read_inputs()
  print '[case directory] %s' % args.case_directory

  # get time-steps to plot
  if any(args.time_steps):
    time_steps = range(args.time_steps[0], 
                       args.time_steps[1]+1, 
                       args.time_steps[2])
  else:
    time_steps = sorted(int(folder) for folder in os.listdir(args.case_directory)
                                    if folder[0] == '0')
 
  # create directory where images will be saved
  images_directory = '%s/images' % args.case_directory
  print ('[images directory] %s' % images_directory)
  if not os.path.isdir(images_directory):
    os.makedirs(images_directory)

  # read the grid nodes
  with open('%s/grid.txt' % args.case_directory, 'r') as infile:
    nx, ny = [int(v) for v in infile.readline().strip().split()]
    grid = numpy.loadtxt(infile, dtype=float)
  x, y = grid[:nx+1], grid[nx+1:]

  # mask coordinates outside the plot limits
  args.bottom_left = max(args.bottom_left, [x[0], y[0]])
  args.top_right = min(args.top_right, [x[-1], y[-1]])
  mask_x = numpy.where(numpy.logical_and(x >= args.bottom_left[0],
                                         x <= args.top_right[0]))[0]
  mask_y = numpy.where(numpy.logical_and(y >= args.bottom_left[1],
                                         y <= args.top_right[1]))[0]
  x, y = x[mask_x], y[mask_y]

  # calculate cell-widths
  dx, dy = x[1:]-x[:-1], y[1:]-y[:-1]

  # coordinates of velocities, vorticity and pressure nodes
  xu, yu = x[1:-1], 0.5*(y[:-1]+y[1:])
  xv, yv = 0.5*(y[:-1]+y[1:]), y[1:-1]
  xw, yw = x[1:-1], y[1:-1]
  xp, yp = 0.5*(x[:-1]+x[1:]), 0.5*(y[:-1]+y[1:])

  # number of velocity nodes in each direction (depends on type of bc)
  nxu, nyu = (nx if 'x' in args.periodic else nx-1), ny
  nxv, nyv = nx, (ny if 'y' in args.periodic else ny-1)

  for time_step in time_steps:
    print 'generating plots at time-step %d ...' % time_step
    
    # calculate u-velocity field
    qx = PetscBinaryIO.PetscBinaryIO().readBinaryFile('%s/%07d/qx.dat' 
                                                      % (args.case_directory,
                                                         time_step))[0]
    qx = qx.reshape((nyu, nxu))
    qx = numpy.array([qx[j][mask_x[:-2]] for j in mask_y[:-1]])
    u = qx / numpy.outer(dy, numpy.ones(dx.size-1))
    # generate u-velocity plot
    X, Y = numpy.meshgrid(xu, yu)
    image_path = '{}/uVelocity{:0>7}.png'.format(images_directory, time_step)
    plotField(X, Y, u, args.velocity_limits, 'u-velocity', 
              args.bottom_left, args.top_right, image_path)
    
    # calculate v-velocity field
    qy = PetscBinaryIO.PetscBinaryIO().readBinaryFile('%s/%07d/qy.dat' 
                                                      % (args.case_directory,
                                                         time_step))[0]
    qy = qy.reshape((nyv, nxv))
    qy = numpy.array([qy[j][mask_x[:-1]] for j in mask_y[:-2]])
    v = qy / numpy.outer(numpy.ones(dy.size-1), dx)
    # generate v-velocity plot
    X, Y = numpy.meshgrid(xv, yv)
    image_path = '{}/vVelocity{:0>7}.png'.format(images_directory, time_step)
    plotField(X, Y, v, args.velocity_limits, 'v-velocity', 
              args.bottom_left, args.top_right, image_path)

    # calculate vorticity field
    w = ( (v[:,1:]-v[:,:-1]) 
        / numpy.outer(numpy.ones(dy.size-1), 0.5*(dx[:-1]+dx[1:])) 
        - (u[1:,:]-u[:-1,:]) 
        / numpy.outer(0.5*(dy[:-1]+dy[1:]), numpy.ones(dx.size-1)) )
    # generate vorticity plot
    X, Y = numpy.meshgrid(xw, yw)
    image_path = '{}/vorticity{:0>7}.png'.format(images_directory, time_step)
    plotField(X, Y, w, args.vorticity_limits, 'vorticity', 
              args.bottom_left, args.top_right, image_path)

    # calculate pressure field
    p = PetscBinaryIO.PetscBinaryIO().readBinaryFile('%s/%07d/phi.dat' 
                                                      % (args.case_directory,
                                                         time_step))[0]
    p = p.reshape((ny, nx))
    p = numpy.array([p[j][mask_x[:-1]] for j in mask_y[:-1]])
    # generate pressure plot
    X, Y = numpy.meshgrid(xp, yp)
    image_path = '{}/pressure{:0>7}.png'.format(images_directory, time_step)
    plotField(X, Y, p, args.pressure_limits, 'pressure', 
              args.bottom_left, args.top_right, image_path)


if __name__ == '__main__':
  main()