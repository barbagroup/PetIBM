#!/usr/bin/env python

# file: cavity.py
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: Plots the centerline velocities for lid-driven cavity flow case
#              and compare with experimental data from Ghia et al. (1982).


import os
import sys
import argparse

import numpy
from matplotlib import pyplot
sys.path.append(os.path.join(os.environ['PETSC_DIR'], 'bin', 'pythonscripts'))
import PetscBinaryIO


def read_inputs():
  """Parses the command-line."""
  # create parser
  parser = argparse.ArgumentParser(description='Plots centerline velocities '
                                               'for lid-driven cavity flow',
                        formatter_class= argparse.ArgumentDefaultsHelpFormatter)
  # fill parser with arguments
  parser.add_argument('--case', dest='case_directory', type=str, 
                      default=os.getcwd(),
                      help='directory of the case containing the solution')
  parser.add_argument('--Re', '-Re', dest='Re', type=str, default='100',
                      help='Reynolds number of the flow')
  parser.add_argument('--time-step', '-t', dest='time_step', type=int, 
                      nargs='+', default=None,
                      help='time-step to plot')
  parser.add_argument('--periodic', dest='periodic', type=str, nargs='+',
                      default=[], help='direction(s) (x and/or y) with '
                                       'periodic boundary conditions')
  return parser.parse_args()


def main():
  """Plots and writes the velocity components at the centerline of the cavity
  and compares with experimental results form Ghia et al. (1982).
  """
  # column indices in file with experimental results
  cols = {'100': {'u': 1, 'v': 7},
          '1000': {'u':2, 'v': 8},
          '3200': {'u':3, 'v': 9},
          '5000': {'u':4, 'v': 10},
          '10000': {'u': 5, 'v': 11}}
  experimental_file = ( '%s/scripts/validation/data/cavity-GGS82.txt' 
                        % os.environ['PETIBM_DIR'] )

  # parse command-line
  args = read_inputs()
  print '[case directory] %s' % args.case_directory
  Re = args.Re

  if not args.time_step:
    args.time_step = sorted(int(folder) for folder in os.listdir(args.case_directory)
                                        if folder[0] == '0')[-1]

  # read the grid nodes
  with open('%s/grid.txt' % args.case_directory, 'r') as infile:
    nx, ny = [int(v) for v in infile.readline().strip().split()]
    grid = numpy.loadtxt(infile, dtype=float)
  x, y = grid[:nx+1], grid[nx+1:]

  if nx%2 == 1 or ny%2 == 1:
    print ( 'Please redo simulation with an even number of grid points '
            'along each direction' )
    sys.exit()

  # create directory where images will be saved
  images_directory = '%s/images' % args.case_directory
  print ('[images directory] %s' % images_directory)
  if not os.path.isdir(images_directory):
    os.makedirs(images_directory)

  # create directory where data will be saved
  data_directory = '%s/data' % args.case_directory
  print ('[data directory] %s' % data_directory)
  if not os.path.isdir(data_directory):
    os.makedirs(data_directory)

  # cell-widths
  dx, dy = x[1:]-x[:-1], y[1:]-y[:-1]

  # velocity nodes
  xu, yu = x[1:-1], 0.5*(y[:-1]+y[1:])
  xv, yv = 0.5*(y[:-1]+y[1:]), y[1:-1]

  # number of velocity nodes in each direction (depends on type of bc)
  nxu, nyu = (nx if 'x' in args.periodic else nx-1), ny
  nxv, nyv = nx, (ny if 'y' in args.periodic else ny-1)

  # calculate u-velocity field
  print 'reading fluxes in x-direction ...'
  qx = PetscBinaryIO.PetscBinaryIO().readBinaryFile('%s/%07d/qx.dat' 
                                                      % (args.case_directory,
                                                         args.time_step))[0]
  qx = qx.reshape((nyu, nxu))[:ny, :nx-1]
  u = qx / numpy.outer(dy, numpy.ones(dx.size-1))

  # calculate v-velocity field
  print 'reading fluxes in y-direction ...'
  qy = PetscBinaryIO.PetscBinaryIO().readBinaryFile('%s/%07d/qy.dat' 
                                                      % (args.case_directory,
                                                         args.time_step))[0]
  qy = qy.reshape((nyv, nxv))[:ny-1, :nx]
  v = qy / numpy.outer(numpy.ones(dy.size-1), dx)

  # plot and write u-velocity along y-centerline
  print 'plotting u-velocity along vertical centerline ...'
  image_path = '%s/uVelocityCenterlineRe%s_%dx%d.png' % (images_directory, Re, 
                                                         nx, ny)
  data_path = '%s/uVelocityCenterlineRe%s_%dx%d.txt' % (data_directory, Re,
                                                        nx, ny)
  pyplot.figure()
  pyplot.grid(True)
  pyplot.xlabel('y', fontsize=16)
  pyplot.ylabel('u-velocity along vertical centerline', fontsize=16)
  pyplot.plot(yu, u[:, nx/2-1], label='PetIBM', 
              color='k', linestyle='-', linewidth=1)
  if Re in list(cols.keys()):
    print '\tcompare with Ghia et al. (1982)'
    with open(experimental_file, 'r') as infile:
      y_exp, u_exp = numpy.loadtxt(infile, dtype=float, 
                                   usecols= (0, cols[Re]['u']), unpack=True)
    pyplot.plot(y_exp, u_exp, label='Ghia et al. (1982)', 
                color='k', linestyle='--', linewidth=1)
  pyplot.legend(loc='best', prop={'size': 16})
  pyplot.title('Lid-driven cavity flow at Re=%s (mesh: %dx%d)' % (Re, nx, ny),
               fontsize=16)
  pyplot.savefig(image_path)
  pyplot.clf()
  pyplot.close()
  print 'writing u-velocity along vertical centerline ...'
  with open(data_path, 'w') as outfile:
    numpy.savetxt(outfile, numpy.c_[yu, u[:, nx/2-1]], fmt='%.6f', delimiter='\t')

  # plot and write v-velocity along x-centerline
  print 'plotting v-velocity along horizontal centerline ...'
  image_path = '%s/vVelocityCenterlineRe%s_%dx%d.png' % (images_directory, Re, 
                                                         nx, ny)
  data_path = '%s/vVelocityCenterlineRe%s_%dx%d.txt' % (data_directory, Re,
                                                        nx, ny)
  pyplot.figure()
  pyplot.grid(True)
  pyplot.xlabel('x', fontsize=16)
  pyplot.ylabel('v-velocity along horizontal centerline', fontsize=16)
  pyplot.plot(xv, v[ny/2-1, :], label='PetIBM', 
              color='k', linestyle='-', linewidth=1)
  if Re in list(cols.keys()):
    print '\tcompare with Ghia et al. (1982)'
    with open(experimental_file, 'r') as infile:
      x_exp, u_exp = numpy.loadtxt(infile, dtype=float, 
                                   usecols= (6, cols[Re]['v']), unpack=True)
    pyplot.plot(x_exp, u_exp, label='Ghia et al. (1982)', 
                color='k', linestyle='--', linewidth=1)
  pyplot.legend(loc='best', prop={'size': 16})
  pyplot.title('Lid-driven cavity flow at Re=%s (mesh: %dx%d)' % (Re, nx, ny),
               fontsize=16)
  pyplot.savefig(image_path)
  pyplot.clf()
  pyplot.close()
  print 'writing v-velocity along horizontal centerline ...'
  with open(data_path, 'w') as outfile:
    numpy.savetxt(outfile, numpy.c_[xv, v[ny/2-1,:]], fmt='%.6f', delimiter='\t')


if __name__ == "__main__":
  main()