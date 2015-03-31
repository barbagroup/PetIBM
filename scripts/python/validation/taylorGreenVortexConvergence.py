#!/usr/bin/env python

# file: taylorGreenVortexConvergence.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Plots the grid-convergence for the Taylor-Green vortex case.


import os
import sys
import argparse
import math

import numpy
from matplotlib import pyplot, cm

sys.path.append('{}/scripts/python'.format(os.environ['PETIBM_DIR']))
import ioPetIBM


def read_inputs():
  """Parses the command-line."""
  # create parser
  parser = argparse.ArgumentParser(description='Convergence for the '
                                               'Taylor-Green vortex case',
                        formatter_class= argparse.ArgumentDefaultsHelpFormatter)
  # fill parser with arguments
  parser.add_argument('--directory', dest='directory', type=str,
                      default=os.getcwd(),
                      help='directory containing the simulation folders')
  parser.add_argument('--Re', '-Re', dest='Re', type=float, default=100.0,
                      help='Reynolds number of the simulation')
  parser.add_argument('--time', '-t', dest='time', type=float, default=0.5,
                      help='time at which the error will be computed')
  parser.add_argument('--time-step', '-ts', dest='time_step', type=int, 
                      default=1000,
                      help='time-step at which the solution will be read')
  parser.add_argument('--amplitude', '-a', dest='amplitude', type=float, 
                      default=1.0, help='amplitude of the Taylor-Green vortex')
  parser.add_argument('--no-save', dest='save', action='store_false',
                      help='does not save the figure')
  parser.add_argument('--output', '-o', dest='output', type=str, 
                      default='grid_convergence',
                      help='name of the .png file saved')
  parser.add_argument('--no-show', dest='show', action='store_false',
                      help='does not display the figure')
  parser.add_argument('--plot', dest='plot', action='store_true',
                      help='plots the field difference between the numerical '
                           'and analytical solutions')
  parser.set_defaults(save=True, show=True)
  # parse command-line
  return parser.parse_args()


def compute_order(ratio, coarse, medium, fine):
  """Computes the observed order of convergence 
  using the solution on three grids.

  Arguments
  ---------
  ratio -- grid-refinement ratio
  coarse, medium, fine -- solutions on three consecutive grids 
                          restricted on the coarsest grid
  """
  return math.log(l2_norm(medium-coarse)/l2_norm(fine-medium)) / math.log(ratio)


def l2_norm(x):
  """Return the discrete L2 norm of x."""
  return math.sqrt(numpy.sum(x**2)/x.size)


def restriction(fine, coarse):
  """Restriction of the solution from a fine grid onto a coarse grid.

  Arguments
  ---------
  fine, coarse -- fine and coarse numerical solutions
  """
  def intersection(a, b, tolerance=1.0E-06):
    return numpy.any(numpy.abs(a-b[:, numpy.newaxis]) <= tolerance, axis=0)
  mask_x = intersection(fine['x'], coarse['x'])
  mask_y = intersection(fine['y'], coarse['y'])
  return {'x': fine['x'][mask_x],
          'y': fine['y'][mask_y],
          'values': numpy.array([fine['values'][j][mask_x] 
                                 for j in xrange(fine['y'].size) if mask_y[j]])}


def taylor_green_vortex(x, y, V=1.0, time=0.0, Re=100.0):
  """Computes the analytical solution of the 2D Taylor-Green vortex.

  Arguments
  ---------
  x, y -- coordinates in the x- and y- directions
  V -- amplitude of the sinusoidal velocity field (default 1.0)
  time -- time at which the solution is computed (default 0.0)
  Re -- Reynolds number of the flow (default 100.0)
  """
  X1, X2 = 0.0, 2.0*math.pi
  x = X1 + (X2-X1)*(x-x[0])/(x[-1]-x[0])
  y = X1 + (X2-X1)*(y-y[0])/(y[-1]-y[0])
  X, Y = numpy.meshgrid(x, y)
  u = -V*numpy.cos(X)*numpy.sin(Y)*math.exp(-2.0*time)
  v = +V*numpy.sin(X)*numpy.cos(Y)*math.exp(-2.0*time)
  p = -0.25*(numpy.cos(2.0*X)+numpy.cos(2.0*Y))*math.exp(-4.0*time)
  w = 2.0*numpy.sin(X)*numpy.sin(Y)*math.exp(-2.0*time)
  return u, v, p, w


def plot_field(x, y, u, name, image_path):
  """Plots the two-dimensional field.

  Arguments
  ---------
  x, y -- x- and y- coordinates
  u -- field to plot
  name -- description of the field variable
  image_path -- path of the file to save
  """
  fig, ax = pyplot.subplots()
  pyplot.xlabel('$x$')
  pyplot.ylabel('$y$')
  X, Y = numpy.meshgrid(x, y)
  levels = numpy.linspace(u.min(), u.max(), 10)
  cont = pyplot.contourf(X, Y, u, 
                         levels=levels, extend='both', cmap=cm.jet)
  cont_bar = fig.colorbar(cont, label='{}'.format(name), 
                          fraction=0.046, pad=0.04)
  ax.axis([x.min(), x.max(), y.min(), y.max()])
  ax.set_aspect('equal')
  pyplot.savefig(image_path)
  pyplot.close()


def main():
  """Plots the grid convergence for the Taylor-Green vortex case."""
  # parse command-line
  parameters = read_inputs()

  # initialization
  simulations = sorted(int(directory) 
                       for directory in os.listdir(parameters.directory)
                       if os.path.isdir('/'.join([parameters.directory, directory])))
  cases = numpy.empty(len(simulations), dtype=dict) 
  for i, case in enumerate(cases):
    cases[i] = {'directory': '{}/{}'.format(parameters.directory, simulations[i]),
                'grid-size': '{0}x{0}'.format(simulations[i])}

  for i, case in enumerate(cases):
    print('\n[case] grid-size: {}'.format(case['grid-size']))
    # read mesh grid
    x, y = ioPetIBM.read_grid(case['directory'])
    cases[i]['grid-spacing'] = (x[-1]-x[0])/(x.size-1)
    # read velocity and pressure fields
    cases[i]['u'], cases[i]['v'] = ioPetIBM.read_velocity(case['directory'], 
                                                          parameters.time_step, 
                                                          [x, y], 
                                                          periodic=['x', 'y'])
    cases[i]['p'] = ioPetIBM.read_pressure(case['directory'], 
                                           parameters.time_step, 
                                           [x, y])
    # compute analytical solution
    cases[i]['u']['analytical'], _, _, _ = taylor_green_vortex(case['u']['x'], 
                                                               case['u']['y'], 
                                                               V=parameters.amplitude, 
                                                               time=parameters.time , 
                                                               Re=parameters.Re)
    _, cases[i]['v']['analytical'], _, _ = taylor_green_vortex(case['v']['x'], 
                                                               case['v']['y'], 
                                                               V=parameters.amplitude, 
                                                               time=parameters.time , 
                                                               Re=parameters.Re)
    _, _, cases[i]['p']['analytical'], _ = taylor_green_vortex(case['p']['x'], 
                                                               case['p']['y'], 
                                                               V=parameters.amplitude, 
                                                               time=parameters.time , 
                                                               Re=parameters.Re)
    # compute L2-norm error
    cases[i]['u']['error'] = l2_norm(case['u']['values']-case['u']['analytical'])
    cases[i]['v']['error'] = l2_norm(case['v']['values']-case['v']['analytical'])
    cases[i]['p']['error'] = l2_norm(case['p']['values']-case['p']['analytical'])
    if parameters.plot:
      print('\nPlot the field difference between numerical and analytical ...')
      # create directory where images will be saved
      images_directory = '{}/images/differences'.format(case['directory'])
      if not os.path.isdir(images_directory):
        os.makedirs(images_directory)
      # load default style
      pyplot.style.use('{}/scripts/python/style/'
                       'style_PetIBM.mplstyle'.format(os.environ['PETIBM_DIR']))
      # u-velocity
      image_path = '{}/uVelocity{:0>7}_numerical.png'.format(images_directory, 
                                                             parameters.time_step)
      plot_field(case['u']['x'], case['u']['y'], case['u']['values'], 
                 'u-velocity', image_path)
      image_path = '{}/uVelocity{:0>7}_analytical.png'.format(images_directory,
                                                              parameters.time_step)
      plot_field(case['u']['x'], case['u']['y'], case['u']['analytical'], 
                 'u-velocity (analytical)', image_path)
      image_path = '{}/uVelocity{:0>7}_difference.png'.format(images_directory, 
                                                              parameters.time_step)
      plot_field(case['u']['x'], case['u']['y'], 
                 case['u']['values']-case['u']['analytical'], 
                 'difference in u-velocity', image_path)
      # v-velocity
      image_path = '{}/vVelocity{:0>7}_numerical.png'.format(images_directory, 
                                                             parameters.time_step)
      plot_field(case['v']['x'], case['v']['y'], case['v']['values'], 
                 'v-velocity', image_path)
      image_path = '{}/vVelocity{:0>7}_analytical.png'.format(images_directory,
                                                              parameters.time_step)
      plot_field(case['v']['x'], case['v']['y'], case['v']['analytical'], 
                 'v-velocity (analytical)', image_path)
      image_path = '{}/vVelocity{:0>7}_difference.png'.format(images_directory, 
                                                              parameters.time_step)
      plot_field(case['v']['x'], case['v']['y'], 
                 case['v']['values']-case['v']['analytical'], 
                 'difference in v-velocity', image_path)
      # pressure
      image_path = '{}/pressure{:0>7}_numerical.png'.format(images_directory, 
                                                            parameters.time_step)
      plot_field(case['p']['x'], case['p']['y'], case['p']['values'], 
                 'pressure', image_path)
      image_path = '{}/pressure{:0>7}_analytical.png'.format(images_directory,
                                                             parameters.time_step)
      plot_field(case['p']['x'], case['p']['y'], case['p']['analytical'], 
                 'pressure (analytical)', image_path)
      image_path = '{}/pressure{:0>7}_difference.png'.format(images_directory, 
                                                             parameters.time_step)
      plot_field(case['p']['x'], case['p']['y'], 
                 case['p']['values']-case['p']['analytical'], 
                 'difference in pressure', image_path)

  print('\nObserved order of convergence:')
  last_three = True
  coarse, medium, fine = cases[-3:] if last_three else cases[:3]
  ratio = coarse['grid-spacing']/medium['grid-spacing']
  alpha = {'u': compute_order(ratio,
                              coarse['u']['values'],
                              restriction(medium['u'], coarse['u'])['values'],
                              restriction(fine['u'], coarse['u'])['values']),
           'v': compute_order(ratio,
                              coarse['v']['values'],
                              restriction(medium['v'], coarse['v'])['values'],
                              restriction(fine['v'], coarse['v'])['values']),
           'p': compute_order(ratio,
                              coarse['p']['values'],
                              restriction(medium['p'], coarse['p'])['values'],
                              restriction(fine['p'], coarse['p'])['values'])}
  print('\tu: {}'.format(alpha['u']))
  print('\tv: {}'.format(alpha['v']))
  print('\tp: {}'.format(alpha['p']))

  if parameters.save or parameters.show:
    print('\nPlot the grid convergence ...')
    pyplot.style.use('{}/scripts/python/style/'
                     'style_PetIBM.mplstyle'.format(os.environ['PETIBM_DIR']))
    pyplot.xlabel('cell-width')
    pyplot.ylabel('$L_2$-norm error')
    pyplot.plot([case['grid-spacing'] for case in cases], 
                [case['u']['error'] for case in cases], 
                label='u-velocity', marker='o')
    pyplot.plot([case['grid-spacing'] for case in cases], 
                [case['v']['error'] for case in cases],
                label='v-velocity', marker='o')
    pyplot.plot([case['grid-spacing'] for case in cases], 
                [case['p']['error'] for case in cases], 
                label='pressure', marker='o')
    h = numpy.linspace(cases[0]['grid-spacing'], cases[-1]['grid-spacing'], 101)
    pyplot.plot(h, h, label='$1^{st}$-order convergence', color='k')
    pyplot.plot(h, h**2, label='$2^{nd}$-order convergence', 
                color='k', linestyle='--')
    pyplot.legend()
    pyplot.xscale('log')
    pyplot.yscale('log')
    if parameters.save:
      pyplot.savefig('{}/{}.png'.format(parameters.directory, parameters.output))
    if parameters.show:
      pyplot.show()

  print('\n[{}] DONE'.format(os.path.basename(__file__)))


if __name__ == '__main__':
  main()