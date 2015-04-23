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
  parser.add_argument('--dt', '-dt', dest='dt', type=float, default=5.0E-04,
                      help='time-increment of the simulations')
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


def l2_norm(field):
  """Computes the L2-norm of a 2d array

  Parameters
  ----------
  field: 2D Numpy array
    The numerical solution.

  Returns
  -------
  l2: float
    The L2-norm.
  """
  j_start, j_end, j_stride = 0, field.shape[0]+1, 1
  i_start, i_end, i_stride = 0, field.shape[1]+1, 1
  return numpy.linalg.norm(field[j_start:j_end:j_stride, i_start:i_end:i_stride])


def compute_order(ratio, coarse, medium, fine):
  """Computes the observed order of convergence 
  using the solution on three grids.

  Parameters
  ----------
  ratio: float
    Grid-refinement ratio.
  coarse, medium, fine: Numpy array
    Solutions on three consecutive grids restricted on the coarsest grid.

  Returns
  -------
  alpha: float
    The observed order of convergence.
  """
  assert coarse.shape == medium.shape and coarse.shape == fine.shape
  return ( math.log(l2_norm(medium-coarse)/l2_norm(fine-medium))
           / math.log(ratio) )


def restriction(fine, coarse):
  """Restriction of the solution from a fine grid onto a coarse grid.

  Parameters
  ----------
  fine, coarse: ioPetIBM.Field
    Fine and coarse numerical solutions.

  Returns
  -------
  fine_on_coarse: ioPetIBM.Field
    The solution on the fine grid restricted to the coarse grid.
  """
  def intersection(a, b, tolerance=1.0E-06):
    return numpy.any(numpy.abs(a-b[:, numpy.newaxis]) <= tolerance, axis=0)
  mask_x = intersection(fine.x, coarse.x)
  mask_y = intersection(fine.y, coarse.y)

  fine_on_coarse = ioPetIBM.Field(x=fine.x[mask_x], y=fine.y[mask_y],
                                  values=numpy.array([fine.values[j][mask_x]
                                                      for j in xrange(fine.y.size)
                                                      if mask_y[j]]))
  assert numpy.allclose(coarse.x, fine_on_coarse.x, rtol=1.0E-04)
  assert numpy.allclose(coarse.y, fine_on_coarse.y, rtol=1.0E-04)
  assert coarse.values.shape == fine_on_coarse.values.shape
  return fine_on_coarse


def taylor_green_vortex(x, y, 
                        x_start=0.0, x_end=1.0, y_start=0.0, y_end=1.0, 
                        V=1.0, time=0.0, Re=100.0):
  """Computes the analytical solution of the 2D Taylor-Green vortex.

  Parameters
  ----------
  x, y: Numpy array
    Coordinates in the x- and y- directions.
  x_start, x_end, y_start, y_end: float
    Limits of the physical domain; default: 0.0, 1.0, 0.0, 1.0.
  V: float
    Amplitude of the sinusoidal velocity field; default: 1.0.
  time: float
    Time at which the solution is computed; default: 0.0.
  Re: float
    Reynolds number of the flow; default: 100.0.

  Returns
  -------
  u, v, p, w: Numpy array
    Analytical solution (velocities, pressure and vorticity).
  """
  X1, X2 = 0.0, 2.0*math.pi
  x = X1 + (X2-X1)*(x-x_start)/(x_end-x_start)
  y = X1 + (X2-X1)*(y-y_start)/(y_end-y_start)
  X, Y = numpy.meshgrid(x, y)
  # u-velocity
  u = -V*numpy.cos(X)*numpy.sin(Y)*math.exp(-2.0*(2.0*math.pi)**2*time/Re)
  # v-velocity
  v = +V*numpy.sin(X)*numpy.cos(Y)*math.exp(-2.0*(2.0*math.pi)**2*time/Re)
  # pressure
  p = -0.25*(numpy.cos(2.0*X)+numpy.cos(2.0*Y))*math.exp(-4.0*(2.0*math.pi)**2*time/Re)
  # z-vorticity
  w = 2.0*numpy.sin(X)*numpy.sin(Y)*math.exp(-2.0*(2.0*math.pi)**2*time/Re)
  return u, v, p, w


def plot_field(x, y, u, name, image_path):
  """Plots the two-dimensional field.

  Parameters
  ----------
  x, y: Numpy array
    x- and y- coordinates.
  u: Numpy array
    Field to plot.
  name: str
    Description of the field variable.
  image_path: str
    Path of the file to save.
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
  args = read_inputs()

  # initialization
  simulations = sorted(int(directory) 
                       for directory in os.listdir(args.directory)
                       if os.path.isdir('/'.join([args.directory, directory])))
  cases = numpy.empty(len(simulations), dtype=dict) 
  for i, case in enumerate(cases):
    cases[i] = {'directory': '{}/{}'.format(args.directory, simulations[i]),
                'grid-size': '{0}x{0}'.format(simulations[i])}

  for i, case in enumerate(cases):
    print('\n[case] grid-size: {}'.format(case['grid-size']))
    # read mesh grid
    x, y = ioPetIBM.read_grid(case['directory'])
    cases[i]['grid-spacing'] = (x[-1]-x[0])/(x.size-1)
    # read velocity and pressure fields
    cases[i]['u'], cases[i]['v'] = ioPetIBM.read_velocity(case['directory'], 
                                                          args.time_step, 
                                                          [x, y], 
                                                          periodic=['x', 'y'])
    cases[i]['p'] = ioPetIBM.read_pressure(case['directory'], args.time_step, [x, y])
    # compute analytical solution
    cases[i]['u'].exact, _, _, _ = taylor_green_vortex(case['u'].x, 
                                                       case['u'].y, 
                                                       V=args.amplitude, 
                                                       time=args.time_step*args.dt, 
                                                       Re=args.Re)
    _, cases[i]['v'].exact, _, _ = taylor_green_vortex(case['v'].x, 
                                                       case['v'].y, 
                                                       V=args.amplitude, 
                                                       time=args.time_step*args.dt, 
                                                       Re=args.Re)
    _, _, cases[i]['p'].exact, _ = taylor_green_vortex(case['p'].x, 
                                                       case['p'].y, 
                                                       V=args.amplitude, 
                                                       time=args.time_step*args.dt, 
                                                       Re=args.Re)
    # compute L2-norm error
    for field in ['u', 'v', 'p']:
      cases[i][field].error = (l2_norm(case[field].values-case[field].exact)
                               / l2_norm(case[field].exact))

    if args.plot:
      print('\nPlot the field difference between numerical and analytical ...')
      # create directory where images will be saved
      images_directory = '{}/images/differences'.format(case['directory'])
      if not os.path.isdir(images_directory):
        os.makedirs(images_directory)
      # load default style
      pyplot.style.use('{}/scripts/python/style/'
                       'style_PetIBM.mplstyle'.format(os.environ['PETIBM_DIR']))
      # set parameters of the plots
      cases[i]['u'].label = 'u-velocity'
      cases[i]['u'].file_name = 'uVelocity'
      cases[i]['v'].label = 'v-velocity'
      cases[i]['v'].file_name = 'vVelocity'
      cases[i]['p'].label = 'pressure'
      cases[i]['p'].file_name = 'pressure'
      # plot velocity fields and pressure field
      for field in ['u', 'v', 'p']:
        image_path = '{}/{}{:0>7}_numerical.png'.format(images_directory,
                                                        case[field].file_name, 
                                                        args.time_step)
        plot_field(case[field].x, case[field].y, case[field].values,
                   case[field].label, image_path)
        image_path = '{}/{}{:0>7}_analytical.png'.format(images_directory,
                                                         case[field].file_name, 
                                                         args.time_step)
        plot_field(case[field].x, case[field].y, case[field].exact,
                   case[field].label, image_path)
        image_path = '{}/{}{:0>7}_difference.png'.format(images_directory,
                                                         case[field].file_name, 
                                                         args.time_step)
        plot_field(case[field].x, case[field].y, case[field].values-case[field].exact,
                   case[field].label, image_path)

  print('\nObserved order of convergence:')
  last_three = True
  coarse, medium, fine = cases[-3:] if last_three else cases[:3]
  ratio = coarse['grid-spacing']/medium['grid-spacing']
  alpha = {'u': compute_order(ratio,
                              coarse['u'].values,
                              restriction(medium['u'], coarse['u']).values,
                              restriction(fine['u'], coarse['u']).values),
           'v': compute_order(ratio,
                              coarse['v'].values,
                              restriction(medium['v'], coarse['v']).values,
                              restriction(fine['v'], coarse['v']).values),
           'p': compute_order(ratio,
                              coarse['p'].values,
                              restriction(medium['p'], coarse['p']).values,
                              restriction(fine['p'], coarse['p']).values)}
  print('\tu: {}'.format(alpha['u']))
  print('\tv: {}'.format(alpha['v']))
  print('\tp: {}'.format(alpha['p']))
  # write orders of convergence into file
  file_path = '{}/orders_of_convergence.dat'.format(args.directory)
  with open(file_path, 'w') as outfile:
    outfile.write('u: {}\n'.format(alpha['u']))
    outfile.write('v: {}\n'.format(alpha['v']))
    outfile.write('p: {}\n'.format(alpha['p']))

  if args.save or args.show:
    print('\nPlot the grid convergence ...')
    pyplot.style.use('{}/scripts/python/style/'
                     'style_PetIBM.mplstyle'.format(os.environ['PETIBM_DIR']))
    pyplot.xlabel('grid-spacing')
    pyplot.ylabel('$L_2$-norm error')
    # plot errors in u-velocity
    pyplot.plot([case['grid-spacing'] for case in cases],
                [case['u'].error for case in cases],
                label='u-velocity', marker='o')
    # plot errors in v-velocity
    pyplot.plot([case['grid-spacing'] for case in cases],
                [case['v'].error for case in cases],
                label='v-velocity', marker='o')
    # plot errors in pressure
    pyplot.plot([case['grid-spacing'] for case in cases],
                [case['p'].error for case in cases],
                label='pressure', marker='o')
    # plot convergence-gauge for 1st- and 2nd- orders
    h = numpy.linspace(cases[0]['grid-spacing'], cases[-1]['grid-spacing'], 101)
    pyplot.plot(h, h, label='$1^{st}$-order convergence', color='k')
    pyplot.plot(h, h**2, label='$2^{nd}$-order convergence', color='k', linestyle='--')
    pyplot.legend()
    pyplot.xscale('log')
    pyplot.yscale('log')
    if args.save:
      pyplot.savefig('{}/{}.png'.format(args.directory, args.output))
    if args.show:
      pyplot.show()


if __name__ == '__main__':
  print('\n[{}] START\n'.format(os.path.basename(__file__)))
  main()
  print('\n[{}] END\n'.format(os.path.basename(__file__)))