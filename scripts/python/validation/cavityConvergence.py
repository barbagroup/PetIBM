#!/usr/bin/env python

# file: cavityConvergence.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Plots the grid-convergence for the lid-driven cavity case.


import os
import sys
import argparse
import math

import numpy
from matplotlib import pyplot

sys.path.append('{}/scripts/python'.format(os.environ['PETIBM_DIR']))
import ioPetIBM


def read_inputs():
  """Parses the command-line."""
  # create parser
  parser = argparse.ArgumentParser(description='Convergence for the '
                                               'lid-driven cavity case',
                        formatter_class= argparse.ArgumentDefaultsHelpFormatter)
  # fill parser with arguments
  parser.add_argument('--directory', dest='directory', type=str,
                      default=os.getcwd(),
                      help='directory containing the simulation folders')
  parser.add_argument('--time-step', '-ts', dest='time_step', type=int, 
                      default=500,
                      help='time-step at which the solution will be read')
  parser.add_argument('--no-save', dest='save', action='store_false',
                      help='does not save the figure')
  parser.add_argument('--output', '-o', dest='output', type=str, 
                      default='grid_convergence',
                      help='name of the .png file saved')
  parser.add_argument('--no-show', dest='show', action='store_false',
                      help='does not display the figure')
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


def main():
  """Plots the grid convergence for the lid-driven cavity case."""
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
    grid = ioPetIBM.read_grid(case['directory'])
    cases[i]['grid-spacing'] = (grid[0][-1]-grid[0][0])/(grid[0].size-1)
    # read velocity components
    cases[i]['u'], cases[i]['v'] = ioPetIBM.read_velocity(case['directory'], 
                                                          args.time_step, 
                                                          grid)
    # pressure
    cases[i]['p'] = ioPetIBM.read_pressure(case['directory'], 
                                           args.time_step, 
                                           grid)

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

  # grid convergence, comparison with finest grid
  fine = cases[-1]
  for i, case in enumerate(cases[:-1]):
    u_fine = restriction(fine['u'], case['u'])
    cases[i]['u'].error = (l2_norm(case['u'].values-u_fine.values)
                           / l2_norm(u_fine.values))
    v_fine = restriction(fine['v'], case['v'])
    cases[i]['v'].error = (l2_norm(case['v'].values-v_fine.values)
                           / l2_norm(v_fine.values))
    p_fine = restriction(fine['p'], case['p'])
    cases[i]['p'].error = (l2_norm(case['p'].values-p_fine.values)
                           / l2_norm(p_fine.values))

  if args.save or args.show:
    print('\nPlot the grid convergence ...')
    pyplot.style.use('{}/scripts/python/style/'
                     'style_PetIBM.mplstyle'.format(os.environ['PETIBM_DIR']))
    pyplot.xlabel('grid-spacing')
    pyplot.ylabel('$L_2$-norm error')
    # plot errors in u-velocity
    pyplot.plot([case['grid-spacing'] for case in cases[:-1]], 
                [case['u'].error for case in cases[:-1]], 
                label='u-velocity', marker='o')
    # plot errors in v-velocity
    pyplot.plot([case['grid-spacing'] for case in cases[:-1]], 
                [case['v'].error for case in cases[:-1]],
                label='v-velocity', marker='o')
    # plot errors in pressure
    pyplot.plot([case['grid-spacing'] for case in cases[:-1]], 
                [case['p'].error for case in cases[:-1]], 
                label='pressure', marker='o')
    h = numpy.linspace(cases[0]['grid-spacing'], cases[-1]['grid-spacing'], 101)
    # plot convergence-gauge for 1st- and 2nd- order
    pyplot.plot(h, h, label='$1^{st}$-order convergence', color='k')
    pyplot.plot(h, h**2, label='$2^{nd}$-order convergence', 
                color='k', linestyle='--')
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