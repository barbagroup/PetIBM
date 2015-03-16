#!/usr/bin/env python

# file: cavityConvergence.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Plots the grid-convergence for the lid-driven cavity case.


import os
import sys
import argparse
import math

import numpy

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
                      default=1000,
                      help='time-step at which the solution will be read')
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
  return ( math.log(numpy.linalg.norm(medium-coarse)
                    /numpy.linalg.norm(fine-medium))
           /math.log(ratio) )


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


def main():
  """Plots the grid convergence for the lid-driven cavity case."""
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
    # read velocity components
    cases[i]['u'], cases[i]['v'] = ioPetIBM.read_velocity(case['directory'], 
                                                          parameters.time_step, 
                                                          [x, y])
    # pressure
    cases[i]['p'] = ioPetIBM.read_pressure(case['directory'], 
                                           parameters.time_step, 
                                           [x, y])

  print('\nObserved order of convergence:')
  last_three = True
  coarse, medium, fine = ([cases[-3], cases[-2], cases[-1]] 
                          if last_three 
                          else [cases[0], cases[1], cases[2]])
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

  print('\n[{}] DONE'.format(os.path.basename(__file__)))


if __name__ == '__main__':
  main()