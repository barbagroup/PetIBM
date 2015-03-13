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
                'grid-size': '{0}x{0}'.format(simulations[i]),
                'n': simulations[i]}

  for i, case in enumerate(cases):
    print('\n[case] grid-size: {}'.format(case['grid-size']))
    ratio = case['n']/cases[0]['n']
    # read grid nodes
    [x, y], [dx, dy] = ioPetIBM.readGrid(case['directory'])
    nx, ny = dx.size, dy.size
    cases[i]['grid-spacing'] = (x[-1]-x[0])/nx
    # read velocity components
    u, v = ioPetIBM.readVelocity(case['directory'], parameters.time_step, [dx, dy])
    cases[i]['u'] = u[ratio-1::ratio, ratio-1::ratio]
    cases[i]['v'] = v[ratio-1::ratio, ratio-1::ratio]
    # pressure
    p = ioPetIBM.readPressure(case['directory'], parameters.time_step, [nx, ny])
    cases[i]['p'] = p[ratio-1::ratio, ratio-1::ratio]

  print('Orders of convergence:')
  def compute_alpha(v):
    return ( math.log(numpy.linalg.norm(cases[-2][v]-cases[-3][v])
                      /numpy.linalg.norm(cases[-1][v]-cases[-2][v]))
             /math.log(cases[-1]['n']/cases[-2]['n']) ) 
  alpha = {'u': compute_alpha('u'),
           'v': compute_alpha('v'),
           'p': compute_alpha('p')} 
  print('\tu: {}'.format(alpha['u']))
  print('\tv: {}'.format(alpha['v']))
  print('\tp: {}'.format(alpha['p']))


if __name__ == '__main__':
  main()