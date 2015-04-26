#!/usr/bin/env python

# file: restartFromSolution.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# brief: Uses the solution from a reference simulation as the initial condition.


import os
import sys
import argparse
import shutil

import numpy

sys.path.append('{}/scripts/python'.format(os.environ['PETIBM_DIR']))
import ioPetIBM


def read_inputs():
  """Parses the command-line."""
  # create parser
  parser = argparse.ArgumentParser(description='Interpolates the solution from '
                                               'a grid to another',
                        formatter_class= argparse.ArgumentDefaultsHelpFormatter)
  # fill parser with arguments
  parser.add_argument('--case', dest='case_directory', type=str, 
                      default=os.getcwd(),
                      help='directory of the new simulation')
  parser.add_argument('--reference', dest='reference_directory', type=str,
                      required=True,
                      help='reference directory where solution will be read')
  parser.add_argument('--time-step', '-n', dest='time_step', type=int, 
                      default=None,
                      help='Time-step solution to be read')
  parser.add_argument('--periodic', '-p', dest='periodic', type=str, nargs='+', 
                      default=[],
                      help='direction(s) (x, y, and/or z) with periodic '
                           'boundary conditions')
  parser.add_argument('--same', dest='same_grid', action='store_true',
                      help='skips the interpolation and copies solution')
  # parse command-line
  return parser.parse_args()

def main():
  """Interpolates the solution from a grid to another."""
  args = read_inputs()

  # reference solution
  reference = {}
  if not args.time_step:
    args.time_step = sorted(int(folder) 
                            for folder in os.listdir(args.reference_directory)
                            if folder[0]=='0')[-1]
  reference['input'] = '{}/{:0>7}'.format(args.reference_directory, 
                                          args.time_step)

  # new simulation
  case = {}
  case['output'] = '{}/0000000'.format(args.case_directory)
  if args.same_grid:
    print('Same grid, data are simply copied ...')
    if os.path.isdir(case['output']):
      shutil.rmtree(case['output'])
    shutil.copytree(reference['input'], case['output'])
    return

  # read reference solution
  reference['grid'] = ioPetIBM.read_grid(args.reference_directory)
  reference['u'], reference['v'], reference['w'] = ioPetIBM.read_velocity(args.reference_directory,
                                                                          args.time_step,
                                                                          reference['grid'],
                                                                          periodic=args.periodic)
  reference['p'] = ioPetIBM.read_pressure(args.reference_directory,
                                          args.time_step,
                                          reference['grid'])


if __name__ == '__main__':
  print('\n[{}] START\n'.format(os.path.basename(__file__)))
  main()
  print('\n[{}] END\n'.format(os.path.basename(__file__)))