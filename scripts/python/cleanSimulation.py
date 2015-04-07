#!/usr/bin/env/ python

# file: cleanSimulation.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Clean a PetIBM simulation.


import os
import argparse


def read_inputs():
  """Parses the command-line."""
  # create parser
  parser = argparse.ArgumentParser(description='Clean PetIBM case',
                        formatter_class= argparse.ArgumentDefaultsHelpFormatter)
  # fill parser with arguments
  parser.add_argument('--case', dest='case_directory', type=str, 
                      default=os.getcwd(),
                      help='directory of the PetIBM simulation')
  parser.add_argument('--no-images', dest='images', action='store_false',
                      help='does not remove the images folder')
  parser.add_argument('--no-data', dest='data', action='store_false',
                      help='does not remove the data folder')
  parser.add_argument('--no-grid', dest='grid', action='store_false',
                      help='does not remove the grid file')
  parser.add_argument('--no-solutions', dest='solutions', action='store_false',
                      help='does not remove the numrical solution folders')
  parser.add_argument('--no-forces', dest='forces', action='store_false',
                      help='does not remove the forces data file')
  parser.add_argument('--no-vtk', dest='vtk_files', action='store_false',
                      help='does not remove .vtk_files folder')
  parser.add_argument('--no-logs', dest='logs', action='store_false',
                      help='does not remove log files '
                           '(iterationCount, performanceSummary)')
  parser.set_defaults(images=True, data=True, grid=True, solutions=True, 
                      forces=True, vtk_files=True, logs=True)
  return parser.parse_args()


def main():
  """Cleans a PetIBM simulation."""
  # parse command-line
  args = read_inputs()
  # get different paths to delete
  paths = {}
  if args.images:
    paths['images'] = '{}/images'.format(args.case_directory)
  if args.data:
    paths['data'] = '{}/data'.format(args.case_directory)
  if args.grid:
    paths['grid'] = '{}/grid.txt'.format(args.case_directory)
  if args.solutions:
    paths['solutions'] = '{}/0*'.format(args.case_directory)
  if args.forces:
    paths['forces'] = '{}/forces.txt'.format(args.case_directory)
  if args.vtk_files:
    paths['vtk_files'] = '{}/vtk_files'.format(args.case_directory)
  if args.logs:
    paths['logs'] = ('{0}/iterationCount.txt '
                     '{0}/performanceSummary.txt'.format(args.case_directory))
  # delete appropriate files/folders
  print('[case-directory] {}'.format(args.case_directory))
  for key, path in paths.iteritems():
    print('\t-> removing {} ...'.format(key))
    os.system('rm -rf {}'.format(path))


if __name__ == '__main__':
  print('\n[{}] START\n'.format(os.path.basename(__file__)))
  main()
  print('\n[{}] END\n'.format(os.path.basename(__file__)))