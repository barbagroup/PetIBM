"""
Cleans a PetIBM simulation (deleting the numerical solution).
"""

import os
import shutil
import argparse


def parse_command_line():
  """
  Parses the command-line.

  Returns
  -------
  args: namespace
    The arguments parsed from the command-line.
  """
  formatter_class = argparse.ArgumentDefaultsHelpFormatter
  parser = argparse.ArgumentParser(description='Cleans a PetIBM simulation',
                                   formatter_class=formatter_class)
  parser.add_argument('--directory',
                      dest='directory',
                      type=str,
                      default=os.getcwd(),
                      help='directory of the PetIBM simulation')
  parser.add_argument('--no-images',
                      dest='images',
                      action='store_false',
                      help='do not remove the images folder')
  parser.add_argument('--no-data',
                      dest='data',
                      action='store_false',
                      help='do not remove the data folder')
  parser.add_argument('--no-grid',
                      dest='grid',
                      action='store_false',
                      help='do not remove the file grid.dat')
  parser.add_argument('--no-solutions',
                      dest='solutions',
                      action='store_false',
                      help='do not remove the numerical solution folders')
  parser.add_argument('--no-forces',
                      dest='forces',
                      action='store_false',
                      help='do not remove the forces.txt')
  parser.add_argument('--no-vtk',
                      dest='vtk_files',
                      action='store_false',
                      help='do not remove .vtk_files folder')
  parser.add_argument('--no-iters',
                      dest='iters',
                      action='store_false',
                      help='do not remove the file iterationCounts.txt')
  parser.add_argument('--no-tensors',
                      dest='tensors',
                      action='store_false',
                      help='do not remove the folder containing the tensors')
  parser.set_defaults(images=True, data=True, grid=True, solutions=True,
                      forces=True, vtk_files=True, iters=True, tensors=True)
  return parser.parse_args()


def main(args):
  """
  Cleans a PetIBM simulation.

  Parameters
  ----------
  args: namespace
    Arguments parsed from the command-line.
  """
  def remove_file(path):
    if os.path.isfile(path):
      os.remove(path)

  def remove_folder(path):
    if os.path.isdir(path):
      shutil.rmtree(path)

  if args.images:
    remove_folder(os.path.join(args.directory, 'images'))
  if args.data:
    remove_folder(os.path.join(args.directory, 'data'))
  if args.grid:
    remove_file(os.path.join(args.directory, 'grid.dat'))
    remove_file(os.path.join(args.directory, 'grid.txt'))
  if args.solutions:
    folders = [os.path.join(args.directory, folder)
               for folder in os.listdir(args.directory)
               if folder.startswith('0')]
    for folder in folders:
      remove_folder(folder)
  if args.forces:
    remove_file(os.path.join(args.directory, 'forces.txt'))
  if args.vtk_files:
    remove_folder(os.path.join(args.directory, 'vtk_files'))
  if args.iters:
    remove_file(os.path.join(args.directory, 'iterationCounts.txt'))
  if args.tensors:
    remove_folder(os.path.join(args.directory, 'tensors'))


if __name__ == '__main__':
  print('\n[{}] START\n'.format(os.path.basename(__file__)))
  args = parse_command_line()
  main(args)
  print('\n[{}] END\n'.format(os.path.basename(__file__)))
