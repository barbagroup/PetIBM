#!/usr/bin/env/ python

# file: regressionAnalysis.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Performs regression tests for PetIBM.


import os
import sys
import time
import shutil
import argparse

import numpy

sys.path.append('{}/scripts/python'.format(os.environ['PETIBM_DIR']))
import ioPetIBM


def read_inputs():
  """Parses the command-line."""
  # create parser
  parser = argparse.ArgumentParser(description='Executes a regression-test',
                        formatter_class= argparse.ArgumentDefaultsHelpFormatter)
  # fill parser with arguments
  parser.add_argument('--build', dest='build_directory', type=str, 
                      default=os.getcwd(),
                      help='directory of the PetIBM build')
  parser.add_argument('--save', dest='save', action='store_true', 
                      help='saves the new numerical solutions')
  parser.add_argument('--no-compile', dest='compile', action='store_false',
                      help='skips PetIBM compilation')
  parser.add_argument('--no-run', dest='run', action='store_false',
                      help='does not run the test-cases')
  parser.set_defaults(save=False, compile=True, run=True)
  # parse command-line
  return parser.parse_args()


def define_test_cases(build_directory):
  """Defines the list of test-case to execute.

  Parameters
  ----------
  build_directory: str
    Directory of the PetIBM build.

  Returns
  -------
  tests: list(TestCase)
    List of test-cases.
  """
  tests = []
  tests.append(TestCase(description='2d lid-driven cavity flow at Re=100',
                        directory='{}/examples/2d/lidDrivenCavity/Re100'.format(build_directory),
                        command='make lidDrivenCavity2dRe100Serial'))
  tests.append(TestCase(description='2d cylinder flow at Re=40',
                        directory='{}/examples/2d/cylinder/Re40'.format(build_directory),
                        command='make cylinder2dRe40'))
  tests.append(TestCase(description='2d cylinder flow at Re=40 with y-periodic boundary conditions',
                        directory='{}/examples/2d/cylinder/Re40PeriodicDomain'.format(build_directory),
                        command='make cylinder2dRe40Periodic',
                        periodic=['y']))
  tests.append(TestCase(description='3d cavity flow at Re=100 with x-periodic boundary conditions',
                        directory='{}/examples/3d/lidDrivenCavity/Re100PeriodicX'.format(build_directory),
                        command='make lidDrivenCavity3dRe100PeriodicX',
                        periodic=['x']))
  return tests


class TestCase(object):
  """Contains information about a test-case."""
  def __init__(self, description, directory, command, periodic=[]):
    """Initializes the test-case.

    Parameters
    ----------
    description: str
      Description of the test-case.
    directory: str
      Directory of the test-case.
    command: str
      Command-line to execute to run the test-case.
    periodic: list(str)
      Directions with periodic boundary conditions.
    """
    self.description = description
    self.directory = directory
    self.reference = None
    self.command = command
    self.periodic = periodic
    self.differences = []
    self.passed = True
    self.saved = False

  def print_info(self):
    """Prints some info about the test-case."""
    print('\n----------')
    print(  'Test-case: {}'.format(self.description))
    print(  '----------')
    print('Directory: {}'.format(self.directory))
    print('Command-line: {}'.format(self.command))

  def run(self):
    """Runs the test-case."""
    os.system(self.command)

  def compare(self,save=False):
    """Compares the numerical solution with a reference.

    Parameters
    ----------
    save: bool
      Save the new numerical solution; default: False.
    """
    self.reference = self.directory.replace('examples', 'regressionAnalysis')
    if not os.path.isdir(self.reference):
      print('\nWARNING: no reference available. Skipping comparison.')
      self.passed = False
      self.differences.append('reference solution not available (new one will be saved)')
      self.save()
    else:
      print('\nReference: {}'.format(self.reference))
      print('Comparing numerical solution to reference...')
      self.compare_grid()
      self.compare_velocity()
      self.compare_pressure()
      self.compare_forces()
      if self.passed and save:
        self.save()

  def compare_arrays(self, array1, array2, tag):
    """Performs element-wise comparison of two given arrays.

    Parameters
    ----------
    array1, array2: numpy.ndarray
      The two arrays to be compared.
    tag: str
      A description of the arrays to be compared.
    """
    if not numpy.allclose(array1, array2, atol=1.0E-06):
      self.passed = False
      self.differences.append('difference in {}'.format(tag))

  def compare_grid(self):
    """Compares the mesh-grid with a reference one."""
    grid = ioPetIBM.read_grid(self.directory)
    grid_reference = ioPetIBM.read_grid(self.reference)
    for i, direction in enumerate(grid):
      self.compare_arrays(direction, grid_reference[i], 
                          tag='grid[{}]'.format(i))

  def compare_velocity(self):
    """Compares the velocity field node by node with a reference."""
    time_step = ioPetIBM.get_time_steps(self.directory)[-1]
    grid = ioPetIBM.read_grid(self.directory)
    velocity = ioPetIBM.read_velocity(self.directory, time_step, grid, 
                                      periodic=self.periodic)
    grid = ioPetIBM.read_grid(self.reference)
    velocity_reference = ioPetIBM.read_velocity(self.reference, time_step, grid, 
                                                periodic=self.periodic)
    for i, component in enumerate(velocity):
      self.compare_arrays(component.x, velocity_reference[i].x, 
                          tag='velocity[{}]: x-nodes'.format(i))
      self.compare_arrays(component.y, velocity_reference[i].y, 
                          tag='velocity[{}]: y-nodes'.format(i))
      try:
        self.compare_arrays(component.z, velocity_reference[i].z, 
                            tag='velocity[{}]: z-nodes'.format(i))
      except:
        pass
      self.compare_arrays(component.values, velocity_reference[i].values, 
                          tag='velocity[{}]: values'.format(i))

  def compare_pressure(self):
    """Compares the pressure field node by node with a reference."""
    time_step = ioPetIBM.get_time_steps(self.directory)[-1]
    grid = ioPetIBM.read_grid(self.directory)
    pressure = ioPetIBM.read_pressure(self.directory, time_step, grid)
    grid = ioPetIBM.read_grid(self.reference)
    pressure_reference = ioPetIBM.read_pressure(self.reference, time_step, grid)
    self.compare_arrays(pressure.x, pressure_reference.x, 
                        tag='pressure: x-nodes')
    self.compare_arrays(pressure.y, pressure_reference.y, 
                        tag='pressure: y-nodes')
    try:
      self.compare_arrays(pressure.z, pressure_reference.z, 
                          tag='pressure: z-nodes')
    except:
      pass
    self.compare_arrays(pressure.values, pressure_reference.values, 
                        tag='pressure: values')

  def compare_forces(self):
    """Compares forces acting on immersed boundaries (if applicable)."""
    try:
      with open('{}/forces.txt'.format(self.directory), 'r') as infile:
        forces = numpy.loadtxt(infile, dtype=float)
      with open('{}/forces.txt'.format(self.reference), 'r') as infile:
        forces_reference = numpy.loadtxt(infile, dtype=float)
      self.compare_arrays(forces, forces_reference, tag='forces')
    except:
      pass
    
  def save(self):
    """Saves the numerical solution into a folder."""
    print('Copy numerical solution of {} into {}\n'.format(self.directory,
                                                           self.reference))
    if os.path.isdir(self.reference):
      shutil.rmtree(self.reference)
    shutil.copytree(self.directory, self.reference)
    self.saved = True

  def write(self, file_path):
    """Writes the results of the regression analysis into a file.

    Parameters
    ----------
    file_path: str
      Path of the file where to write.
    """
    with open(file_path, 'a') as outfile:
      outfile.write('\n----------\n')
      outfile.write(  'Test-case: {}\n'.format(self.description))
      outfile.write(  '----------\n')
      outfile.write('directory: {}\n'.format(self.directory))
      outfile.write('reference: {}\n'.format(self.reference))
      outfile.write('passed: {}\n'.format('yes' if self.passed else 'no'))
      if not self.passed:
        outfile.write('reason(s): {}\n'.format('\n\t'.join(self.differences)))
      outfile.write('saved: {}\n'.format('yes' if self.saved else 'no'))


def print_configuration(build_directory):
  """Prints the configuration used to build PetIBM.

  Parameters
  ----------
  build_directory: str
    Directory of the PetIBM build.
  """
  print('\n=============')
  print(  'Configuration')
  print(  '=============\n')
  print('Build directory: {}'.format(build_directory))
  print('PETSC_DIR: {}'.format(os.environ['PETSC_DIR']))
  print('PETSC_ARCH: {}'.format(os.environ['PETSC_ARCH']))


def compile_PetIBM(build_directory):
  """Compiles PetIBM.

  Parameters
  ----------
  build_directory: str
    Directory of the PetIBM build.
  """
  print('\n==============')
  print(  'Compile PetIBM')
  print(  '==============\n')
  os.chdir(build_directory)
  os.system('make clean')
  os.system('make all')
  os.system('make check')
  # check existence of executables
  for executable in ['petibm2d', 'petibm3d']:
    if not os.path.exists('{}/src/{}'.format(build_directory, executable)):
      print('ERROR: {} does not exist.'.format(executable))
      sys.exit()


def perform_regression_analysis(build_directory, tests, run, save):
  """Runs test-cases and performs regression analysis.

  Parameters
  ----------
  build_directory: str
    Directory of the PetIBM build.
  tests: dict(TestCase)
    Dictionary containing all the test-cases to run.
  run: bool
    Do you want to run the test-cases?
  save: bool
    Do you want to save the numerical results for future regression analysis?
  """
  print('\n===================')
  print(  'Regression analysis')
  print(  '===================\n')
  os.chdir('{}/examples'.format(build_directory))
  os.system('make examples')
  # create regressionAnalysis folder if need
  regression_directory = '{}/regressionAnalysis'.format(build_directory)
  if not os.path.isdir(regression_directory):
    os.makedirs(regression_directory)
  # write intro to summary file
  summary_path = '{}/summary.txt'.format(regression_directory)
  with open(summary_path, 'w') as outfile:
    outfile.write('Regression analysis performed on {}\n'.format(time.strftime('%m/%d/%Y')))
  # run test-cases
  print('Looping over the test-cases...')
  global_passed = True
  for test in tests:
    test.print_info()
    if run:
      test.run()
    test.compare(save=save)
    test.write(summary_path)
    if not test.passed:
      global_passed = False

  print('\nPassed: {}'.format('yes' if global_passed else 'no'))
  if not global_passed:
    print('Check {} for more info.'.format(summary_path))


def main():
  """Cleans, compiles PetIBM, then runs test-cases to ensure that the output
  matches the solutions of the previous version."""
  # parse command-line
  args = read_inputs()
  print_configuration(args.build_directory)
  if args.compile:
    compile_PetIBM(args.build_directory)
  tests = define_test_cases(args.build_directory)
  perform_regression_analysis(args.build_directory, tests, args.run, args.save)


if __name__ == '__main__':
  print('\n[{}] START\n'.format(os.path.basename(__file__)))
  main()
  print('\n[{}] END\n'.format(os.path.basename(__file__)))
