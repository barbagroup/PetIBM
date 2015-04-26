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
  parser.add_argument('--save', dest='save', action='store_true', 
                      help='saves the new numerical solutions')
  parser.add_argument('--no-compile', dest='compile', action='store_false',
                      help='skips PetIBM compilation')
  parser.add_argument('--no-run', dest='run', action='store_false',
                      help='does not run the test-cases')
  parser.set_defaults(compile=True, run=True)
  # parse command-line
  return parser.parse_args()


def define_test_cases():
  """Defines the test-cases to be run."""
  # dictionary that will contains test-cases info
  tests = {}

  # test-case: 2d lid-driven cavity flow (Re=100)
  case = '2d lid-driven cavity (Re=100)'
  tests[case] = TestCase('{}/cases/2d/lidDrivenCavity/Re100'.format(os.environ['PETIBM_DIR']))
  tests[case].petibmexec = os.environ['PETIBM2D']
  tests[case].mpiexec = os.environ['MPIEXEC']
  tests[case].n = 1
  tests[case].arguments = ('-sys2_pc_type gamg -sys2_pc_gamg_type agg '
                           '-sys2_pc_gamg_agg_nsmooths 1')
  # test-case: 2d cylinder (Re=40)
  case = '2d cylinder (Re=40)'
  tests[case] = TestCase('{}/cases/2d/cylinder/Re40'.format(os.environ['PETIBM_DIR']))
  tests[case].petibmexec = os.environ['PETIBM2D']
  tests[case].mpiexec = os.environ['MPIEXEC']
  tests[case].n = 2
  tests[case].arguments = ('-sys2_pc_type gamg -sys2_pc_gamg_type agg '
                           '-sys2_pc_gamg_agg_nsmooths 1')
  # test-case: 3d lid-driven cavity flow (Re=100, periodic=['x'])
  case = '3d lid-driven cavity (Re=100)'
  tests[case] = TestCase('{}/cases/3d/lidDrivenCavity/Re100PeriodicX'.format(os.environ['PETIBM_DIR']))
  tests[case].periodic = ['x']
  tests[case].petibmexec = os.environ['PETIBM3D']
  tests[case].mpiexec = os.environ['MPIEXEC']
  tests[case].n = 2
  tests[case].arguments = ('-sys2_pc_type gamg -sys2_pc_gamg_type agg '
                           '-sys2_pc_gamg_agg_nsmooths 1')
  # test-case: 3d sphere (Re=300)
  case = '3d sphere (Re=300)'
  tests[case] = TestCase('{}/cases/3d/sphere/Re300'.format(os.environ['PETIBM_DIR']))
  tests[case].petibmexec = os.environ['PETIBM3D']
  tests[case].mpiexec = os.environ['MPIEXEC']
  tests[case].n = 4
  tests[case].arguments = ('-sys2_pc_type gamg -sys2_pc_gamg_type agg '
                           '-sys2_pc_gamg_agg_nsmooths 1')
  return tests


class TestCase(object):
  """Contains information about a test-case."""
  def __init__(self, directory, 
               petibmexec=None, mpiexec=None, n=1, arguments=None, periodic=[]):
    """Initializes the test-case with arguments for the command-line.

    Parameters
    ----------
    directory: str
      Directory of the simulation.
    petibmexec: None or str
      Path of the PetIBM binary executable.
    mpiexec: None or str
      Path of the MPI executable.
    n: int
      Number of processes to use; default: 1.
    arguments: None or list(str)
      PETSc command-line arguments.
    periodic: list(str)
      Directions with periodic boundary conditions.
    """
    self.directory = directory
    self.petibmexec = petibmexec
    self.mpiexec = mpiexec
    self.n = n
    self.arguments = arguments
    self.periodic = periodic
    self.differences = []
    self.passed = True

  def run(self):
    """Runs the test-case."""
    print('\nRun test-case: {}\n'.format(self.directory))
    os.system('{} -n {} {} -caseFolder {} {}'.format(self.mpiexec, self.n, 
                                                     self.petibmexec,
                                                     self.directory,
                                                     self.arguments))

  def regression_analysis(self, directory):
    """Compares the numerical solution with the solution from previous version.

    Parameters
    ----------
    directory: str
      Directory where is stored the reference case folder.
    """
    self.basename = os.path.relpath(self.directory, os.environ['PETIBM_DIR'])
    self.reference = '{}/{}'.format(directory, self.basename)
    if not os.path.isdir(self.reference):
      print('\nWARNING: no numerical solution to compare with\n')
      self.save(directory)
      self.passed = False
      self.differences.append('no numerical solution to compare with')
      return
    self.compare_grid()
    self.compare_velocity()
    self.compare_pressure()
    self.compare_forces()

  def compare_arrays(self, array1, array2, tag):
    """Performs element-wise comparison of two given arrays.

    Parameters
    ----------
    array1, array2: numpy.ndarray
      The two arrays to be compared.
    tag: str
      A description of the arrays to be compared.
    """
    if not numpy.allclose(array1, array2, rtol=1.0E-06):
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
    
  def save(self, directory):
    """Saves the numerical solution into a folder.

    Parameters
    ----------
    directory: str
      Directory where the numerical solution will be saved.
    """
    print('\nCopy numerical solution of {} into {}\n'.format(self.basename,
                                                             directory))
    destination = '{}/{}'.format(directory, self.basename)
    if os.path.isdir(destination):
      shutil.rmtree(destination)
    shutil.copytree(self.directory, destination)

  def write(self, file_path):
    """Writes the results of the regression analysis into a file.

    Parameters
    ----------
    file_path: str
      Path of the file.
    """
    with open(file_path, 'a') as outfile:
      message = ('OK' if self.passed else 'NOT OK')
      outfile.write('\n\n{}: {}\n'.format(self.basename, message))
      outfile.write('\t\n'.join(self.differences))


def main():
  """Cleans, compiles PetIBM, then runs test-cases to ensure that the output
  matches the solutions of the previous version."""

  # parse command-line
  args = read_inputs()

  os.environ['MPIEXEC'] = '{}/arch-linux2-c-opt/bin/mpiexec'.format(os.environ['PETSC_DIR'])

  # clean and compile PetIBM
  if args.compile:
    print('\nClean and compile PetIBM\n')
    os.chdir(os.environ['PETIBM_DIR'])
    os.system('make cleanall')
    os.system('make')
  # check existence of executables
  os.environ['PETIBM2D'] = '{}/bin/PetIBM2d'.format(os.environ['PETIBM_DIR'])
  os.environ['PETIBM3D'] = '{}/bin/PetIBM3d'.format(os.environ['PETIBM_DIR'])
  for executable in [os.environ['PETIBM2D'], os.environ['PETIBM3D']]:
    if not os.path.exists(executable):
      print('ERROR: could not build {}'.format(executable))
      sys.exit()

  tests = define_test_cases()

  # create regression directory
  regression_directory = '{}/regression_analysis'.format(os.environ['PETIBM_DIR'])
  if not os.path.isdir(regression_directory):
    os.makedirs(regression_directory)
  regression_file_path = '{}/summary.txt'.format(regression_directory)
  with open(regression_file_path, 'w') as outfile:
    outfile.write('Regression analysis performed on {}\n\n'.format(time.strftime('%m/%d/%Y')))

  # run test-cases
  for test in tests.itervalues():
    if args.run:
      test.run()
    test.regression_analysis(regression_directory)
    test.write(regression_file_path)
    if test.passed and args.save:
      test.save(regression_directory)


if __name__ == '__main__':
  print('\n[{}] START\n'.format(os.path.basename(__file__)))
  main()
  print('\n[{}] END\n'.format(os.path.basename(__file__)))