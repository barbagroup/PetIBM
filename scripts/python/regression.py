#!/usr/bin/env/ python

# file: regression.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Performs regression tests for PetIBM.


import os
import sys
import time
import shutil

sys.path.append('{}/scripts/python'.format(os.environ['PETIBM_DIR']))
import ioPetIBM


class TestCase(object):
  """Contains information about a test-case."""
  def __init__(self, directory, 
               petibmexec=None, mpiexec=None, n=None, arguments=None,
               periodic=[]):
    """Initializes the test-case with arguments for the command-line."""
    self.directory = directory
    self.basename = os.path.basename(directory)
    self.petibmexec = petibmexec
    self.mpiexec = mpiexec
    self.n = n
    self.arguments = arguments
    self.periodic = periodic
    self.differences = []

  def run(self):
    """Run the test-case."""
    print('\nRun test-case: {}\n'.format(self.directory))
    os.system('{} -n {} {} -caseFolder {} {}'.format(self.mpiexec, self.n, 
                                                     self.petibmexec,
                                                     self.directory,
                                                     self.arguments))

  def regression(self, regression_directory):
    """Compares the numerical solution with the solution from previous version."""
    previous_solution = '{}/{}'.format(regression_directory, self.basename)
    if not os.path.isdir(previous_solution):
      print('WARNING: no numerical solution to compare with')
      self.save(regression_directory)
      self.passed = False
      self.differences.append('\tno previous numerical solution to compare with\n')
      return
    self.compare_grid(previous_solution)
    self.compare_velocity(previous_solution)
    self.compare_pressure(previous_solution)
    self.passed = (True if not self.differences else False)

  def compare_arrays(self, array1, array2, tag):
    same = (array1 == array2).all()
    if not same:
      self.differences.append(tag)

  def compare_grid(self, comparator):
    grid = ioPetIBM.read_grid(self.directory)
    previous_grid = ioPetIBM.read_grid(comparator)
    self.compare_arrays(grid, previous_grid, 'grid')

  def compare_velocity(self, comparator):
    time_step = ioPetIBM.get_time_steps(self.directory)[-1]
    grid = ioPetIBM.read_grid(self.directory)
    velocity = ioPetIBM.read_velocity(self.directory, time_step, grid, 
                                      periodic=self.periodic)
    grid = ioPetIBM.read_grid(comparator)
    previous_velocity = ioPetIBM.read_velocity(comparator, time_step, grid, 
                                               periodic=self.periodic)
    for i, component in enumerate(velocity):
      self.compare_arrays(component.x, previous_velocity[i].x, 
                          'velocity[{}]: x-nodes'.format(i))
      self.compare_arrays(component.y, previous_velocity[i].y, 
                          'velocity[{}]: y-nodes'.format(i))
      if component.z:
        self.compare_arrays(component.z, previous_velocity[i].z, 
                            'velocity[{}]: z-nodes'.format(i))
      self.compare_arrays(component.values, previous_velocity[i].values, 
                          'velocity[{}]: values'.format(i))

  def compare_pressure(self, comparator):
    time_step = ioPetIBM.get_time_steps(self.directory)[-1]
    grid = ioPetIBM.read_grid(self.directory)
    pressure = ioPetIBM.read_pressure(self.directory, time_step, grid)
    grid = ioPetIBM.read_grid(comparator)
    previous_pressure = ioPetIBM.read_pressure(comparator, time_step, grid)
    self.compare_arrays(pressure.x, previous_pressure.x, 'pressure: x-nodes')
    self.compare_arrays(pressure.y, previous_pressure.y, 'pressure: y-nodes')
    if pressure.z:
      self.compare_arrays(pressure.z, previous_pressure.z, 'pressure: z-nodes')
    self.compare_arrays(pressure.values, previous_pressure.values, 'pressure: values')

  def save(self, regression_directory):
    """Saves the numerical solution into the regression directory."""
    print('\nCopy numerical solution of {} into {}\n'.format(self.basename,
                                                             regression_directory))
    destination = '{}/{}'.format(regression_directory, self.basename)
    if os.path.isdir(destination):
      shutil.rmtree(destination)
    shutil.copytree(self.directory, destination)

  def write(self, file_path):
    with open(file_path, 'a') as outfile:
      if self.passed:
        outfile.write('\n{}: OK\n'.format(self.basename))
      else:
        outfile.write('\n{}: NOT OK\n'.format(self.basename))
        outfile.write('\t\n'.join(self.differences))    


def main():
  """Cleans, compiles PetIBM, then runs test-cases to ensure that the output
  matches the solutions of the previous version."""

  os.environ['MPIEXEC'] = '{}/arch-linux2-c-opt/bin/mpiexec'.format(os.environ['PETSC_DIR'])

  # clean and compile PetIBM
  print('\nClean and compile PetIBM\n')
  os.chdir(os.environ['PETIBM_DIR'])
  #os.system('make cleanall')
  #os.system('make')
  # check existence of executables
  os.environ['PETIBM2D'] = '{}/bin/PetIBM2d'.format(os.environ['PETIBM_DIR'])
  os.environ['PETIBM3D'] = '{}/bin/PetIBM3d'.format(os.environ['PETIBM_DIR'])
  for executable in [os.environ['PETIBM2D'], os.environ['PETIBM3D']]:
    if not os.path.exists(executable):
      print('ERROR: could not build {}'.format(executable))
      sys.exit()

  tests = {}

  # test-case: 2d lid-driven cavity flow (Re=100)
  case = '2d lid-driven cavity (Re=100)'
  tests[case] = TestCase('/home/mesnardo/tests_PetIBM/lidDrivenCavity2dRe100')
  tests[case].petibmexec = os.environ['PETIBM2D']
  tests[case].mpiexec = os.environ['MPIEXEC']
  tests[case].n = 1
  tests[case].arguments = '-sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1'

  # test-case: 2d cylinder (Re=40)
  case = '2d cylinder (Re=40)'
  tests[case] = TestCase('/home/mesnardo/tests_PetIBM/cylinder2dRe40')
  tests[case].petibmexec = os.environ['PETIBM2D']
  tests[case].mpiexec = os.environ['MPIEXEC']
  tests[case].n = 2
  tests[case].arguments = '-sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1'

  # test-case: 3d sphere (Re=300)
  case = '3d sphere (Re=300)'
  # tests[case] = TestCase('/home/mesnardo/tests_PetIBM/sphere3dRe300')
  # tests[case].petibmexec = os.environ['PETIBM3D']
  # tests[case].mpiexec = os.environ['MPIEXEC']
  # tests[case].n = 4
  # tests[case].arguments = '-sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1'

  # create regression directory
  regression_directory = '{}/regression_results'.format(os.environ['PETIBM_DIR'])
  if not os.path.isdir(regression_directory):
    os.makedirs(regression_directory)
  regression_file_path = '{}/results.txt'.format(regression_directory)
  with open(regression_file_path, 'w') as outfile:
    outfile.write('Regression performed on {}\n\n'.format(time.strftime('%m/%d/%Y')))

  # run test-cases
  for test in tests.itervalues():
    test.run()
    test.regression(regression_directory)
    test.write(regression_file_path)
    if test.passed:
      test.save(regression_directory)


  # do comparison

if __name__ == '__main__':
  print('\n[{}] START\n'.format(os.path.basename(__file__)))
  main()
  print('\n[{}] END\n'.format(os.path.basename(__file__)))