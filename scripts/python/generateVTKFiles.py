"""
Writes .vtk files out of the saved PETSc solutions.
"""

import os
import argparse

import ioPetIBM


def parse_command_line():
  """
  Parses the command-line.

  Returns
  -------
  args: namespace
    Database with arguments parsed from the command-line.
  """
  formatter_class = argparse.ArgumentDefaultsHelpFormatter
  parser = argparse.ArgumentParser(description='Converts PETSc output to VTK '
                                               'format for 3D case',
                                   formatter_class=formatter_class)
  parser.add_argument('--directory',
                      dest='directory',
                      type=str,
                      default=os.getcwd(),
                      help='directory of the simulation')
  parser.add_argument('--variables', '-v',
                      dest='variables',
                      type=str,
                      nargs='+',
                      default=['velocity', 'pressure'],
                      help='list of variables to generate '
                           '(velocity, pressure)')
  parser.add_argument('--bottom-left', '-bl',
                      dest='bottom_left',
                      type=float,
                      nargs='+',
                      default=[float('-inf'), float('-inf'), float('-inf')],
                      help='coordinates of the bottom-left corner')
  parser.add_argument('--top-right', '-tr',
                      dest='top_right',
                      type=float,
                      nargs='+',
                      default=[float('inf'), float('inf'), float('inf')],
                      help='coordinates of the top-right corner')
  parser.add_argument('--time-steps', '-t',
                      dest='time_steps',
                      type=int,
                      nargs='+',
                      default=[],
                      help='time-steps to convert (start, end, increment)')
  parser.add_argument('--stride', '-s',
                      dest='stride',
                      type=int,
                      default=1,
                      help='stride at which vector are written')
  parser.add_argument('--periodic', '-p',
                      dest='periodic',
                      type=str,
                      nargs='+',
                      default=[],
                      help='direction(s) (x and/or y and/or z) '
                            'with periodic boundary conditions')
  return parser.parse_args()


def interpolate_cell_centers(velocity):
  """
  Interpolates the velocity field at the cell-centers.

  Parameters
  ----------
  velocity: list of Field objects
    Velocity field on a staggered grid.

  Returns
  -------
  velocity: list of Field objects
    Velocity field at cell-centers.
  """
  dim3 = (True if len(velocity) == 3 else False)
  x_centers, y_centers = velocity[1].x[1:-1], velocity[0].y[1:-1]
  u, v = velocity[0].values, velocity[1].values
  if dim3:
    z_centers = velocity[0].z[1:-1]
    w = velocity[2].values
    u = 0.5 * (u[1:-1, 1:-1, :-1] + u[1:-1, 1:-1, 1:])
    v = 0.5 * (v[1:-1, :-1, 1:-1] + v[1:-1:, 1:, 1:-1])
    w = 0.5 * (w[:-1, 1:-1, 1:-1] + w[1:, 1:-1, 1:-1])
    # tests
    assert (z_centers.size, y_centers.size, x_centers.size) == u.shape
    assert (z_centers.size, y_centers.size, x_centers.size) == v.shape
    assert (z_centers.size, y_centers.size, x_centers.size) == w.shape
    return [ioPetIBM.Field(x=x_centers, y=y_centers, z=z_centers, values=u),
            ioPetIBM.Field(x=x_centers, y=y_centers, z=z_centers, values=v),
            ioPetIBM.Field(x=x_centers, y=y_centers, z=z_centers, values=w)]
  else:
    u = 0.5 * (u[1:-1, :-1] + u[1:-1, 1:])
    v = 0.5 * (v[:-1, 1:-1] + v[1:, 1:-1])
    # tests
    assert (y_centers.size, x_centers.size) == u.shape
    assert (y_centers.size, x_centers.size) == v.shape
    return [ioPetIBM.Field(x=x_centers, y=y_centers, values=u),
            ioPetIBM.Field(x=x_centers, y=y_centers, values=v)]


def main(args):
  """
  Generates .vtk files out of the saved PETSc solutions.

  Parameters
  ----------
  args: namespace
    Databased with arguments parsed from the command-line.
  """
  print('[directory] {}'.format(args.directory))
  print('[variables] {}'.format(args.variables))

  # list of time-steps to post-process
  time_steps = ioPetIBM.get_time_steps(args.directory, args.time_steps)

  # read mesh grid
  grid = ioPetIBM.read_grid(args.directory)

  for time_step in time_steps:
    if 'velocity' in args.variables:
      velocity = ioPetIBM.read_velocity(args.directory, time_step, grid,
                                        periodic=args.periodic)
      # need to get velocity at cell-centers, not staggered arrangement
      velocity = interpolate_cell_centers(velocity)
      ioPetIBM.write_vtk(velocity, args.directory, time_step,
                         name='velocity',
                         view=[args.bottom_left, args.top_right],
                         stride=args.stride)
    if 'pressure' in args.variables:
      pressure = ioPetIBM.read_pressure(args.directory, time_step, grid)
      ioPetIBM.write_vtk(pressure, args.directory, time_step,
                         name='pressure',
                         view=[args.bottom_left, args.top_right],
                         stride=args.stride)


if __name__ == '__main__':
  print('\n[{}] START\n'.format(os.path.basename(__file__)))
  args = parse_command_line()
  main(args)
  print('\n[{}] END\n'.format(os.path.basename(__file__)))
