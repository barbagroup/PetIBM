"""
Plots the centerline velocities for lid-driven cavity flow case
and compare with experimental data from Ghia et al. (1982).
"""

import os
import sys
import argparse

import numpy
from matplotlib import pyplot

sys.path.append(os.path.join(os.environ['PETIBM_DIR'], 'scripts', 'python'))
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
  parser = argparse.ArgumentParser(description='Plots centerline velocities '
                                               'for lid-driven cavity flow',
                                   formatter_class=formatter_class)
  parser.add_argument('--directory',
                      dest='directory',
                      type=str,
                      default=os.getcwd(),
                      help='directory of the simulation')
  parser.add_argument('--Re', '-Re',
                      dest='Re',
                      type=str,
                      default='100',
                      help='Reynolds number of the flow')
  parser.add_argument('--time-step', '-t',
                      dest='time_step',
                      type=int,
                      default=None,
                      help='time-step to plot')
  return parser.parse_args()


def main(args):
  """
  Plots and writes the velocity components at the centerline of the cavity
  and compares with experimental results form Ghia et al. (1982).

  Parameters
  ----------
  args: namespace
    Database with arguments parsed from the command-line.
  """
  # column indices in file with experimental results
  cols = {'100': {'u': 1, 'v': 7},
          '1000': {'u': 2, 'v': 8},
          '3200': {'u': 3, 'v': 9},
          '5000': {'u': 4, 'v': 10},
          '10000': {'u': 5, 'v': 11}}
  experimental_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                   'data',
                                   'cavity-GGS82.txt')

  print('[case directory] {}'.format(args.directory))
  Re = args.Re

  if not args.time_step:
    args.time_step = ioPetIBM.get_time_steps(directory=args.directory)[-1]

  # read mesh grid
  grid = ioPetIBM.read_grid(directory=args.directory)
  nx, ny = grid[0].size - 1, grid[1].size - 1

  # create directory where images will be saved
  images_directory = os.path.join(args.directory, 'images')
  print('[images directory] {}'.format(images_directory))
  if not os.path.isdir(images_directory):
    os.makedirs(images_directory)

  # create directory where data will be saved
  data_directory = os.path.join(args.directory, 'data')
  print('[data directory] {}'.format(data_directory))
  if not os.path.isdir(data_directory):
    os.makedirs(data_directory)

  # get velocity field
  u, v = ioPetIBM.read_velocity(args.time_step, grid,
                                directory=args.directory)

  # load default style for figures
  pyplot.style.use(os.path.join(os.environ['PETIBM_DIR'],
                                'scripts',
                                'python',
                                'style',
                                'style_PetIBM.mplstyle'))

  # plot and write u-velocity along y-centerline
  print('\nPlot u-velocity along vertical centerline ...')
  pyplot.grid(True)
  pyplot.xlabel('y',)
  pyplot.ylabel('u-velocity along vertical centerline')
  if nx % 2 == 0:
    u_centerline = u.values[:, nx / 2 - 1]
  else:
    u_centerline = 0.5 * (u.values[:, nx / 2 - 1] + u.values[:, nx / 2])
  pyplot.plot(u.y, u_centerline, label='PetIBM',
              color='k', linestyle='-', linewidth=1)
  if Re in list(cols.keys()):
    print('\nComparison with Ghia et al. (1982)')
    with open(experimental_file, 'r') as infile:
      y_exp, u_exp = numpy.loadtxt(infile, dtype=float,
                                   usecols=(0, cols[Re]['u']), unpack=True)
    pyplot.plot(y_exp, u_exp, label='Ghia et al. (1982)',
                linewidth=0, marker='o')
  pyplot.legend()
  pyplot.title('Lid-driven cavity flow at Re={} (mesh: {}x{})'
               ''.format(Re, nx, ny))
  pyplot.savefig(os.path.join(images_directory,
                              'uVelocityCenterlineRe{}_{}x{}.png'
                              .format(Re, nx, ny)))
  pyplot.close()
  print('\nWrite u-velocity along vertical centerline ...')
  data_path = os.path.join(data_directory,
                           'uVelocityCenterlineRe{}_{}x{}.txt'
                           .format(Re, nx, ny))
  with open(data_path, 'w') as outfile:
    numpy.savetxt(outfile, numpy.c_[u.y, u_centerline],
                  fmt='%.6f', delimiter='\t',
                  header='u-velocity along vertical centerline\n'
                         'Re={}, mesh: {}x{}\n[y]\t[u]'.format(Re, nx, ny))

  # plot and write v-velocity along x-centerline
  print('\nPlot v-velocity along horizontal centerline ...')
  pyplot.grid(True)
  pyplot.xlabel('x')
  pyplot.ylabel('v-velocity along horizontal centerline')
  if ny % 2 == 0:
    v_centerline = v.values[ny / 2 - 1, :]
  else:
    v_centerline = 0.5 * (v.values[ny / 2 - 1, :] + v.values[ny / 2, :])
  pyplot.plot(v.x, v_centerline, label='PetIBM',
              color='k', linestyle='-', linewidth=1)
  if Re in list(cols.keys()):
    print('Comparison with Ghia et al. (1982)')
    with open(experimental_file, 'r') as infile:
      x_exp, u_exp = numpy.loadtxt(infile, dtype=float,
                                   usecols=(6, cols[Re]['v']), unpack=True)
    pyplot.plot(x_exp, u_exp, label='Ghia et al. (1982)',
                linewidth=0, marker='o')
  pyplot.legend()
  pyplot.title('Lid-driven cavity flow at Re={} (mesh: {}x{})'
               ''.format(Re, nx, ny))
  pyplot.savefig(os.path.join(images_directory,
                              'vVelocityCenterlineRe{}_{}x{}.png'
                              .format(Re, nx, ny)))
  pyplot.close()
  print('\nwrite v-velocity along horizontal centerline ...')
  data_path = os.path.join(data_directory,
                           'vVelocityCenterlineRe{}_{}x{}.txt'
                           .format(Re, nx, ny))
  with open(data_path, 'w') as outfile:
    numpy.savetxt(outfile, numpy.c_[v.x, v_centerline],
                  fmt='%.6f', delimiter='\t',
                  header='v-velocity along horizontal centerline\n'
                         'Re={}, mesh: {}x{}\n[x]\t[v]'.format(Re, nx, ny))


if __name__ == '__main__':
  print('\n[{}] START\n'.format(os.path.basename(__file__)))
  args = parse_command_line()
  main(args)
  print('\n[{}] END\n'.format(os.path.basename(__file__)))
