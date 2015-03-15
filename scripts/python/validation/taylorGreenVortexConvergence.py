#!/usr/bin/env python

# file: taylorGreenVortexConvergence.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Plots the grid-convergence for the Taylor-Green vortex case.


import os
import sys
import argparse
import math

import numpy
from matplotlib import pyplot

sys.path.append('{}/scripts/python'.format(os.environ['PETIBM_DIR']))
import ioPetIBM


def read_inputs():
  """Parses the command-line."""
  # create parser
  parser = argparse.ArgumentParser(description='Convergence for the '
                                               'Taylor-Green vortex case',
                        formatter_class= argparse.ArgumentDefaultsHelpFormatter)
  # fill parser with arguments
  parser.add_argument('--directory', dest='directory', type=str,
                      default=os.getcwd(),
                      help='directory containing the simulation folders')
  parser.add_argument('--Re', '-Re', dest='Re', type=float, default=100.0,
                      help='Reynolds number of the simulation')
  parser.add_argument('--time', '-t', dest='time', type=float, default=0.5,
                      help='time at which the error will be computed')
  parser.add_argument('--time-step', '-ts', dest='time_step', type=int, 
                      default=1000,
                      help='time-step at which the solution will be read')
  parser.add_argument('--amplitude', '-a', dest='amplitude', type=float, 
                      default=1.0, help='amplitude of the Taylor-Green vortex')
  parser.add_argument('--no-save', dest='save', action='store_false',
                      help='does not save the figure')
  parser.add_argument('--output', '-o', dest='output', type=str, 
                      default='grid_convergence',
                      help='name of the .png file saved')
  parser.add_argument('--no-show', dest='show', action='store_false',
                      help='does not display the figure')
  parser.set_defaults(save=True, show=True)
  # parse command-line
  return parser.parse_args()


def taylor_green_vortex(x, y, V=1.0, time=0.0, Re=100.0, name=None):
  """Computes the analytical solution of the 2D Taylor-Green vortex.

  Arguments
  ---------
  x, y -- coordinates in the x- and y- directions
  V -- amplitude of the sinusoidal velocity field (default 1.0)
  time -- time at which the solution is computed (default 0.0)
  Re -- Reynolds number of the flow (default 100.0)
  name -- name of the variable to compute (default None)
  """
  x = math.pi - 2.0*math.pi*(x-x[0])/(x[-1]-x[0])
  y = math.pi - 2.0*math.pi*(y-y[0])/(y[-1]-y[0])
  X, Y = numpy.meshgrid(x, y)
  if name == 'u-velocity':
    return V*numpy.sin(X)*numpy.cos(Y)*math.exp(-2.0*time/Re)
  elif name == 'v-velocity':
    return -V*numpy.cos(X)*numpy.sin(Y)*math.exp(-2.0*time/Re)
  elif name == 'vorticity':
    return 2.0*numpy.sin(X)*numpy.sin(Y)*math.exp(-2.0*time/Re)
  elif name == 'pressure':
    return 0.25*(numpy.cos(2.0*X)+numpy.sin(2.0*Y))*math.exp(-4.0*time/Re)
  else:
    return 'Error: variable name unknown'


def main():
  """Plots the grid convergence for the Taylor-Green vortex case."""
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
    # read mesh grid
    x, y = ioPetIBM.read_grid(case['directory'])
    cases[i]['grid-spacing'] = (x[-1]-x[0])/(x.size-1)
    # velocity components
    u, v = ioPetIBM.read_velocity(case['directory'], parameters.time_step, 
                                  [x, y], periodic=['x', 'y'])
    u_analytical = taylor_green_vortex(u['x'], u['y'], 
                                       V=parameters.amplitude, 
                                       time=parameters.time , 
                                       Re=parameters.Re, 
                                       name='u-velocity')
    v_analytical = taylor_green_vortex(v['x'], v['y'], 
                                       V=parameters.amplitude, 
                                       time=parameters.time , 
                                       Re=parameters.Re, 
                                       name='v-velocity')
    cases[i]['u'] = u['values'][ratio-1::ratio, ratio-1::ratio]
    cases[i]['u-error'] = numpy.linalg.norm(u['values']-u_analytical)/u['values'].size
    cases[i]['v'] = v['values'][ratio-1::ratio, ratio-1::ratio]
    cases[i]['v-error'] = numpy.linalg.norm(v['values']-v_analytical)/v['values'].size
    # pressure
    p = ioPetIBM.read_pressure(case['directory'], parameters.time_step, [x, y])
    p_analytical = taylor_green_vortex(p['x'], p['y'], 
                                       V=parameters.amplitude, 
                                       time=parameters.time , 
                                       Re=parameters.Re, 
                                       name='pressure')
    cases[i]['p'] = p['values'][ratio-1::ratio, ratio-1::ratio]
    cases[i]['p-error'] = numpy.linalg.norm(p['values']-p_analytical)/p['values'].size

  print('\nOrders of convergence:')
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

  if parameters.save or parameters.show:
    print('\nPlot the grid convergence ...')
    pyplot.style.use('{}/scripts/python/style/'
                     'style_PetIBM.mplstyle'.format(os.environ['PETIBM_DIR']))
    pyplot.xlabel('cell-width')
    pyplot.ylabel('$L_2$-norm error')
    pyplot.plot([case['grid-spacing'] for case in cases], 
                [case['u-error'] for case in cases], 
                label='u-velocity', marker='o')
    pyplot.plot([case['grid-spacing'] for case in cases], 
                [case['v-error'] for case in cases],
                label='v-velocity', marker='o')
    pyplot.plot([case['grid-spacing'] for case in cases], 
                [case['p-error'] for case in cases], 
                label='pressure', marker='o')
    h = numpy.linspace(cases[0]['grid-spacing'], cases[-1]['grid-spacing'], 101)
    pyplot.plot(h, h, label='$1^{st}$-order convergence', color='k')
    pyplot.plot(h, h**2, label='$2^{nd}$-order convergence', 
                color='k', linestyle='--')
    pyplot.legend()
    pyplot.xscale('log')
    pyplot.yscale('log')
    if parameters.save:
      pyplot.savefig('{}/{}.png'.format(parameters.directory, parameters.output))
    if parameters.show:
      pyplot.show()


if __name__ == '__main__':
  main()