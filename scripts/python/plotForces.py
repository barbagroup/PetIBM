#!/usr/bin/env python

# file: plotForces.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Plots the instantaneous forces acting on multiple bodies.


import os
import argparse

import numpy
from matplotlib import pyplot


def read_inputs():
  """Parses the command-line."""
  # create parser
  parser = argparse.ArgumentParser(description='Plots the instantaneous forces',
                        formatter_class= argparse.ArgumentDefaultsHelpFormatter)
  # fill parser with arguments
  parser.add_argument('--case', dest='case_directory', type=str, 
                      default=os.getcwd(), help='directory of the simulation')
  parser.add_argument('--dimensions', '-d', dest='dimensions', type=int,
                      default=2,
                      help='number of dimensions')
  parser.add_argument('--times', '-t', dest='time_limits', type=float, 
                      nargs='+', default=[0.0, float('inf')],
                      help='temporal interval to consider (initial, final)')
  parser.add_argument('--stride', '-s', dest='stride', type=int, default=1,
                      help='stride to read the values of time and forces')
  parser.add_argument('--average', '-a', dest='average_limits', type=float, 
                      nargs='+', default=[0.0, float('inf')],
                      help='temporal limits to consider to average forces')
  parser.add_argument('--limits', '-l', dest='plot_limits', type=float, 
                      nargs='+', default=[None, None, None, None],
                      help='limits of the plot')
  parser.add_argument('--bodies', dest='body_names', type=str, nargs='+',
                      help='name of each immersed body')
  parser.add_argument('--output', dest='image_name', type=str, default='forces',
                      help='name of the .png file (without extension)')
  parser.add_argument('--no-drag', dest='drag', action='store_false',
                      help='does not display the force in the x-direction')
  parser.add_argument('--no-lift', dest='lift', action='store_false',
                      help='does not display the force in the y-direction')
  parser.add_argument('--no-sideforce', dest='sideforce', action='store_false',
                      help='does not display the force in the z-direction')
  parser.add_argument('--no-show', dest='show', action='store_false',
                      help='does not display the figure')
  parser.add_argument('--no-save', dest='save', action='store_false',
                      help='does not save the figure')
  parser.set_defaults(drag=True, lift=True, sideforce=True, 
                      show=True, save=True)
  return parser.parse_args()


class Force(object):
  """Contains info about a force."""
  def __init__(self, f):
    """Initializes the force."""
    self.values = f
    self.min, self.max = f.min(), f.max()

  def mean(self, mask):
    """Computes the mean force.

    Arguments
    ---------
    mask -- range of indices to consider
    """
    self.mean = numpy.mean(self.values[mask])


class Body(object):
  """Contains forces applied to one immersed boundary."""
  def __init__(self, name, fx, fy, fz=numpy.empty(0)):
    """Initializes the forces.

    Arguments
    ---------
    name - name of the body
    fx, fy -- forces in the x- and y- directions
    """
    self.name = name
    self.fx, self.fy = Force(fx), Force(fy)
    if fz.any():
      self.fz = Force(fz)

  def means(self, mask):
    """Computes the mean forces.

    Arguments
    ---------
    mask -- range of indices to consider
    """
    self.fx.mean(mask)
    self.fy.mean(mask)
    print 'Body: %s' % self.name
    print '\t-> <fx> = %.6f' % self.fx.mean
    print '\t-> <fy> = %.6f' % self.fy.mean
    if hasattr(self, 'fz'):
      self.fz.mean(mask)
      print '\t-> <fz> = %.6f' % self.fz.mean


class Case(object):
  """Contains info about forces of a simulation."""
  def __init__(self, parameters):
    """Stores the parameters.

    Arguments
    ---------
    parameters -- arguments from parser
    """
    print '\nCase: %s' % parameters.case_directory
    self.parameters = parameters

  def read_forces(self):
    """Reads forces from file."""
    print '\nReading forces ...'
    forces_path = '%s/forces.txt' % self.parameters.case_directory
    with open(forces_path, 'r') as infile:
      data = numpy.loadtxt(infile, dtype=float).transpose()
    t, forces = data[0], data[1:]
    # create time-limits mask and unicity mask
    if self.parameters.time_limits[0] == 0.0:
      self.parameters.time_limits[0] = t[1]
    mask = numpy.where(numpy.logical_and(t >= self.parameters.time_limits[0],
              t <= self.parameters.time_limits[1]))[0][::self.parameters.stride]
    _, mask_unique = numpy.unique(t, return_index=True)
    mask = list(set(mask) & set(mask_unique))
    self.times = t[mask]
    n_bodies = forces.shape[0]/self.parameters.dimensions
    # give default name to bodies
    if not self.parameters.body_names:
      self.parameters.body_names = ['body %d' % i for i in xrange(n_bodies)]
    # create Body objects that store the forces
    if self.parameters.dimensions == 3:
      self.bodies = [Body(self.parameters.body_names[i], 
                          forces[3*i][mask], forces[3*i+1][mask], forces[3*i+2][mask])
                     for i in xrange(n_bodies)]
    else:
      self.bodies = [Body(self.parameters.body_names[i], 
                          forces[2*i][mask], forces[2*i+1][mask]) 
                     for i in xrange(n_bodies)]

  def get_mean_values(self):
    """Computes the averaged forces"""
    # create average-limits mask
    lower_limit = self.parameters.average_limits[0]
    upper_limit = self.parameters.average_limits[1]
    mask = numpy.where(numpy.logical_and(self.times >= lower_limit,
                                         self.times <= upper_limit))[0]
    # compute means for each immersed body
    print '\nAveraging forces between %g and %g:' % (self.times[mask[0]],
                                                     self.times[mask[-1]])
    for i in xrange(len(self.bodies)):
      self.bodies[i].means(mask)

  def plot_forces(self):
    """Displays the forces into a figure."""
    print '\nPlotting forces ...'
    pyplot.figure()
    pyplot.grid(True)
    pyplot.xlabel('times', fontsize=16)
    pyplot.ylabel('forces', fontsize=16)
    # only support 4 immersed bodies currently
    colors = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a']
    # plot forces for each immersed body
    for i, body in enumerate(self.bodies):
      if self.parameters.drag:
        pyplot.plot(self.times, body.fx.values, label='%s - fx' % body.name,
                    color=colors[i], linestyle='-', linewidth=2)
      if self.parameters.lift:
        pyplot.plot(self.times, body.fy.values, label='%s - fy' % body.name,
                    color=colors[i], linestyle='--', linewidth=2)
      if self.parameters.dimensions == 3 and self.parameters.sideforce:
        pyplot.plot(self.times, body.fz.values, label='%s - fz' % body.name,
                    color=colors[i], linestyle=':', linewidth=2)
    pyplot.legend(loc='best', prop={'size': 18}, bbox_to_anchor=(1.0, 1.0))
    # get default plot-limits if not specified
    if not any(self.parameters.plot_limits):
      self.parameters.plot_limits[0] = self.times[0]
      self.parameters.plot_limits[1] = self.times[-1]
      min_fx, max_fx = float('inf'), float('-inf')
      if self.parameters.drag:
        min_fx = min([body.fx.min for body in self.bodies])
        max_fx = max([body.fx.max for body in self.bodies])
      min_fy, max_fy = float('inf'), float('-inf')
      if self.parameters.lift:
        min_fy = min([body.fy.min for body in self.bodies])
        max_fy = max([body.fy.max for body in self.bodies])
      min_fz, max_fz = float('inf'), float('-inf')
      if self.parameters.dimensions == 3 and self.parameters.sideforce:
        min_fz = min([body.fz.min for body in self.bodies])
        max_fz = max([body.fz.max for body in self.bodies])
      minimum, maximum = min(min_fx, min_fy, min_fz), max(max_fx, max_fy, max_fz)
      self.parameters.plot_limits[2] = minimum - 0.1*abs(minimum)
      self.parameters.plot_limits[3] = maximum + 0.1*abs(maximum)
    pyplot.axis([self.parameters.plot_limits[0], self.parameters.plot_limits[1], 
                 self.parameters.plot_limits[2], self.parameters.plot_limits[3]])
    # save the image as .png file
    if self.parameters.save:
      images_directory = '%s/images' % self.parameters.case_directory
      if not os.path.isdir(images_directory):
        os.makedirs(images_directory)
      pyplot.savefig('%s/%s.png' % (images_directory, 
                                    self.parameters.image_name), 
                     bbox_inches='tight')
      print '\nFigure saved in folder %s' % os.path.basename(images_directory)
    # dislay the figure
    if self.parameters.show:
      pyplot.show()
    pyplot.close()


def main():
  """Plots the instantaneous force coefficients."""
  parameters = read_inputs()
  case = Case(parameters)
  case.read_forces()
  case.get_mean_values()
  if parameters.show or parameters.save:
    case.plot_forces()


if __name__ == '__main__':
  main()