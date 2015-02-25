#!/usr/bin/env python

# file: plotForces.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Plot the instantaneous forces.


import os
import sys
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
  parser.add_argument('--name', dest='image_name', type=str, default='forces',
                      help='name of the .png file (without extension)')
  parser.add_argument('--no-show', dest='show', action='store_false',
                      help='does not display the figure')
  parser.add_argument('--no-save', dest='save', action='store_false',
                      help='does not save the figure')
  parser.set_defaults(show=True, save=True)
  return parser.parse_args()


class Force(object):
  """Contains info about a force."""
  def __init__(self):
    """Initializes the force."""
    self.values = numpy.empty(0)

  def mean(self, mask):
    """Computes the mean force.

    Arguments
    ---------
    mask -- range of indices to consider
    """
    self.mean = numpy.mean(self.values[mask])

class Case(object):
  """Contains info about forces of a simulation."""
  def __init__(self, parameters):
    """Reads and plots the forces.

    Arguments
    ---------
    parameters -- arguments from parser
    """
    print 'Case: %s' % parameters.case_directory
    self.parameters = parameters
    self.legend = os.path.basename(self.parameters.case_directory)
    self.read_forces()
    self.get_mean_values()
    self.plot_forces()

  def read_forces(self):
    """Reads forces from file."""
    print 'Reading forces ...'
    self.fx, self.fy, self.fz = Force(), Force(), Force()
    forces_path = '%s/forces.txt' % self.parameters.case_directory
    with open(forces_path, 'r') as infile:
      t, f_x, f_y, f_z = numpy.loadtxt(infile, dtype=float, delimiter='\t', 
                                       unpack=True)

    t = numpy.unique(t)
    mask = numpy.where(numpy.logical_and(t >= self.parameters.time_limits[0],
              t <= self.parameters.time_limits[1]))[0][::self.parameters.stride]
    self.times = t[mask]
    self.fx.values = f_x[mask]
    self.fy.values = f_y[mask]
    self.fz.values = f_z[mask]

  def get_mean_values(self):
    """Computes the averaged forces"""
    lower_limit = self.parameters.average_limits[0]
    upper_limit = self.parameters.average_limits[1]
    mask = numpy.where(numpy.logical_and(self.times >= lower_limit,
                                         self.times <= upper_limit))[0]
    self.fx.mean(mask)
    self.fy.mean(mask)
    self.fz.mean(mask)
    print '\nAveraging forces between %g and %g:' % (self.times[mask[0]],
                                                     self.times[mask[-1]])
    print '\tfx = %g' % self.fx.mean
    print '\tfy = %g' % self.fy.mean
    print '\tfz = %g' % self.fz.mean
    print '\n'

  def plot_forces(self):
    """Displays the forces into a figure."""
    if not (self.parameters.show or self.parameters.save):
      return
    print 'Plotting forces ...'
    pyplot.figure()
    pyplot.grid(True)
    pyplot.xlabel('times', fontsize=16)
    pyplot.ylabel('forces', fontsize=16)
    pyplot.plot(self.times, self.fx.values, label=r'$f_x$', 
                color='b', linestyle='-', linewidth=1)
    pyplot.plot(self.times, self.fy.values, label=r'$f_y$', 
                color='b', linestyle='--', linewidth=1)
    pyplot.plot(self.times, self.fz.values, label=r'$f_z$', 
                color='k', linestyle='-', linewidth=1)
    pyplot.legend(loc='best', prop={'size': 18}, bbox_to_anchor=(1.0, 1.0))

    if not self.parameters.plot_limits[0]:
      self.parameters.plot_limits[0] = self.times[0]
      self.parameters.plot_limits[1] = self.times[-1]
      start = 10
      minimum = min(self.fx.values[start:].min(), 
                    self.fy.values[start:].min(), 
                    self.fz.values[start:].min())
      maximum = max(self.fx.values[start:].max(), 
                    self.fy.values[start:].max(), 
                    self.fz.values[start:].max())
      self.parameters.plot_limits[2] = minimum - 0.5
      self.parameters.plot_limits[3] = maximum + 0.5
    pyplot.axis([self.parameters.plot_limits[0], self.parameters.plot_limits[1], 
                 self.parameters.plot_limits[2], self.parameters.plot_limits[3]])
    if self.parameters.save:
      images_directory = '%s/images' % self.parameters.case_directory
      if not os.path.isdir(images_directory):
        os.makedirs(images_directory)
      pyplot.savefig('%s/%s.png' % (images_directory, 
                                    self.parameters.image_name), 
                     bbox_inches='tight')
      print 'Figure saved in folder %s' % os.path.basename(images_directory)
    if self.parameters.show:
      pyplot.show()
    pyplot.close()


def main():
  """Plots the instantaneous force coefficients."""
  args = read_inputs()
  Case(args)


if __name__ == '__main__':
  main()