#!/usr/bin/env python

# file: plotForces.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Plots the instantaneous forces acting on multiple bodies.


import os
import argparse

import numpy
from scipy import signal
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
                      help='number of dimensions of the problem')
  parser.add_argument('--coefficient', dest='coefficient', type=float, 
                      default=1.0,
                      help='force to force coefficient converter')
  parser.add_argument('--times', '-t', dest='time_limits', type=float, 
                      nargs='+', default=[0.0, float('inf')],
                      help='temporal interval to consider (initial, final)')
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
  parser.add_argument('--extrema', dest='extrema', action='store_true',
                      help='displays the forces extrema')
  parser.add_argument('--order', dest='order', type=int, default=5,
                      help='number of side-points used for comparison '
                           'to get extrema')
  parser.add_argument('--no-show', dest='show', action='store_false',
                      help='does not display the figure')
  parser.add_argument('--no-save', dest='save', action='store_false',
                      help='does not save the figure')
  parser.set_defaults(drag=True, lift=True, sideforce=True, 
                      show=True, save=True)
  return parser.parse_args()


class ForceCoefficient(object):
  """Contains info about a force."""
  def __init__(self, values):
    """Initializes the force coefficient.

    Parameters
    ----------
    values: Numpy array
      Instantaneous values of the force coefficient.
    """
    self.values = values

  def get_mean(self, mask):
    """Computes the mean force coefficient.

    Parameters
    ----------
    mask: Numpy array
      Range of indices to consider.
    """
    self.mean = numpy.mean(self.values[mask])

  def get_fluctuations(self, mask):
    """Computes the fluctuations around the mean value.

    Parameters
    ----------
    mask: Numpy array
      Rane of indices to consider.
    """
    self.fluctuations = numpy.absolute([self.values[mask].min() - self.mean, 
                                        self.values[mask].max() - self.mean])

  def get_extrema(self, order=5):
    """Computes extrema indices (minima and maxima) of the force coefficient.

    Parameters
    ----------
    order: int
      Number of points on each side to use for comparison; default: 5.
    """
    minima = signal.argrelextrema(self.values, numpy.less_equal, order=order)[0][:-1]
    maxima = signal.argrelextrema(self.values, numpy.greater_equal, order=order)[0][:-1]
    # remove indices that are too close
    self.minima = minima[numpy.append(True, minima[1:]-minima[:-1] > order)]
    self.maxima = maxima[numpy.append(True, maxima[1:]-maxima[:-1] > order)]


class Body(object):
  """Contains information about forces acting on a body."""
  def __init__(self, name, forces, coefficient=1.0):
    """Initializes the force coefficients.

    Parameters
    ----------
    name: str
      Name of the body.
    forces: list(Numpy array)
      Force in the x-, y-  and z-(if applicable) directions.
    coefficient: float
      Number to convert a force into a force coefficient; default: 1.0.
    """
    self.name = name
    self.force_coefficients = []
    for force in forces:
      self.force_coefficients.append(ForceCoefficient(coefficient*force))

  def get_means(self, mask):
    """Computes the mean forces.

    Parameters
    ----------
    mask: Numpy array
      Range of indices to consider.
    """
    for i in xrange(len(self.force_coefficients)):
      self.force_coefficients[i].get_mean(mask)

  def get_fluctuations(self, mask):
    """Computes the fluctuations around the mean force over a given range.

    Parameters
    ----------
    mask: Numpy array
      Range of indices to consider.
    """
    for i in xrange(len(self.force_coefficients)):
      self.force_coefficients[i].get_fluctuations(mask)

  def get_extrema(self, order=5):
    """Computes the extrema indices of each forces acting on the body.

    Parameters
    ----------
    times: Numpy array
      Discrete time values.
    mask: Numpy array
      Range of indices to consider.
    order: float
      Parameters used for comparison to define an extremum; default: 5.
    """
    for i in xrange(len(self.force_coefficients)):
      self.force_coefficients[i].get_extrema(order=order)

  def get_strouhal(self, times):
    """Computes the Strouhal number based on the time-variational lift 
    (using the minima of the lift).

    Parameters
    ----------
    times: Numpy array
      Discrete time values.
    """
    cl = self.force_coefficients[1]
    if times[cl.minima[-1]] > times[cl.maxima[-1]]:
      if cl.minima.size > 1:
        self.strouhal = 1.0/(times[cl.minima[-1]]-times[cl.minima[-2]])
      else:
        self.strouhal = float('inf')
    else:
      if cl.maxima.size > 1:
        self.strouhal = 1.0/(times[cl.maxima[-1]]-times[cl.maxima[-2]])
      else:
        self.strouhal = float('inf')

  def print_info(self):
    """Prints the time-averaged forces and force-coefficients if necessary."""
    print('\nBody: {}'.format(self.name))
    print('\tSt = {:0.3f}'.format(self.strouhal))
    names = ['cd', 'cl', 'cm']
    for i, coefficient in enumerate(self.force_coefficients):
      print('\t<{}> = {:0.3f} \t[{:0.4f}, {:0.4f}]'.format(names[i],
                                                           coefficient.mean,
                                                           coefficient.fluctuations[0],
                                                           coefficient.fluctuations[1]))


class Case(object):
  """Contains info about forces of a simulation."""
  def __init__(self, parameters):
    """Stores the parameters.

    Parameters
    ----------
    parameters: ArgumentPaser instance
      Contains command-line arguments.
    """
    print('\nCase: {}'.format(parameters.case_directory))
    self.parameters = parameters

  def read_forces(self):
    """Reads forces from file."""
    print('\nReading forces ...')
    forces_path = '{}/forces.txt'.format(self.parameters.case_directory)
    with open(forces_path, 'r') as infile:
      data = numpy.loadtxt(infile, dtype=float).transpose()
    t, forces = data[0], data[1:]
    n_bodies = forces.shape[0]/self.parameters.dimensions
    # keep unique values between time-limits of consideration
    mask = numpy.where(numpy.logical_and(t >= self.parameters.time_limits[0],
                                         t <= self.parameters.time_limits[1]))[0]
    _, mask_uniqueness = numpy.unique(t, return_index=True)
    mask = list(set(mask) & set(mask_uniqueness))
    self.times = t[mask]
    # give default name to bodies if necessary
    if not self.parameters.body_names:
      self.parameters.body_names = ['body {}'.format(i+1) for i in xrange(n_bodies)]
    # create Body objects that store the forces of each body
    print('\nInfo: ForceCoefficient = {} x Force'.format(self.parameters.coefficient))
    print('If not correct, use command-line argument --coefficient')
    if self.parameters.dimensions == 3:
      self.bodies = [Body(self.parameters.body_names[i], 
                          [forces[3*i][mask], forces[3*i+1][mask], forces[3*i+2][mask]],
                          coefficient=self.parameters.coefficient)
                     for i in xrange(n_bodies)]
    else:
      self.bodies = [Body(self.parameters.body_names[i], 
                          [forces[2*i][mask], forces[2*i+1][mask]],
                          coefficient=self.parameters.coefficient) 
                     for i in xrange(n_bodies)]

  def get_properties(self):
    """Computes various properties such as the average, the fluctuations,
    the Stouhal number and the extrema indices, for each body."""
    # neglect high value at beginning
    if self.parameters.average_limits[0] == 0.0:
      self.parameters.average_limits[0] = self.times[1]
    # get mask
    lower_limit, upper_limit = self.parameters.average_limits
    mask = numpy.where(numpy.logical_and(self.times >= lower_limit,
                                         self.times <= upper_limit))[0]
    # computes proerties for each body
    print('\nConsidering force coefficients between {} and {}:'.format(self.times[mask[0]],
                                                                       self.times[mask[-1]]))
    for i, body in enumerate(self.bodies):
      self.bodies[i].get_means(mask)
      self.bodies[i].get_fluctuations(mask)
      self.bodies[i].get_extrema(order=self.parameters.order)
      self.bodies[i].get_strouhal(self.times)
      self.bodies[i].print_info()

  def plot_force_coefficients(self):
    """Displays the force coefficients into a figure."""
    print('\nPlotting forces ...')
    pyplot.style.use('{}/style/'
                     'style_PetIBM.mplstyle'.format(os.path.dirname(__file__)))
    fig, ax = pyplot.subplots(figsize=(8, 6))
    color_cycle = ax._get_lines.color_cycle
    pyplot.grid(True, zorder=0)
    pyplot.xlabel('time')
    pyplot.ylabel('force coefficients')
    # plot forces for each immersed body
    for i, body in enumerate(self.bodies):
      color = next(color_cycle)
      if self.parameters.drag:
        cd = body.force_coefficients[0]
        pyplot.plot(self.times, cd.values,  
                    label='{} - $C_D$'.format(body.name),
                    color=color, linestyle='-', zorder=10)
        if self.parameters.extrema:
          pyplot.scatter(self.times[cd.minima], cd.values[cd.minima],
                         c=color, marker='o', zorder=10)
          pyplot.scatter(self.times[cd.maxima], cd.values[cd.maxima],
                         c=color, marker='o', zorder=10)
      if self.parameters.lift:
        cl = body.force_coefficients[1]
        pyplot.plot(self.times, cl.values,
                    label='{} - $C_L$'.format(body.name), 
                    color=color, linestyle='--', zorder=10)
        if self.parameters.extrema:
          pyplot.scatter(self.times[cl.minima], cl.values[cl.minima],
                         c=color, marker='s', zorder=10)
          pyplot.scatter(self.times[cl.maxima], cl.values[cl.maxima],
                         c=color, marker='s', zorder=10)
      if self.parameters.dimensions == 3 and self.parameters.sideforce:
        cz = body.force_coefficients[2]
        pyplot.plot(self.times, cz.values, 
                    label='{} - $C_Z$'.format(body.name), 
                    color=color, linestyle=':', zorder=10)
        if self.parameters.extrema:
          pyplot.scatter(self.times[cz.minima], cz.values[cz.minima],
                         c=color, marker='*', zorder=10)
          pyplot.scatter(self.times[cz.maxima], cz.values[cz.maxima],
                         c=color, marker='*', zorder=10)
    pyplot.legend()
    # get default plot-limits if not specified
    if not any(self.parameters.plot_limits):
      self.parameters.plot_limits[0] = self.times[0]
      self.parameters.plot_limits[1] = self.times[-1]
      min_fx, max_fx = float('inf'), float('-inf')
      minima, maxima = [], []
      for body in self.bodies:
        for coefficient in body.force_coefficients:
          minima.append(coefficient.values[2:].min())
          maxima.append(coefficient.values[2:].max())
      minimum, maximum = min(minima), max(maxima)
      self.parameters.plot_limits[2] = minimum - 0.1*abs(minimum)
      self.parameters.plot_limits[3] = maximum + 0.1*abs(maximum)
    pyplot.axis([self.parameters.plot_limits[0], self.parameters.plot_limits[1], 
                 self.parameters.plot_limits[2], self.parameters.plot_limits[3]])
    # save the image as .png file
    if self.parameters.save:
      images_directory = '{}/images'.format(self.parameters.case_directory)
      if not os.path.isdir(images_directory):
        os.makedirs(images_directory)
      pyplot.savefig('{}/{}.png'.format(images_directory, self.parameters.image_name))
      print('\nFigure saved in folder {}'.format(os.path.basename(images_directory)))
    # dislay the figure
    if self.parameters.show:
      pyplot.show()
    pyplot.close()


def main():
  """Plots the instantaneous force coefficients."""
  args = read_inputs()
  case = Case(args)
  case.read_forces()
  case.get_properties()
  if args.show or args.save:
    case.plot_force_coefficients()


if __name__ == '__main__':
  print('\n[{}] START\n'.format(os.path.basename(__file__)))
  main()
  print('\n[{}] END\n'.format(os.path.basename(__file__)))