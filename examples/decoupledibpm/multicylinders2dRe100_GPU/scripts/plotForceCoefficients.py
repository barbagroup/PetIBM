"""
Plot the force coefficients for all immersed bodies.
"""

import os
import numpy
from matplotlib import pyplot


def get_force_coefficients(filepath, label=None, coeff=1.0, usecols=None):
  """
  Returns the force coefficients after reading the forces from file.

  Parameters
  ----------
    filepath: string
      Path of the file to read.
    label: string, optional
      Name of the body;
      default: None.
    coeff: float, optional
      Coefficient to convert force into force coefficient;
      default: 1.0.
    usecols: tuple of integers, optional
      Columns to read in the file;
      default: None.

  Returns
  -------
    dict: dictionary
      Contains the label of the body, the time values,
      and the force coefficients.
  """
  with open(filepath, 'r') as infile:
    data = numpy.loadtxt(infile,
                         dtype=numpy.float64, usecols=usecols, unpack=True)
  return {'label': label, 't': data[0],
          'cd': coeff * data[1], 'cl': coeff * data[2]}


def get_time_mask(data, time_limits=(-numpy.inf, numpy.inf)):
  """
  Get a mask given the time limits.

  Parameters
  ----------
    data: 1D Numpy array of floats
      The time values

  Returns
  -------
    mask: 1D Numpy array of integers
      The mask
  """
  mask = numpy.where(numpy.logical_and(time_limits[0] <= data,
                                       data <= time_limits[1]))[0]
  return mask


# Set up root directory
script_dir = os.path.dirname(os.path.realpath(__file__))
root_dir = os.sep.join(script_dir.split(os.sep)[:-1])

# Get force coefficients for all immersed bodies
filepath = os.path.join(root_dir, 'forces.txt')
bodies = []
bodies.append(get_force_coefficients(filepath, label='Body 1',
                                     coeff=2.0, usecols=(0, 1, 2)))
bodies.append(get_force_coefficients(filepath, label='Body 2',
                                     coeff=2.0, usecols=(0, 3, 4)))

# Compute the time-averaged force coefficients and min/max lift coefficient.
time_limits = (125.0, 200.0)
for body in bodies:
  mask = get_time_mask(body['t'], time_limits=time_limits)
  cd, cl = body['cd'], body['cl']
  cd_mean, cl_mean = cd[mask].mean(), cl[mask].mean()
  cl_min, cl_max = cl[mask].min(), cl[mask].max()
  print('{}:'.format(body['label']))
  print('\t<Cd> = {:0.4f}'.format(cd_mean))
  print('\t<Cl> = {:0.4f} ([{:0.4f}, {:0.4f}])'
        .format(cl_mean, cl_min, cl_max))

# Plot the force coefficients.
pyplot.style.use('seaborn-dark')
fig, ax = pyplot.subplots(2, figsize=(10.0, 6.0), sharex=True)
ax[0].grid(zorder=0)
ax[0].set_ylabel('$C_D$', fontname='DejaVu Serif', fontsize=16)
for body in bodies:
  ax[0].plot(body['t'], body['cd'],
             label=body['label'], linewidth=1, zorder=10)
ax[0].set_ylim(1.6, 1.8)
ax[1].grid(zorder=0)
ax[1].set_xlabel('non-dimensional time unit',
                 fontname='DejaVu Serif', fontsize=16)
ax[1].set_ylabel('$C_L$', fontname='DejaVu Serif', fontsize=16)
for body in bodies:
  ax[1].plot(body['t'], body['cl'],
             label=body['label'], linewidth=1, zorder=10)
ax[1].set_xlim(0.0, 200.0)
ax[1].set_ylim(-0.6, 0.6)
for a in ax:
  for method in ['get_xticklabels', 'get_yticklabels']:
    for label in getattr(a, method)():
      label.set_fontname('DejaVu Serif')
      label.set_fontsize(14)
handles, labels = ax[0].get_legend_handles_labels()
fig.legend(handles, labels,
           ncol=2, loc='center', prop={'family': 'serif', 'size': 14},
           frameon=False, bbox_to_anchor=(0.54, 0.53))
fig.tight_layout()

# Save the figure.
figures_dir = os.path.join(root_dir, 'figures')
if not os.path.isdir(figures_dir):
  os.makedirs(figures_dir)
filepath = os.path.join(figures_dir, 'forceCoefficients.png')
fig.savefig(filepath)

pyplot.show()
