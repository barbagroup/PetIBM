"""
Plot the force coefficients for all immersed bodies.
"""

import pathlib
import numpy
from matplotlib import pyplot


def get_force_coefficients(filepath, data, label=None,
                           coeff=1.0, usecols=None):
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
        t, fx, fy = numpy.loadtxt(infile, dtype=numpy.float64,
                                  usecols=usecols, unpack=True)
    data[label] = {'t': t, 'cd': coeff * fx, 'cl': coeff * fy}
    return


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
simu_dir = pathlib.Path(__file__).absolute().parents[1]
data_dir = simu_dir / 'output'

# Get force coefficients for all immersed bodies
filepath = data_dir / 'forces-0.txt'
data = {}
get_force_coefficients(filepath, data, label='Body 1',
                       coeff=2.0, usecols=(0, 1, 2))
get_force_coefficients(filepath, data, label='Body 2',
                       coeff=2.0, usecols=(0, 3, 4))

# Compute the time-averaged force coefficients and min/max lift coefficient.
time_limits = (125.0, 200.0)
for label, body in data.items():
    mask = get_time_mask(body['t'], time_limits=time_limits)
    cd, cl = body['cd'], body['cl']
    cd_mean, cl_mean = cd[mask].mean(), cl[mask].mean()
    cl_min, cl_max = cl[mask].min(), cl[mask].max()
    print('{}:'.format(label))
    print('\t<Cd> = {:0.4f}'.format(cd_mean))
    print('\t<Cl> = {:0.4f} ([{:0.4f}, {:0.4f}])'
          .format(cl_mean, cl_min, cl_max))

pyplot.rc('font', family='serif', size=16)

# Plot the force coefficients.
fig, ax = pyplot.subplots(nrows=2, figsize=(10.0, 6.0), sharex=True)
ax[0].grid()
ax[0].set_ylabel('Drag coefficient')
for label, body in data.items():
    ax[0].plot(body['t'], body['cd'], label=label)
ax[0].set_ylim(1.6, 1.8)
ax[1].grid()
ax[1].set_xlabel('Non-dimensional time')
ax[1].set_ylabel('Lift coefficient')
for label, body in data.items():
    ax[1].plot(body['t'], body['cl'], label=label)
ax[1].set_xlim(0.0, 200.0)
ax[1].set_ylim(-0.6, 0.6)
handles, labels = ax[0].get_legend_handles_labels()
fig.legend(handles, labels,
           ncol=2, loc='center', frameon=False, bbox_to_anchor=(0.54, 0.53))
fig.tight_layout()

pyplot.show()

# Save figure.
fig_dir = simu_dir / 'figures'
fig_dir.mkdir(parents=True, exist_ok=True)
filepath = fig_dir / 'forceCoefficients.png'
fig.savefig(str(filepath), dpi=300)
