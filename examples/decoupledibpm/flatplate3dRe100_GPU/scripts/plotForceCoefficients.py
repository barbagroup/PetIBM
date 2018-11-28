"""
Plots the steady-state force coefficients of an inclined flat-plate with
aspect-ratio 2 at Reynolds number 100 for angles of attack between 0 and 90
degrees.
Compares with experimental results reported in Taira et al. (2007).
_References:_
* Taira, K., Dickson, W. B., Colonius,
  T., Dickinson, M. H., & Rowley, C. W. (2007).
  Unsteadiness in flow over a flat plate at angle-of-attack at low Reynolds
  numbers.
  AIAA Paper, 710, 2007.
"""

import os
import pathlib
import numpy
from matplotlib import pyplot


root_dir = pathlib.Path(__file__).absolute().parents[1]

# Read forces and computes mean values for each angle of inclination.
time_limits = (15.0, 20.0)
angles = numpy.arange(0, 90 + 1, 10, dtype=numpy.int32)
cd = numpy.zeros_like(angles, dtype=numpy.float64)
cl = numpy.zeros_like(angles, dtype=numpy.float64)
for i, angle in enumerate(angles):
    filepath = root_dir / 'AoA{}'.format(angle) / 'forces-0.txt'
    with open(filepath, 'r') as infile:
        data = numpy.loadtxt(infile, dtype=numpy.float64, unpack=True)
    mask = numpy.where(numpy.logical_and(data[0] >= time_limits[0],
                                         data[0] <= time_limits[1]))
    cd[i], cl[i] = data[1][mask].mean(), data[2][mask].mean()

# Read experimental data from Taira et al. (2007).
examples_dir = os.environ.get('PETIBM_EXAMPLES')
if not examples_dir:
    examples_dir = root_dir.parents[1]
data_dir = examples_dir / 'data'
taira = {'cd': {'aoa': None, 'values': None,
                'filename': 'taira_et_al_2007_flatPlateRe100AR2_CdvsAoA.dat'},
         'cl': {'aoa': None, 'values': None,
                'filename': 'taira_et_al_2007_flatPlateRe100AR2_ClvsAoA.dat'}}
for key in taira.keys():
    filepath = data_dir / taira[key]['filename']
    with open(filepath, 'r') as infile:
        data = numpy.loadtxt(infile, dtype=numpy.float64, unpack=True)
    taira[key]['aoa'], taira[key]['values'] = data[0], data[1]

pyplot.rc('font', family='serif', size=16)

# Plots the force coefficients versus the angle-of-attack and compares with
# experimental results reported in Taira et al. (2007).
fig, ax = pyplot.subplots(nrows=2, figsize=(6.0, 6.0), sharex=True)
ax[0].grid()
ax[0].set_ylabel('Drag coefficient')
ax[0].scatter(angles, cd,
              label='PetIBM',
              marker='x', s=40,
              facecolors='black', edgecolors='none')
ax[0].scatter(taira['cd']['aoa'], taira['cd']['values'],
              label='Taira et al. (2007)',
              marker='o', s=40,
              facecolors='none', edgecolors='#1B9E77')
ax[0].set_ylim(0.0, 2.0)
ax[1].grid()
ax[1].set_xlabel('Angle of attack (deg)')
ax[1].set_ylabel('Lift coefficient')
ax[1].scatter(angles, cl,
              label='PetIBM',
              marker='x', s=40,
              facecolors='black', edgecolors='none')
ax[1].scatter(taira['cl']['aoa'], taira['cl']['values'],
              label='Taira et al. (2007)',
              marker='o', s=40,
              facecolors='none', edgecolors='#1B9E77')
ax[1].set_xlim(0.0, 90.0)
ax[1].set_ylim(0.0, 2.0)
handles, labels = ax[0].get_legend_handles_labels()
fig.legend(handles, labels,
           ncol=2, loc='center',
           frameon=False, bbox_to_anchor=(0.54, 0.53))
fig.tight_layout()

pyplot.show()

# Save figure.
fig_dir = root_dir / 'figures'
fig_dir.mkdir(parents=True, exist_ok=True)
filepath = fig_dir / 'forceCoefficients.png'
fig.savefig(str(filepath), dpi=300)
