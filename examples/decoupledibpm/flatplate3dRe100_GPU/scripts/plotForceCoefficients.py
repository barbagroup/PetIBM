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
import numpy
from matplotlib import pyplot


if not os.environ.get('PETIBM_EXAMPLES'):
  raise KeyError('Environment variable PETIBM_EXAMPLES is not set; '
                 'Set PETIBM_EXAMPLES as the root directory of the examples.')

script_dir = os.path.dirname(os.path.realpath(__file__))
root_dir = os.sep.join(script_dir.split(os.sep)[:-1])

# Read forces and computes mean values for each angle of inclination.
time_limits = (15.0, 20.0)
angles = numpy.arange(0, 90 + 1, 10, dtype=numpy.int32)
cd = numpy.zeros_like(angles, dtype=numpy.float64)
cl = numpy.zeros_like(angles, dtype=numpy.float64)
for i, angle in enumerate(angles):
  filepath = os.path.join(root_dir, 'AoA{}'.format(angle), 'forces.txt')
  with open(filepath, 'r') as infile:
    data = numpy.loadtxt(infile, dtype=numpy.float64, unpack=True)
  mask = numpy.where(numpy.logical_and(data[0] >= time_limits[0],
                                       data[0] <= time_limits[1]))[0]
  cd[i], cl[i] = data[1][mask].mean(), data[2][mask].mean()

# Read experimental data from Taira et al. (2007).
directory = os.path.join(os.environ['PETIBM_EXAMPLES'], 'data')
taira = {'cd': {'aoa': None, 'values': None,
                'filename': 'taira_et_al_2007_flatPlateRe100AR2_CdvsAoA.dat'},
         'cl': {'aoa': None, 'values': None,
                'filename': 'taira_et_al_2007_flatPlateRe100AR2_ClvsAoA.dat'}}
for key in taira.keys():
  filepath = os.path.join(directory, taira[key]['filename'])
  with open(filepath, 'r') as infile:
    data = numpy.loadtxt(infile, dtype=numpy.float64, unpack=True)
  taira[key]['aoa'], taira[key]['values'] = data[0], data[1]

# Plots the force coefficients versus the angle-of-attack and compares with
# experimental results reported in Taira et al. (2007).
pyplot.style.use('seaborn-dark')
fig, ax = pyplot.subplots(2, figsize=(6.0, 6.0), sharex=True)
ax[0].grid(zorder=0)
ax[0].set_ylabel('$C_D$', fontname='DejaVu Serif', fontsize=16)
ax[0].scatter(angles, cd,
              label='PetIBM',
              marker='x', s=40,
              facecolors='black', edgecolors='none',
              zorder=4)
ax[0].scatter(taira['cd']['aoa'], taira['cd']['values'],
              label='Taira et al. (2007)',
              marker='o', s=40,
              facecolors='none', edgecolors='#1B9E77',
              zorder=3)
ax[0].set_ylim(0.0, 2.0)
ax[1].grid(zorder=0)
ax[1].set_xlabel('angle of attack (deg)',
                 fontname='DejaVu Serif', fontsize=16)
ax[1].set_ylabel('$C_L$', fontname='DejaVu Serif', fontsize=16)
ax[1].scatter(angles, cl,
              label='PetIBM',
              marker='x', s=40,
              facecolors='black', edgecolors='none',
              zorder=4)
ax[1].scatter(taira['cl']['aoa'], taira['cl']['values'],
              label='Taira et al. (2007)',
              marker='o', s=40,
              facecolors='none', edgecolors='#1B9E77',
              zorder=3)
ax[1].set_xlim(0.0, 90.0)
ax[1].set_ylim(0.0, 2.0)
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

figures_dir = os.path.join(root_dir, 'figures')
if not os.path.isdir(figures_dir):
  os.makedirs(figures_dir)
filepath = os.path.join(figures_dir, 'forceCoefficients.png')
fig.savefig(filepath)

pyplot.show()
