#!/usr/bin/env python

# file: runConvectiveTermTest.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# brief: Runs the unit-test for the convective term in 2d and 3d 
#        and plot the spatial convergence.


import os

import numpy
from matplotlib import pyplot


def plot_spatial_convergence(file_path, legend, output_file):
  """Plots the relative errors in the convective term versus the grid-spacing.

  Parameters
  ----------
  file_path: str
    Path of the file with relative errors and grid-spacings.
  legend: str
    Legend of the plot.
  output_file: str
    Path of the image to be saved.
  """
  with open(file_path, 'r') as infile:
    h, errors = numpy.loadtxt(infile, dtype=float, unpack=True)
  pyplot.figure(figsize=(6, 6))
  pyplot.grid(True)
  pyplot.xlabel('grid-spacing')
  pyplot.ylabel('relative error in convective term ($L_2$-norm)')
  pyplot.plot(h, errors, label=legend,
              marker='o', markersize=6)
  h_gauge = numpy.linspace(h[0], h[-1], 101)
  pyplot.plot(h_gauge, errors[0]*(h_gauge/h[0]),
              label='$1^{st}$-order convergence',
              color='k', ls='-')
  pyplot.plot(h_gauge, errors[0]*(h_gauge/h[0])**2,
              label='$2^{nd}$-order convergence',
              color='k', ls='--')
  pyplot.legend(loc='best', frameon=False)
  pyplot.xscale('log')
  pyplot.yscale('log')
  pyplot.savefig(output_file)

def main():
  """Builds the unit-test, runs it for 2d grids and 3d grids 
  and plots the spatial convergence for the 2d and 3d cases.
  """
  test_directory = os.getcwd()

  # plot and save spatial convergence
  pyplot.style.use('{}/scripts/python/style/'
                   'style_PetIBM.mplstyle'.format(os.environ['PETIBM_DIR']))
  # 2d cases
  plot_spatial_convergence('{}/data/relativeErrors2d.dat'.format(test_directory), 
                           'convective term (2d)', 
                           '{}/data/convergenceConvectiveTerm2d.png'.format(test_directory))
  # 3d cases
  plot_spatial_convergence('{}/data/relativeErrors3d.dat'.format(test_directory), 
                           'convective term (3d)', 
                           '{}/data/convergenceConvectiveTerm3d.png'.format(test_directory))


if __name__ == '__main__':
  main()