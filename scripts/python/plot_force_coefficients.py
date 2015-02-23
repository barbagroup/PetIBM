#!/usr/bin/env python

# file: plot_force_coefficients.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Plot the instantaneous force coefficients.


import os
import sys
import argparse

import numpy
from matplotlib import pyplot


def read_inputs():
  """Parses the command-line."""
  # create parser
  parser = argparse.ArgumentParser(description='Plots the instantaneous '
                                               'force coefficients',
                        formatter_class= argparse.ArgumentDefaultsHelpFormatter)
  # fill parser with arguments
  parser.add_argument('--case', dest='case_directory', type=str, 
                      default=os.getcwd(), help='directory of the simulation')
  parser.add_argument('--times', '-t', dest='times', type=float, nargs='+',
                      default=[0.0, float('inf'), None],
                      help='times to display (initial, final, stride)')
  return parser.parse_args()


class ForceCoefficient(object):
  """Contains info about a force coefficient."""
  def __init__(self, name):
    self.name = name
    self.values = numpy.empty(0)

class Case(object):
  def __init__(self, directory):
    self.directory = directory


def main():
  """Plots the instantaneous force coefficients."""


if __name__ == '__main__':
  main()