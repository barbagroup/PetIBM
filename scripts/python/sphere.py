#!/usr/bin/env python

# file: sphere.py
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: Generates coordinates of a sphere.


import argparse
import os
import math

import numpy


def read_inputs():
  """Creates the parser and parses the command-line."""
  # create parser
  parser = argparse.ArgumentParser(description='Generates a spherical body',
                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  # add arguments to parser
  parser.add_argument('--radius', '-r', dest='radius', type=float, default=0.5,
            help='radius of the sphere')
  parser.add_argument('--center', '-c', dest='center', type=float, nargs='+',
                      default=[0.0, 0.0, 0.0],
                      help='coordinates of the center of the sphere')
  parser.add_argument('-ds', dest='ds', type=float, default='0.015',
                      help='mesh-spacing')
  parser.add_argument('--name', dest='file_name', type=str, default=None,
                      help='name of the output file (without extension)')
  parser.add_argument('--save-dir', dest='save_dir', type=str, 
                      default=os.getcwd(),
                      help='directory where file will be saved')
  return parser.parse_args()

def main():
  """Generates the coordinates of a sphere and writes into file."""
  # parse the command-line
  args = read_inputs()

  # parameters
  R = args.radius
  xc, yc, zc = args.center
  h = args.ds

  if not args.file_name:
    args.file_name = 'sphere_%g' % args.ds

  n_phi = int(math.ceil(math.pi*R/h))
  phi = numpy.linspace(0.0, math.pi, n_phi)[1:-1]
  n_theta = int(math.ceil(2.0*math.pi*R/h))
  theta = numpy.linspace(0.0, 2.0*math.pi, n_theta)[1:-1]

  x = xc + R*numpy.outer(numpy.sin(phi), numpy.cos(theta)).flatten()
  x = numpy.insert(x, 0, xc)
  x = numpy.insert(x, -1, xc)
  y = yc + R*numpy.outer(numpy.sin(phi), numpy.sin(theta)).flatten()
  #y = numpy.ravel(yc + R*numpy.outer(numpy.sin(phi), numpy.sin(theta)))
  y = numpy.insert(y, 0, yc)
  y = numpy.insert(y, -1, yc)
  z = zc + R*numpy.outer(numpy.cos(phi), numpy.ones(theta.size)).flatten()
  #z = numpy.ravel(zc + R*numpy.outer(numpy.cos(phi), numpy.ones(theta.size)))
  z = numpy.insert(z, 0, R+zc)
  z = numpy.insert(z, -1, -R+zc)

  with open('%s/%s.body' % (args.save_dir, args.file_name), 'w') as outfile:
    outfile.write('%d\n' % x.size)
    numpy.savetxt(outfile, numpy.c_[x, y, z], fmt='%.6f', delimiter='\t')


if __name__=="__main__":
  main()