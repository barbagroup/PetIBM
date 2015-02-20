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

  # phi between 0 and \pi, poles are separately calculated
  n_phi = int(math.ceil(math.pi*R/h))+1
  phi = numpy.linspace(0.0, math.pi, n_phi)[1:-1]

  # north pole
  x, y, z = xc, yc, R+zc
  for phi in phi:
    # theta between 0 and 2\pi
    n_theta = int(math.ceil(2.0*math.pi*R*math.sin(phi)/h))+1
    theta = numpy.linspace(0.0, 2.0*math.pi, n_theta)[:-1]
    x = numpy.append(x, xc + R*math.sin(phi)*numpy.cos(theta))
    y = numpy.append(y, yc + R*math.sin(phi)*numpy.sin(theta))
    z = numpy.append(z, zc + R*math.cos(phi)*numpy.ones(theta.size))
  # south pole
  x, y, z = numpy.append(x, xc), numpy.append(y, yc), numpy.append(z, -R+zc)

  # write coordinates into file
  if not args.file_name:
    args.file_name = 'sphere_%g' % args.ds
  with open('%s/%s.body' % (args.save_dir, args.file_name), 'w') as outfile:
    outfile.write('%d\n' % x.size)
    numpy.savetxt(outfile, numpy.c_[x, y, z], fmt='%.6f', delimiter='\t')


if __name__=="__main__":
  main()