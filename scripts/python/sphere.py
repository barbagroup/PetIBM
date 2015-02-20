#!/usr/bin/env python

# file: sphere.py
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: Generates coordinates of a sphere.


import argparse
import os

import numpy


def read_inputs():
  """Creates the parser and parses the command-line."""
  # create parser
  parser = argparse.ArgumentParser(description='Generates a spherical body',
                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  # add arguments to parser
  parser.add_argument('--radius', '-r', dest='radius', type=float, default=0.5,
            help='radius of the sphere')
  parser.add_argument('--center', '-c', dest=center, type=float, nargs='+',
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

  # `total` keeps track of the number of points on the body
  # it is initialized with 2 to count the points on the poles
  total = 2

  # calculate angle increments in altitude and azimuth
  nCircHalf = int(np.pi*R/h)+1
  dPhi = np.pi/nCircHalf
  for i in range(1,nCircHalf):
    phi = i*dPhi
    r = R*np.sin(phi)
    nCircAzim = int(2*np.pi*r/h)+1
    total += nCircAzim

  # write body data to file
  f = open(filename, 'w')
  f.write("%d\n" % total)
  f.write("%f\t%f\t%f\n" % (x0, y0, R+z0))
  for i in range(1,nCircHalf):
    phi = i*dPhi
    z = R*np.cos(phi)
    r = R*np.sin(phi)
    nCircAzim = int(2*np.pi*r/h)+1
    dTheta = 2*np.pi/nCircAzim
    for j in range(nCircAzim):
      theta = j*dTheta
      x = r*np.cos(theta)
      y = r*np.sin(theta)
      f.write("%f\t%f\t%f\n" % (x+x0, y+y0, z+z0))
  f.write("%f\t%f\t%f\n" % (x0, y0, -R+z0))
  f.close()


if __name__=="__main__":
  main()