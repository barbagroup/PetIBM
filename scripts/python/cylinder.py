#!/usr/bin/env python

# file: cylinder.py
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: Generates coordinates of a 3D cylinder.


import argparse
import os

import numpy


def read_inputs():
  """Creates the parser and parses the command-line."""
  # create parser
  parser = argparse.ArgumentParser(description='Generates a cylindrical body',
                        formatter_class= argparse.ArgumentDefaultsHelpFormatter)
  # add arguments to parser
  parser.add_argument('--radius', '-r', dest='radius', type=float, default=0.5,
                      help='radius of the cylinder')
  parser.add_argument('--limits', '-l', dest='limits', type=float, nargs='+',
                      default=[-3.0, 3.0],
                      help='longitudinal limits of the cylinder')
  parser.add_argument('--center', '-c', dest='center', type=float, nargs='+',
                      default=[0.0, 0.0],
                      help='coordinates of the center of the cylinder')
  parser.add_argument('-ds', dest='ds', type=float, default='0.015',
                      help='mesh-spacing')
  parser.add_argument('--name', dest='file_name', type=str, default=None,
                      help='name of the output file (without extension)')
  parser.add_argument('--save-dir', dest='save_dir', type=str, 
                      default=os.getcwd(),
                      help='directory where file will be saved')
  return parser.parse_args()

def main():
  """Generates the coordinates of a 3D cylinder and writes into file."""
  # parse command-line
  args = read_inputs()
  
  R = args.radius
  zmin = args.limits[0]
  zmax = args.limits[1]
  ds = args.ds
  
  length = zmax-zmin
  nz = int(round(length/ds))
  if abs(nz-length/ds) > 1.0E-08:
    print ('Choose a mesh spacing such that the cylinder length is an '
           'integral multiple of it')
    print ('%s-direction: length l=%g \t spacing ds=%g \t l/ds=%g' 
           % (length, ds, length/ds))
    sys.exit()

  nc = int(numpy.ceil(2.0*numpy.pi*R/ds))

  if not args.file_name:
    args.file_name = 'cylinder_%g.body' % args.ds

  with open(args.save_dir+'/'+args.file_name, 'w') as outfile:
    outfile.write('%d\n' % (nc*nz))
    i = numpy.arange(nc)
    x, y = R*numpy.cos(2.0*numpy.pi*i/nc), R*numpy.sin(2.0*numpy.pi*i/nc)
    for k in xrange(nz):
      z = zmin + (k+0.5)*ds * numpy.ones(nc)
      numpy.savetxt(outfile, numpy.c_[x, y, z], fmt='%.6f', delimiter='\t')


if __name__ == "__main__":
  main()