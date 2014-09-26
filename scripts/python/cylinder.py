#!/usr/bin/env python

import argparse
import numpy as np

def read_inputs():
	"""Creates the parser and parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='Generates a cylindrical body ',
									 formatter_class= \
									 argparse.ArgumentDefaultsHelpFormatter)
	# add arguments to the parser
	parser.add_argument('-R', dest='R', type=float, default=0.5,
						help='radius of the cylinder')
	parser.add_argument('-zmin', dest='zmin', type=float, default=-3.,
						help='minimum z-coordinate')
	parser.add_argument('-zmax', dest='zmax', type=float, default=3.,
						help='maximum z-coordinate')
	parser.add_argument('-x0', dest='x0', type=float, default=0.0,
						help='x-coordinate of the center of the cylinder')
	parser.add_argument('-y0', dest='y0', type=float, default=0.0,
						help='y-coordinate of the center of the cylinder')
	parser.add_argument('-ds', dest='ds', type=float, default=0.015,
						help='mesh spacing')
	parser.add_argument('-f', '--filename', dest='filename', type=str,
						default='cylinder_0.015.body',
						help='name of the file generated')
	return parser.parse_args()

def main():
	args = read_inputs()
	R = args.R
	zmax = args.zmax
	zmin = args.zmin
	ds = args.ds

	nz = int(round((zmax-zmin)/ds))
	if np.abs((zmax-zmin)/ds - nz) > 1e-8:
		print "Choose a mesh spacing such that the cylinder height is an integral multiple of it!"
		print "z-direction: %s/%s = %s" % (str(zmax-zmin), str(ds), str((zmax-zmin)/ds))
	
	nb = int(np.ceil(2*np.pi*R/ds))
	total = nb*nz
	f = open(args.filename, 'w')
	f.write("%d\n" % total)
	for k in range(nz):
		z = zmin+(k+0.5)*ds
		for i in range(nb):
			x = R*np.cos(i*2*np.pi/nb)
			y = R*np.sin(i*2*np.pi/nb)
			f.write("%f\t%f\t%f\n" % (x,y,z))
	f.close()

if __name__ == "__main__":
	main()