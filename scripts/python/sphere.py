#!/usr/bin/env python

import argparse
import numpy as np

def read_inputs():
	"""Creates the parser and parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='Generates a spherical body ',
									 formatter_class= \
									 argparse.ArgumentDefaultsHelpFormatter)
	# add arguments to the parser
	parser.add_argument('-R', dest='R', type=float, default=0.5,
						help='radius of the sphere')
	parser.add_argument('-x0', dest='x0', type=float, default=0.0,
						help='x-coordinate of the center of the sphere')
	parser.add_argument('-y0', dest='y0', type=float, default=0.0,
						help='y-coordinate of the center of the sphere')
	parser.add_argument('-z0', dest='z0', type=float, default=0.0,
						help='z-coordinate of the center of the sphere')
	parser.add_argument('-ds', dest='ds', type=str, default='0.015',
						help='mesh spacing')
	parser.add_argument('-f', '--filename', dest='filename', type=str,
						default='',
						help='name of the file generated')
	return parser.parse_args()

def main():
	"""Generates a sphere."""
	# parse the command-line
	args = read_inputs()

	# parameters of the sphere
	R = args.R
	x0, y0, z0 = args.x0, args.y0, args.z0
	h = float(args.ds)
	if args.filename:
		filename = args.filename
	else:
		filename = 'sphere_' + args.ds + '.body'

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
