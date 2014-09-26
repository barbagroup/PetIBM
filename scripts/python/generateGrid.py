#!/usr/bin/env python

import argparse
import numpy as np
import sys
import os
import re

def read_inputs():
	# create parser
	parser = argparse.ArgumentParser(description="Generates a cartesian mesh with a uniform region surrounded by a stretched grid", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	
	# add arguments to parser
	parser.add_argument("-output", dest="output", help="name of file generated", default="domain.yaml")
	parser.add_argument("-input", dest="input", help="name of input file", default="gridOptions")

	return parser.parse_args()

def calculate_ratios(axis, min, max, u_min, u_max, h, maxAR, g):
	precision = 2
	totalCells = 0

	g.write("- direction: %s\n" % axis)
	g.write("  start: %s\n" % str(min))
	g.write("  subDomains:\n")

	L = u_min - min
	current_precision = 1
	next_ratio = 2.
	while current_precision <= precision:
		ratio = next_ratio 
		n = int(round(np.log(L*(ratio-1)/h + 1.)/np.log(ratio)))
		AR = ratio**(n-1)
		if AR < maxAR:
			next_ratio += (0.1)**current_precision
			current_precision+=1
		else:
			next_ratio -= (0.1)**current_precision
	totalCells += n

	g.write("    - end: %s\n" % str(u_min))
	g.write("      cells: %d\n" % n)
	g.write("      stretchRatio: %s\n" % str(1.0/ratio))

	n = int(round((u_max-u_min)/h))
	if np.abs((u_max-u_min)/h - n) > 1e-8:
		print "Choose a mesh spacing such that the uniform region is an integral multiple of it!"
		print "%s-direction: %s/%s = %s" % (axis, str(u_max-u_min), str(h), str((u_max-u_min)/h))
		sys.exit()
	totalCells += n

	g.write("    - end: %s\n" % str(u_max))
	g.write("      cells: %d\n" % n)
	g.write("      stretchRatio: 1.0\n")

	L = max - u_max
	current_precision = 1
	next_ratio = 2.
	while current_precision <= precision:
		ratio = next_ratio 
		n = int(round(np.log(L*(ratio-1)/h + 1.)/np.log(ratio)))
		AR = ratio**(n-1)
		if AR < maxAR:
			next_ratio += (0.1)**current_precision
			current_precision+=1
		else:
			next_ratio -= (0.1)**current_precision
	totalCells += n
	
	g.write("    - end: %s\n" % str(max))
	g.write("      cells: %d\n" % n)
	g.write("      stretchRatio: %s\n\n" % str(ratio))

	return totalCells

def generate_grid(inFile, outFile):
	# read the input file
	f = open(inFile, 'r')
	for line in f:
		b = filter(None, re.split('\[|\]|\n|:|,| ', line))
		if b != []:
			if b[0] == 'FinestMeshSpacing':
				h = float(b[1])
			elif b[0] == 'MaxAspectRatio':
				maxAR = float(b[1])
			elif b[0] == 'XRange':
				xmin = float(b[1])
				xmax = float(b[2])
			elif b[0] == 'YRange':
				ymin = float(b[1])
				ymax = float(b[2])
			elif b[0] == 'ZRange':
				zmin = float(b[1])
				zmax = float(b[2])
			elif b[0] == 'UniformRegionXRange':
				u_xmin = float(b[1])
				u_xmax = float(b[2])
			elif b[0] == 'UniformRegionYRange':
				u_ymin = float(b[1])
				u_ymax = float(b[2])
			elif b[0] == 'UniformRegionZRange':
				u_zmin = float(b[1])
				u_zmax = float(b[2])
	f.close()

	# write the output file
	g = open(outFile, 'w')
	nx = calculate_ratios("x", xmin, xmax, u_xmin, u_xmax, h, maxAR, g)
	ny = calculate_ratios("y", ymin, ymax, u_ymin, u_ymax, h, maxAR, g)
	try:
		nz = calculate_ratios("z", zmin, zmax, u_zmin, u_zmax, h, maxAR, g)
	except:
		pass
	g.close()
	print "Grid written to file " + outFile
	print "Size of mesh: %d x %d" % (nx, ny),
	try:
		print "x %d" % nz
	except:
		print " "

def main():
	args = read_inputs()
	generate_grid(args.input, args.output)

if __name__ == "__main__":
	main()