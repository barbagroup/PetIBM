#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import sys
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import errno
sys.path.append(os.path.join(os.environ['PETSC_DIR'],'bin','pythonscripts'))
import PetscBinaryIO

def isNumber(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
        
def mkdir(path, overwrite=False):
	try:
		os.makedirs(path)
	except OSError as exception:
		if exception.errno == errno.EEXIST:
			if not overwrite:
				#print "Path '%s' already exists" % path
				return
		else: raise

if __name__=="__main__":

	# command line arguments
	# read the case folder, the extent of the plotting region, the stride, and the time steps to plot
	parser = argparse.ArgumentParser(description="Converts the PETSc output to VTK format", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-folder", dest="folder", help="Case folder", default="cases/2d/cavityRe100")
	parser.add_argument("-xmin", type=float, dest="xmin", help="lower x-limit of the plotting region", default=float("-inf"))
	parser.add_argument("-xmax", type=float, dest="xmax", help="upper x-limit of the plotting region", default=float("inf"))
	parser.add_argument("-ymin", type=float, dest="ymin", help="lower y-limit of the plotting region", default=float("-inf"))
	parser.add_argument("-ymax", type=float, dest="ymax", help="upper y-limit of the plotting region", default=float("inf"))
	parser.add_argument("-stride", type=int, dest="stride", help="output only the nodes separated by a certain interval", default=1)
	parser.add_argument("-startStep", type=int, dest="startStep", help="start step", default=-1)
	parser.add_argument("-nsave", type=int, dest="nsave", help="nsave", default=-1)
	parser.add_argument("-nt", type=int, dest="nt", help="nsave", default=-1)
	CLargs = parser.parse_args()

	folder = CLargs.folder
	stride = CLargs.stride
		
	print "Case folder: " + folder
	
	infoFile = folder + "/simulationInfo.txt"
	f = open(infoFile, "r")
	args_list = f.read().split()
	f.close()
	
	# read information from the file simulationInfo.txt
	fileParser = argparse.ArgumentParser()
	fileParser.add_argument("-nx", type=int, dest="nx", help="number of cells in x-direction", default=32)
	fileParser.add_argument("-ny", type=int, dest="ny", help="number of cells in y-direction", default=32)
	fileParser.add_argument("-startStep", type=int, dest="startStep", help="start step", default=-1)
	fileParser.add_argument("-nt", type=int, dest="nt", help="number of time steps", default=200)
	fileParser.add_argument("-nsave", type=int, dest="nsave", help="data save stride", default=100)
	fileParser.add_argument("-dt", type=float, dest="dt", help="time increment", default=0.01)
	fileParser.add_argument("-xperiodic", dest="xperiodic", help="periodicity in x-direction", default="False")
	fileParser.add_argument("-yperiodic", dest="yperiodic", help="periodicity in y-direction", default="False")
	args = fileParser.parse_args(args_list)
	
	# if the command line arguments are missing, set the values from the file
	startStep = CLargs.startStep if CLargs.startStep > -1 else args.startStep
	nsave = CLargs.nsave if CLargs.nsave > -1 else args.nsave
	nt = CLargs.nt if CLargs.nt > -1 else args.nt

	# number of cells in the mesh in each direction
	nx = args.nx
	ny = args.ny
	xperiodic = True if args.xperiodic.lower()=='true' else False
	yperiodic = True if args.yperiodic.lower()=='true' else False
	
	# number of u, v variables in each direction
	# this depends on whether the domain in a particular direction is periodic
	Unx, Uny = nx if xperiodic else nx-1, ny
	Vnx, Vny = nx, ny if yperiodic else ny-1
	
	print "Size of U-array: %d x %d" % (Unx, Uny)
	print "Size of V-array: %d x %d" % (Vnx, Vny)
	
	# read the coordinates of the nodes of the grid
	grid = np.loadtxt(folder+"/grid.txt")
	x = grid[:nx+1]
	y = grid[nx+1:nx+1+ny+1]
	
	# create arrays to store the values of the cell widths
	dx = np.zeros(nx)
	dy = np.zeros(ny)
	
	# calculate the cell widths
	dx[0:nx] = x[1:nx+1]-x[0:nx]
	dy[0:ny] = y[1:ny+1]-y[0:ny]
	
	# create arrays to store the coordinates of the locations of velocities
	# these are the locations at which the VTK file stores velocities
	# and they correspond to the centers of the cells of the mesh
	X = np.zeros(Unx-1)
	Y = np.zeros(Vny-1)
	X[0:Unx-1] = 0.5*(x[1:Unx]+x[2:Unx+1])
	Y[0:Vny-1] = 0.5*(y[1:Vny]+y[2:Vny+1])

	# indices to store the start and end locations of the plotting region
	startx, starty = 0, 0
	endx, endy = Unx-1, Vny-1
	
	# calculate the above indices using the command line options
	# both the start and the end index points are inside the plotting region
	for i in xrange(Unx-1):
		if CLargs.xmin < X[i]:
			break
		startx += 1

	for i in reversed(xrange(Unx-1)):
		if CLargs.xmax > X[i]:
			break
		endx -= 1

	for j in xrange(Vny-1):
		if CLargs.ymin < Y[j]:
			break
		starty += 1

	for j in reversed(xrange(Vny-1)):
		if CLargs.ymax > Y[j]:
			break
		endy -= 1

	print startx, endx
	print starty, endy
	print stride
	
	mkdir(folder+"/output")
	
	for n in xrange(startStep+nsave, nt+nsave, nsave):
		# read the fluxes from file and reshape them accordingly
		petscObjs = PetscBinaryIO.PetscBinaryIO().readBinaryFile('%s/%07d/qx.dat' % (folder,n))[0]
		qx = petscObjs.reshape((Uny, Unx))
		
		petscObjs = PetscBinaryIO.PetscBinaryIO().readBinaryFile('%s/%07d/qy.dat' % (folder,n))[0]
		qy = petscObjs.reshape((Vny, Vnx))

		xsize = len(xrange(startx+1, endx+1, stride))
		ysize = len(xrange(starty+1, endy+1, stride))

		# write the VTK file
		outFile = '%s/output/velocity%07d.vtk' % (folder,n)
		g = open(outFile, 'w')	
		g.write('# vtk DataFile Version 3.0\n')
		g.write('Header\n')
		g.write('ASCII\n')
		g.write('DATASET RECTILINEAR_GRID\n')
		g.write('DIMENSIONS %d %d 1\n' % (xsize, ysize))
		g.write('X_COORDINATES %d double\n' % xsize)
		for i in xrange(startx+1, endx+1, stride):
			g.write('%f ' % X[i-1])
		g.write('\n')
		g.write('Y_COORDINATES %d double\n' % ysize)
		for j in xrange(starty+1,endy+1, stride):
			g.write('%f ' % Y[j-1])
		g.write('\n')
		g.write('Z_COORDINATES 1 double\n0.0\n')
	
		g.write("POINT_DATA %d\n" % (xsize*ysize))
		g.write('VECTORS velocity double\n')
		for j in xrange(starty+1,endy+1, stride):
			for i in xrange(startx+1, endx+1, stride):
				g.write( "%f\t%f\t0.0\n" % ( 0.5*(qx[j][i-1]+qx[j][i])/dy[j], 0.5*(qy[j-1][i]+qy[j][i])/dx[i]) )
	
		g.close()
		
		print 'Wrote file ' + outFile + '.'
