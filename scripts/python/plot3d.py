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
	parser.add_argument("-folder", dest="folder", help="Case folder", default="cases/3d/cavityX")
	parser.add_argument("-xmin", type=float, dest="xmin", help="lower x-limit of the plotting region", default=float("-inf"))
	parser.add_argument("-xmax", type=float, dest="xmax", help="upper x-limit of the plotting region", default=float("inf"))
	parser.add_argument("-ymin", type=float, dest="ymin", help="lower y-limit of the plotting region", default=float("-inf"))
	parser.add_argument("-ymax", type=float, dest="ymax", help="upper y-limit of the plotting region", default=float("inf"))
	parser.add_argument("-zmin", type=float, dest="zmin", help="lower z-limit of the plotting region", default=float("-inf"))
	parser.add_argument("-zmax", type=float, dest="zmax", help="upper z-limit of the plotting region", default=float("inf"))
	parser.add_argument("-stride", type=int, dest="stride", help="output only the nodes separated by a certain interval", default=1)
	parser.add_argument("-startStep", type=int, dest="startStep", help="start step", default=-1)
	parser.add_argument("-nsave", type=int, dest="nsave", help="nsave", default=-1)
	parser.add_argument("-nt", type=int, dest="nt", help="nt", default=-1)
	CLargs = parser.parse_args()

	folder = CLargs.folder
	stride = CLargs.stride
	
	print "Case folder:", folder
	
	infoFile = folder + "/simulationInfo.txt"
	f = open(infoFile, "r")
	args_list = f.read().split()
	f.close()

	# read information from the file simulationInfo.txt
	fileParser = argparse.ArgumentParser()
	fileParser.add_argument("-nx", type=int, dest="nx", help="number of cells in x-direction", default=32)
	fileParser.add_argument("-ny", type=int, dest="ny", help="number of cells in y-direction", default=32)
	fileParser.add_argument("-nz", type=int, dest="nz", help="number of cells in z-direction", default=32)
	fileParser.add_argument("-startStep", type=int, dest="startStep", help="start step", default=-1)
	fileParser.add_argument("-nt", type=int, dest="nt", help="number of time steps", default=200)
	fileParser.add_argument("-nsave", type=int, dest="nsave", help="data save stride", default=100)
	fileParser.add_argument("-xperiodic", dest="xperiodic", help="periodicity in x-direction", default="False")
	fileParser.add_argument("-yperiodic", dest="yperiodic", help="periodicity in y-direction", default="False")
	fileParser.add_argument("-zperiodic", dest="zperiodic", help="periodicity in z-direction", default="False")
	args = fileParser.parse_args(args_list)

	# if the command line arguments are missing, set the values from the file
	startStep = CLargs.startStep if CLargs.startStep > -1 else args.startStep
	nsave = CLargs.nsave if CLargs.nsave > -1 else args.nsave
	nt = CLargs.nt if CLargs.nt > -1 else args.nt
	
	# number of cells in the mesh in each direction
	nx = args.nx
	ny = args.ny
	nz = args.nz
	xperiodic = True if args.xperiodic.lower()=='true' else False
	yperiodic = True if args.yperiodic.lower()=='true' else False
	zperiodic = True if args.zperiodic.lower()=='true' else False

	# number of u, v, w variables in each direction
	# this depends on whether the domain in a particular direction is periodic
	Unx, Uny, Unz = nx if xperiodic else nx-1, ny, nz
	Vnx, Vny, Vnz = nx, ny if yperiodic else ny-1, nz
	Wnx, Wny, Wnz = nx, ny, nz if zperiodic else nz-1
	
	print "Size of U-array: %d x %d x %d" % (Unx, Uny, Unz)
	print "Size of V-array: %d x %d x %d" % (Vnx, Vny, Vnz)
	print "Size of W-array: %d x %d x %d" % (Wnx, Wny, Wnz)
	
	# read the coordinates of the nodes of the grid
	grid = np.loadtxt(folder+"/grid.txt")
	x = grid[:nx+1]
	y = grid[nx+1:nx+1+ny+1]
	z = grid[nx+1+ny+1:]
	
	# create arrays to store the values of the cell widths
	dx = np.zeros(nx)
	dy = np.zeros(ny)
	dz = np.zeros(nz)
	
	# calculate the cell widths
	dx[0:nx] = x[1:nx+1]-x[0:nx]
	dy[0:ny] = y[1:ny+1]-y[0:ny]
	dz[0:nz] = z[1:nz+1]-z[0:nz]
	
	# create arrays to store the coordinates of the locations of velocities
	# these are the locations at which the VTK file stores velocities
	# and they correspond to the centers of the cells of the mesh
	X = np.zeros(Unx-1)
	Y = np.zeros(Vny-1)
	Z = np.zeros(Wnz-1)
	X[0:Unx-1] = 0.5*(x[1:Unx]+x[2:Unx+1])
	Y[0:Vny-1] = 0.5*(y[1:Vny]+y[2:Vny+1])
	Z[0:Wnz-1] = 0.5*(z[1:Wnz]+z[2:Wnz+1])

	# indices to store the start and end locations of the plotting region
	startx, starty, startz = 0, 0, 0
	endx, endy, endz = Unx-1, Vny-1, Wnz-1
	
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

	for k in xrange(Wnz-1):
		if CLargs.zmin < Z[k]:
			break
		startz += 1

	for k in reversed(xrange(Wnz-1)):
		if CLargs.zmax > Z[k]:
			break
		endz -= 1

	print startx, endx
	print starty, endy
	print startz, endz
	print stride

	mkdir(folder+"/output")

	for n in xrange(startStep+nsave, nt+nsave, nsave):
		# read the fluxes from file and reshape them accordingly
		petscObjs = PetscBinaryIO.PetscBinaryIO().readBinaryFile('%s/%07d/qx.dat' % (folder,n))[0]
		qx = petscObjs.reshape((Unz, Uny, Unx))
		
		petscObjs = PetscBinaryIO.PetscBinaryIO().readBinaryFile('%s/%07d/qy.dat' % (folder,n))[0]
		qy = petscObjs.reshape((Vnz, Vny, Vnx))

		petscObjs = PetscBinaryIO.PetscBinaryIO().readBinaryFile('%s/%07d/qz.dat' % (folder,n))[0]
		qz = petscObjs.reshape((Wnz, Wny, Wnx))

		xsize = len(xrange(startx+1, endx+1, stride))
		ysize = len(xrange(starty+1, endy+1, stride))
		zsize = len(xrange(startz+1, endz+1, stride))

		# write the VTK file
		outFile = '%s/output/velocity%07d.vtk' % (folder,n)
		g = open(outFile, 'w')
		g.write('# vtk DataFile Version 3.0\n')
		g.write('Header\n')
		g.write('ASCII\n')
		g.write('DATASET RECTILINEAR_GRID\n')
		g.write('DIMENSIONS %d %d %d\n' % (xsize, ysize, zsize))
		g.write('X_COORDINATES %d double\n' % xsize)
		for i in xrange(startx+1, endx+1, stride):
			g.write('%f ' % X[i-1])
		g.write('\n')
		g.write('Y_COORDINATES %d double\n' % ysize)
		for j in xrange(starty+1, endy+1, stride):
			g.write('%f ' % Y[j-1])
		g.write('\n')
		g.write('Z_COORDINATES %d double\n' % zsize)
		for k in xrange(startz+1, endz+1, stride):
			g.write('%f ' % Z[k-1])
		g.write('\n')
	
		g.write("POINT_DATA %d\n" % (xsize*ysize*zsize))
		g.write('VECTORS velocity double\n')
		for k in xrange(startz+1,endz+1, stride):
			for j in xrange(starty+1,endy+1, stride):
				for i in xrange(startx+1,endx+1, stride):
					g.write( "%f\t%f\t%f\n" % ( 0.5*(qx[k][j][i-1]+qx[k][j][i])/(dy[j]*dz[k]), 0.5*(qy[k][j-1][i]+qy[k][j][i])/(dx[i]*dz[k]), 0.5*(qz[k-1][j][i]+qz[k][j][i])/(dx[i]*dy[j]) ) )
	
		g.close()
		
		print 'Wrote file ' + outFile + '.'
