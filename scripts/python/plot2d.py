#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import sys
import numpy as np
import matplotlib.pyplot as plt
import array
import h5py
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

	parser = argparse.ArgumentParser(description="Converts the PETSc output to VTK format", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-folder", dest="folder", help="Case folder", default="cases/2d/cavityRe100")
	parser.add_argument("-xmin", type=float, dest="xmin", help="lower x-limit of the plotting region", default=float("-inf"))
	parser.add_argument("-xmax", type=float, dest="xmax", help="upper x-limit of the plotting region", default=float("inf"))
	parser.add_argument("-ymin", type=float, dest="ymin", help="lower y-limit of the plotting region", default=float("-inf"))
	parser.add_argument("-ymax", type=float, dest="ymax", help="upper y-limit of the plotting region", default=float("inf"))
	CLargs = parser.parse_args()

	folder = CLargs.folder
		
	print "Case folder: " + folder
	
	infoFile = folder + "/simulationInfo.txt"
	
	f = open(infoFile, "r")
	args_list = f.read().split()
	f.close()
	
	fileParser = argparse.ArgumentParser()
	fileParser.add_argument("-nx", type=int, dest="nx", help="number of cells in x-direction", default=32)
	fileParser.add_argument("-ny", type=int, dest="ny", help="number of cells in y-direction", default=32)
	fileParser.add_argument("-nt", type=int, dest="nt", help="number of time steps", default=200)
	fileParser.add_argument("-nsave", type=int, dest="nsave", help="data save stride", default=100)
	fileParser.add_argument("-xperiodic", dest="xperiodic", help="periodicity in x-direction", default="False")
	fileParser.add_argument("-yperiodic", dest="yperiodic", help="periodicity in y-direction", default="False")
	args = fileParser.parse_args(args_list)
	
	nx = args.nx
	ny = args.ny
	xperiodic = True if args.xperiodic.lower()=='true' else False
	yperiodic = True if args.yperiodic.lower()=='true' else False
	Unx, Uny = nx if xperiodic else nx-1, ny
	Vnx, Vny = nx, ny if yperiodic else ny-1
	
	print "Size of U-array: %d x %d" % (Unx, Uny)
	print "Size of V-array: %d x %d" % (Vnx, Vny)
	
	grid = np.loadtxt(folder+"/grid.txt")
	x = grid[:nx+1]
	y = grid[nx+1:nx+1+ny+1]
	
	dx = np.zeros(nx)
	dy = np.zeros(ny)
	
	dx[0:nx] = x[1:nx+1]-x[0:nx]
	dy[0:ny] = y[1:ny+1]-y[0:ny]
	
	X = np.zeros(Unx-1)
	Y = np.zeros(Uny-1)
	X[0:Unx-1] = 0.5*(x[1:Unx]+x[2:Unx+1])
	Y[0:Vny-1] = 0.5*(y[1:Vny]+y[2:Vny+1])

	startx, starty = 0, 0
	endx, endy = Unx-1, Vny-1
	
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
	
	mkdir(folder+"/output")
	
	for n in xrange(args.nsave, args.nt+args.nsave, args.nsave):
		petscObjs = PetscBinaryIO.PetscBinaryIO().readVec('%s/%07d/qx.dat' % (folder,n))[1:]
		qx = petscObjs.reshape((Uny, Unx))
		
		petscObjs = PetscBinaryIO.PetscBinaryIO().readVec('%s/%07d/qy.dat' % (folder,n))[1:]
		qy = petscObjs.reshape((Vny, Vnx))

		outFile = '%s/output/velocity%07d.vtk' % (folder,n)
		g = open(outFile, 'w')
	
		g.write('# vtk DataFile Version 3.0\n')
		g.write('Header\n')
		g.write('ASCII\n')

		g.write('DATASET RECTILINEAR_GRID\n')
		g.write('DIMENSIONS %d %d 1\n' % (endx-startx, endy-starty))
		g.write('X_COORDINATES %d double\n' % (endx-startx))
		for i in xrange(Unx-1):
			g.write('%f ' % X[i])
		g.write('\n')
		g.write('Y_COORDINATES %d double\n' % (endy-starty))
		for j in xrange(Vny-1):
			g.write('%f ' % Y[j])
		g.write('\n')
		g.write('Z_COORDINATES 1 double\n0.0\n')
	
		g.write("POINT_DATA %d\n" % ((endx-startx)*(endy-starty)))
		g.write('VECTORS velocity double\n')
		for j in xrange(starty+1,endy+1):
			for i in xrange(startx+1,endx+1):
				g.write( "%f\t%f\t0.0\n" % ( 0.5*(qx[j][i-1]+qx[j][i])/dy[j], 0.5*(qy[j-1][i]+qy[j][i])/dx[i]) )
	
		g.close()
		
		print 'Wrote file ' + outFile + '.'
