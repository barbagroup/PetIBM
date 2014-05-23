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
				print "Path '%s' already exists" % path
				return
		else: raise

if __name__=="__main__":

	try:
		folder = sys.argv[1]
	except IndexError:
		folder = "cases/3d/cavityRe100"
		
	print "Case folder: " + folder
	
	infoFile = folder + "/simulationInfo.txt"
	
	f = open(infoFile, "r")
	args_list = f.read().split()
	f.close()
	
	fileParser = argparse.ArgumentParser()
	fileParser.add_argument("-nx", type=int, dest="nx", help="number of cells in x-direction", default=32)
	fileParser.add_argument("-ny", type=int, dest="ny", help="number of cells in y-direction", default=32)
	fileParser.add_argument("-nz", type=int, dest="nz", help="number of cells in z-direction", default=32)
	fileParser.add_argument("-nt", type=int, dest="nt", help="number of time steps", default=200)
	fileParser.add_argument("-nsave", type=int, dest="nsave", help="data save stride", default=100)
	fileParser.add_argument("-xperiodic", dest="xperiodic", help="periodicity in x-direction", default="False")
	fileParser.add_argument("-yperiodic", dest="yperiodic", help="periodicity in y-direction", default="False")
	fileParser.add_argument("-zperiodic", dest="zperiodic", help="periodicity in z-direction", default="False")
	args = fileParser.parse_args(args_list)
	
	nx = args.nx
	ny = args.ny
	nz = args.nz
	xperiodic = True if args.xperiodic.lower()=='true' else False
	yperiodic = True if args.yperiodic.lower()=='true' else False
	zperiodic = True if args.zperiodic.lower()=='true' else False
	Unx, Uny, Unz = nx if xperiodic else nx-1, ny, nz
	Vnx, Vny, Vnz = nx, ny if yperiodic else ny-1, nz
	Wnx, Wny, Wnz = nx, ny, nz if zperiodic else nz-1
	
	print "Size of U-array: %d x %d x %d" % (Unx, Uny, Unz)
	print "Size of V-array: %d x %d x %d" % (Vnx, Vny, Vnz)
	print "Size of W-array: %d x %d x %d" % (Wnx, Wny, Wnz)
	
	grid = np.loadtxt(folder+"/grid.txt")
	x = grid[:nx+1]
	y = grid[nx+1:nx+1+ny+1]
	z = grid[nx+1+ny+1:]
	
	dx = np.zeros(nx)
	dy = np.zeros(ny)
	dz = np.zeros(nz)
	
	dx[0:nx] = x[1:nx+1]-x[0:nx]
	dy[0:ny] = y[1:ny+1]-y[0:ny]
	dz[0:nz] = z[1:nz+1]-z[0:nz]
	
	X = np.zeros(Unx-1)
	Y = np.zeros(Uny-1)
	Z = np.zeros(Unz-1)
	X[0:Unx-1] = 0.5*(x[1:Unx]+x[2:Unx+1])
	Y[0:Vny-1] = 0.5*(y[1:Vny]+y[2:Vny+1])
	Z[0:Wnz-1] = 0.5*(z[1:Wnz]+z[2:Wnz+1])
	
	mkdir(folder+"/output")
	
	for n in xrange(args.nsave, args.nt+args.nsave, args.nsave):
		petscObjs = PetscBinaryIO.PetscBinaryIO().readVec('%s/%07d/qx.dat' % (folder,n))[1:]
		U = petscObjs.reshape((Unz, Uny, Unx))
		
		petscObjs = PetscBinaryIO.PetscBinaryIO().readVec('%s/%07d/qy.dat' % (folder,n))[1:]
		V = petscObjs.reshape((Vnz, Vny, Vnx))

		petscObjs = PetscBinaryIO.PetscBinaryIO().readVec('%s/%07d/qz.dat' % (folder,n))[1:]
		W = petscObjs.reshape((Wnz, Wny, Wnx))

		outFile = '%s/output/velocity%07d.vtk' % (folder,n)
		g = open(outFile, 'w')
	
		g.write('# vtk DataFile Version 3.0\n')
		g.write('Header\n')
		g.write('ASCII\n')

		g.write('DATASET RECTILINEAR_GRID\n')
		g.write('DIMENSIONS %d %d %d\n' % (Unx-1, Vny-1, Wnz-1))
		g.write('X_COORDINATES %d double\n' % (Unx-1))
		for i in xrange(Unx-1):
			g.write('%f ' % X[i])
		g.write('\n')
		g.write('Y_COORDINATES %d double\n' % (Vny-1))
		for j in xrange(Vny-1):
			g.write('%f ' % Y[j])
		g.write('\n')
		g.write('Z_COORDINATES %d double\n' % (Wnz-1))
		for k in xrange(Wnz-1):
			g.write('%f ' % Z[k])
		g.write('\n')
	
		g.write("POINT_DATA %d\n" % ((Unx-1)*(Vny-1)*(Wnz-1)))
		g.write('VECTORS velocity double\n')
		for k in xrange(1,Wnz):
			for j in xrange(1,Vny):
				for i in xrange(1,Unx):
					g.write( "%f\t%f\t%f\n" % ( 0.5*(U[k][j][i-1]+U[k][j][i])/(dy[j]*dz[k]), 0.5*(V[k][j-1][i]+V[k][j][i])/(dx[i]*dz[k]), 0.5*(W[k-1][j][i]+W[k][j][i])/(dx[i]*dy[j]) ) )
	
		g.close()
		
		print 'Wrote file ' + outFile + '.'
