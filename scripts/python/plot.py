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

	parser = argparse.ArgumentParser(description="Converts the PETSc output to VTK format", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-folder", dest="folder", help="Case folder", default="cases/2d/cylinderRe40")
	parser.add_argument("-xmin", type=float, dest="xmin", help="lower x-limit of the plotting region", default=float("-inf"))
	parser.add_argument("-xmax", type=float, dest="xmax", help="upper x-limit of the plotting region", default=float("inf"))
	parser.add_argument("-ymin", type=float, dest="ymin", help="lower y-limit of the plotting region", default=float("-inf"))
	parser.add_argument("-ymax", type=float, dest="ymax", help="upper y-limit of the plotting region", default=float("inf"))
	parser.add_argument("-vortlim", type=float, dest="vortlim", help="largest absolute value among contour lines displayed", default=3)
	parser.add_argument("-numlevels", type=float, dest="numlevels", help="number of vortex contour line levels (choose an even number)", default=14)
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
	y = grid[nx+1:]

	xmin = x[0] if CLargs.xmin < x[0] else CLargs.xmin
	xmax = x[-1] if CLargs.xmax > x[-1] else CLargs.xmax
	ymin = y[0] if CLargs.ymin < y[0] else CLargs.ymin
	ymax = y[-1] if CLargs.ymax > y[-1] else CLargs.ymax

	dx = np.zeros(nx)
	dy = np.zeros(ny)
	dx[0:nx] = x[1:nx+1]-x[0:nx]
	dy[0:ny] = y[1:ny+1]-y[0:ny]

	xu = np.zeros(Unx)
	yu = np.zeros(Uny)
	xu[0:Unx] = x[1:Unx+1]
	yu[0:Uny] = 0.5*(y[0:Uny]+y[1:Uny+1])

	xv = np.zeros(Vnx)
	yv = np.zeros(Vny)
	xv[0:Vnx] = 0.5*(x[0:Vnx]+x[1:Vnx+1])
	yv[0:Vny] = y[1:Vny+1]

	xo = np.zeros(nx-1)
	yo = np.zeros(ny-1)
	xo[0:nx-1] = x[1:nx]
	yo[0:ny-1] = y[1:ny]

	Ostartx, Ostarty = 0, 0
	Oendx, Oendy = nx-1, ny-1
	
	for i in xrange(nx-1):
		if xmin < xo[i]:
			break
		Ostartx += 1

	for i in reversed(xrange(nx-1)):
		if xmax > xo[i]:
			break
		Oendx -= 1

	for j in xrange(ny-1):
		if ymin < yo[j]:
			break
		Ostarty += 1

	for j in reversed(xrange(ny-1)):
		if ymax > yo[j]:
			break
		Oendy -= 1

	Omg = np.zeros((Oendy-Ostarty, Oendx-Ostartx))

	xp = np.zeros(nx)
	yp = np.zeros(ny)
	xp[0:nx] = 0.5*(x[0:nx]+x[1:nx+1])
	yp[0:ny] = 0.5*(y[0:ny]+y[1:ny+1])

	Pstartx, Pstarty = 0, 0
	Pendx, Pendy = nx, ny
	
	for i in xrange(nx):
		if xmin < xp[i]:
			break
		Pstartx += 1

	for i in reversed(xrange(nx)):
		if xmax > xp[i]:
			break
		Pendx -= 1

	for j in xrange(ny):
		if ymin < yp[j]:
			break
		Pstarty += 1

	for j in reversed(xrange(ny)):
		if ymax > yp[j]:
			break
		Pendy -= 1

	mkdir(folder+"/output")
	plt.ioff()
	
	for n in xrange(args.nsave, args.nt+args.nsave, args.nsave):
		
		# U
		petscObjs = PetscBinaryIO.PetscBinaryIO().readVec('%s/%07d/qx.dat' % (folder,n))[1:]
		U = petscObjs.reshape((Uny, Unx))
		for j in xrange(Uny):
			U[j,:] = U[j,:]/dy[j]
		
		X, Y = np.meshgrid(xu,yu)
		CS = plt.contour(X, Y, U, levels=np.linspace(-0.5, 1, 16))
		plt.colorbar(CS)
		plt.axis([xmin, xmax, ymin, ymax])
		plt.gca().set_aspect('equal', adjustable='box')
		plt.savefig('%s/output/U%07d.png' % (folder,n))
		plt.clf()

		# V
		petscObjs = PetscBinaryIO.PetscBinaryIO().readVec('%s/%07d/qy.dat' % (folder,n))[1:]
		V = petscObjs.reshape((Vny, Vnx))
		for i in xrange(Vnx):
			V[:,i] = V[:,i]/dx[i]
				
		X, Y = np.meshgrid(xv,yv)
		CS = plt.contour(X, Y, V, levels=np.linspace(-0.5, 0.5, 11))
		plt.colorbar(CS)
		plt.axis([xmin, xmax, ymin, ymax])
		plt.gca().set_aspect('equal', adjustable='box')
		plt.savefig('%s/output/V%07d.png' %(folder,n))
		plt.clf()

		# Vorticity
		X, Y = np.meshgrid(xo[Ostartx:Oendx],yo[Ostarty:Oendy])
		for j in xrange(Ostarty,Oendy):
			for i in xrange(Ostartx,Oendx):
				Omg[j-Ostarty, i-Ostartx] = (V[j, i+1]-V[j, i])/(0.5*(dx[i]+dx[i+1])) - (U[j+1, i]-U[j, i])/(0.5*(dy[j]+dy[j+1]))
		
		CS = plt.contour(X, Y, Omg, levels=np.linspace(-CLargs.vortlim, CLargs.vortlim, CLargs.numlevels))
		plt.colorbar(CS)
		plt.axis([xmin, xmax, ymin, ymax])
		plt.gca().set_aspect('equal', adjustable='box')
		plt.savefig('%s/output/O%07d.png' %(folder,n))
		plt.clf()

		#CS = plt.pcolor(X, Y, Omg, cmap='RdBu', vmin=-CLargs.vlim, vmax=CLargs.vlim)
		#plt.colorbar(CS)
		#plt.axis([xmin, xmax, ymin, ymax])
		#plt.gca().set_aspect('equal', adjustable='box')
		#plt.savefig('%s/output/R%07d.png' %(folder,n))
		#plt.clf()

		# Pressure
		petscObjs = PetscBinaryIO.PetscBinaryIO().readVec('%s/%07d/phi.dat' % (folder,n))[1:]
		P = petscObjs.reshape((ny, nx))

		X, Y = np.meshgrid(xp[Pstartx:Pendx],yp[Pstarty:Pendy])
		CS = plt.contour(X, Y, P[Pstarty:Pendy,Pstartx:Pendx], levels=np.linspace(-1.0, 1.0, 21))
		plt.colorbar(CS)
		plt.axis([xmin, xmax, ymin, ymax])
		plt.gca().set_aspect('equal', adjustable='box')
		plt.savefig('%s/output/P%07d.png' %(folder,n))
		plt.clf()
		
		print 'Plots generated for timestep %d.' % n
	