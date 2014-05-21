#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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
		folder = "cases/cavityRe100"
	
	print "Case folder: " + folder
	
	Unx, Uny = 31, 32
	Vnx, Vny = 32, 31
	
	print "Size of U-array: %d x %d" % (Unx, Uny)
	print "Size of V-array: %d x %d" % (Vnx, Vny)
	
	dx = 1./32
	dy = 1./32
	
	mkdir(folder+"/output")
	plt.ioff()
	
	nsave = 1000
	nt = 1000
	for n in range(nsave, nt+nsave, nsave):
		
		# U
		petscObjs = PetscBinaryIO.PetscBinaryIO().readVec('%s/%07d/qx.dat' % (folder,n))[1:]
		U = petscObjs.reshape((Uny, Unx))/dy
		
		x  = np.linspace(dx, 1.-dx, Unx)
		y  = np.linspace(dy/2, 1.-dy/2, Uny)
		
		X, Y = np.meshgrid(x,y)
		CS = plt.contour(X, Y, U, levels=np.linspace(-0.5, 1, 16))
		plt.colorbar(CS)
		plt.savefig('%s/output/U%07dcontour.png' % (folder,n))
		plt.clf()
		ax = Axes3D(plt.gcf())
		ax.plot_wireframe(X, Y, U)
		plt.savefig('%s/output/U%07dwireframe.png' % (folder,n))
		plt.clf()
		
		# V
		petscObjs = PetscBinaryIO.PetscBinaryIO().readVec('%s/%07d/qy.dat' % (folder,n))[1:]
		V = petscObjs.reshape((Vny, Vnx))/dx
		
		x  = np.linspace(dx/2, 1.-dx/2, Vnx)
		y  = np.linspace(dy, 1-dy, Vny)
		
		X, Y = np.meshgrid(x,y)
		CS = plt.contour(X, Y, V, levels=np.linspace(-0.5, 1, 16))
		plt.colorbar(CS)
		plt.savefig('%s/output/V%07dcontour.png' % (folder,n))
		plt.clf()
		ax = Axes3D(plt.gcf())
		ax.plot_wireframe(X, Y, V)
		plt.savefig('%s/output/V%07dwireframe.png' % (folder,n))
		plt.clf()
		
		print 'Plots generated for timestep %d.' % n
