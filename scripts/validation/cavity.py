#!/usr/bin/env python

import argparse
import numpy as np
import os
import os.path
import sys
import csv
import matplotlib.pyplot as plt
sys.path.append(os.path.join(os.environ['PETSC_DIR'],'bin','pythonscripts'))
import PetscBinaryIO

# Parse command line options
parser = argparse.ArgumentParser(description="Runs the validation case for flow in a lid-driven cavity for a specified Reynolds number", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-Re", dest="Re", help="Reynolds number", default='100')
args = parser.parse_args()
Re = args.Re

if Re=='100':
	uCol           = 1
	vCol           = 7
elif Re=='1000':
	uCol           = 2
	vCol           = 8
elif Re=='3200':
	uCol           = 3
	vCol           = 9
elif Re=='5000':
	uCol           = 4
	vCol           = 10
elif Re=='10000':
	uCol           = 5
	vCol           = 11
else:
	print "Unavailable option for Reynolds number. Choose 100, 1000, 3200, 5000 or 10000."
	sys.exit()

dataFile = 'scripts/validation/data/cavity-GGS82.txt'

columns = []
f=open(dataFile, 'rb')
reader = csv.reader(f, delimiter='\t')
for column in zip(*reader):
	columns.append(column)

folder = 'cases/2d/lidDrivenCavity/Re' + Re

infoFile = folder + "/simulationInfo.txt"
	
f = open(infoFile, "r")
args_list = f.read().split()
f.close()

fileParser = argparse.ArgumentParser()
fileParser.add_argument("-nx", type=int, dest="nx", help="number of cells in x-direction", default=32)
fileParser.add_argument("-ny", type=int, dest="ny", help="number of cells in y-direction", default=32)
fileParser.add_argument("-startStep", type=int, dest="startStep", help="number of time steps", default=0)
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

dx = np.zeros(nx)
dy = np.zeros(ny)
dx[0:nx] = x[1:nx+1]-x[0:nx]
dy[0:ny] = y[1:ny+1]-y[0:ny]

yu = np.zeros(Uny)
yu[0:Uny] = 0.5*(y[0:Uny]+y[1:Uny+1])

xv = np.zeros(Vnx)
xv[0:Vnx] = 0.5*(x[0:Vnx]+x[1:Vnx+1])

# U
petscObjs = PetscBinaryIO.PetscBinaryIO().readVec('%s/%07d/qx.dat' % (folder,args.nt))[1:]
U = petscObjs.reshape((Uny, Unx))
for j in xrange(Uny):
	U[j,:] = U[j,:]/dy[j]

# V
petscObjs = PetscBinaryIO.PetscBinaryIO().readVec('%s/%07d/qy.dat' % (folder,args.nt))[1:]
V = petscObjs.reshape((Vny, Vnx))
for i in xrange(Vnx):
	V[:,i] = V[:,i]/dx[i]

plt.plot(columns[0], columns[uCol], 'o')
plt.plot(yu, U[:,nx/2-1], '-')
plt.axis((0,1,-0.7,1.3))
plt.xlabel('y-coordinate')
plt.ylabel('u-velocity along centerline')
plt.title('Lid-driven cavity at Re=%s' % Re)
plt.savefig("scripts/validation/uRe%s.png" % Re)
plt.clf()

plt.plot(columns[6], columns[vCol], 'o')
plt.plot(xv, V[ny/2-1,:], '-')
plt.axis((0,1,-0.7,1.3))
plt.xlabel('x-coordinate')
plt.ylabel('v-velocity along centerline')
plt.title('Lid-driven cavity at Re=%s' % Re)
plt.savefig("scripts/validation/vRe%s.png" % Re)