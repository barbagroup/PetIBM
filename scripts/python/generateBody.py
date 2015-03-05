#!/usr/bin/env python

# file: generateBody.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Generates a body file.


import os
import sys
import argparse
import math

import numpy
from matplotlib import pyplot


def read_inputs():
  """Parses the command-line."""
  # create parser
  parser = argparse.ArgumentParser(description='Geometry discretization',
                        formatter_class= argparse.ArgumentDefaultsHelpFormatter)
  # fill parser with arguments
  # geometry arguments
  parser.add_argument('--file', '-f', dest='file_path', type=str,
                      help='path of the coordinates file')
  parser.add_argument('--circle', dest='circle', type=float, nargs='+',
                      help='radius and center-coordinates of the circle')
  parser.add_argument('--line', '-l', dest='line', type=float, nargs='+',
                      help='length and starting-point of the line')
  parser.add_argument('--sphere', dest='sphere', type=float, nargs='+',
                      help='radius and center-coordinates of the sphere')
  # discretization arguments
  parser.add_argument('--n', '-n', dest='n', type=int, default=None,
                      help='target-number of discrete points')
  parser.add_argument('--ds', '-ds', dest='ds', type=float, default=None,
                      help='target-distance between two consecutive points')
  # geometry modification arguments
  parser.add_argument('--rotation', '-r', dest='rotation', type=float, 
                      nargs='+', default=[0.0, 0.0, 0.0],
                      help='rotation (angle in degrees, x-center, y-center)')
  parser.add_argument('--translation', '-t', dest='translation', type=float,
                      nargs='+', default=[0.0, 0.0],
                      help='displacement in the x- and y- directions')
  parser.add_argument('--scale', '-s', dest='scale', type=float, default=1.0,
                      help='scaling factor for 2D geometry')
  parser.add_argument('--extrusion', '-e', dest='extrusion', type=float, 
                      nargs='+',
                      help='z-limits of the 3D cylinder')
  # output arguments
  parser.add_argument('--save-name', dest='save_name', type=str, 
                      default='new_body', help='name of the new body file')
  parser.add_argument('--extension', dest='extension', type=str, default='bdy',
                      help='extension of the output file')
  parser.add_argument('--save-dir', dest='save_directory', type=str, 
                      default=os.getcwd(),
                      help='directory where body file will be saved')
  parser.add_argument('--show', dest='show', action='store_true',
                      help='displays the initial and current geometries')
  return parser.parse_args()


class Geometry:
  """Definition of the geometry."""
  def __init__(self):
    pass

  def write(self, file_path, dimensions=2, cylinder=True):
    """Write the coordinates of the geometry into a file.
    
    Arguments
    ---------
    file_path -- path of the output file
    dimensions -- number of dimensions (default 2)
    cylinder -- is body a cylinder (default True)
    """
    file_name = os.path.splitext(os.path.basename(file_path))[0]
    with open(file_path, 'w') as outfile:
      if dimensions == 2:
        outfile.write('%d\n' % self.x.size)
        numpy.savetxt(outfile, numpy.c_[self.x, self.y],
                      fmt='%.6f', delimiter='\t')
      elif dimensions == 3:
        if cylinder:
          outfile.write('%d\n' % (self.x.size*self.z.size))
          for z in self.z:
            numpy.savetxt(outfile, 
                          numpy.c_[self.x, self.y, z*numpy.ones(self.x.size)],
                          fmt='%.6f', delimiter='\t')
        else:
          outfile.write('%d\n' % self.x.size)
          numpy.savetxt(outfile, numpy.c_[self.x, self.y, self.z], 
                        fmt='%.6f', delimiter='\t')

  def get_distance(self, x_start, y_start, x_end, y_end):
    """Returns the distance between two points.
    
    Arguments
    ---------
    x_start, y_end -- coordinates of the first point
    x_end, y_end -- coordinates of the second point
    """
    return math.sqrt((x_end-x_start)**2 + (y_end-y_start)**2)

  def get_perimeter(self):
    """Returns the perimeter of the geometry."""
    x, y = numpy.append(self.x, self.x[0]), numpy.append(self.y, self.y[0])
    return numpy.sum(numpy.sqrt((x[1:]-x[:-1])**2+(y[1:]-y[:-1])**2))

  def get_mass_center(self):
    """Computes the center of mass."""
    self.x_cm = numpy.mean(self.x)
    self.y_cm = numpy.mean(self.y)

  def rotation(self, parameters=[0.0, None, None]):
    """Rotates the geometry.
    
    Arguments
    ---------
    parameters -- angle (degrees) and center of rotation (default: 0.0, None, None)
    """
    angle = -parameters[0]*math.pi/180.
    if len(parameters) == 1 or not(parameters[1:]):
      self.get_mass_center()
      x_rot, y_rot = self.x_cm, self.y_cm
    else:
      x_rot, y_rot = parameters[1:]
    x_tmp = x_rot + (self.x-x_rot)*math.cos(angle) - (self.y-y_rot)*math.sin(angle)
    y_tmp = y_rot + (self.x-x_rot)*math.sin(angle) + (self.y-y_rot)*math.cos(angle)
    self.x, self.y = x_tmp, y_tmp
    self.get_mass_center()

  def translation(self, displacements):
    """Translates the geometry.
    
    Arguments
    ---------
    displacement -- x- and y- displacements
    """
    self.get_mass_center()
    self.x += displacements[0]
    self.y += displacements[1]
    self.get_mass_center()

  def scale(self, ratio):
    """Scales the geometry.
    
    Arguments
    ---------
    ratio -- scaling ratio
    """
    self.get_mass_center()
    self.x = self.x_cm + ratio*(self.x - self.x_cm)
    self.y = self.y_cm + ratio*(self.y - self.y_cm)

  def extrusion(self, z_limits, ds=None):
    """Creates the third direction.

    Arguments
    ---------
    z_limits -- z-limits of the cylinder
    """
    z_start, z_end = z_limits[0], z_limits[1]
    if not self.ds:
      print 'no ds'
      self.ds = self.get_perimeter()/self.x.size
    n = int(math.ceil((z_end-z_start)/self.ds))
    self.z = numpy.linspace(z_start+0.5*self.ds, z_end-0.5*self.ds, n)

  def discretization(self, n=None, ds=None):
    """Discretizes the geometry 
    given a characteristic length or a nmumber of points
    
    Arguments
    ---------
    n -- number of points (default None)
    ds -- segment-length (default None)
    """
    # calculate either the number of points or characteristic length
    if n and not ds:
      ds = self.get_perimeter()/n
    elif ds and not n:
      n = int(math.ceil(self.get_perimeter()/ds))
      ds = self.get_perimeter()/n
    elif not (n and ds):
      return # keep original discretization
    if n == self.x.size:
      return # if same discretization

    # copy coordinates and initialize new ones
    x_old = numpy.append(self.x, self.x[0])
    y_old = numpy.append(self.y, self.y[0])
    x_new, y_new = numpy.empty(n, dtype=float), numpy.empty(n, dtype=float)
    # first element
    x_new[0], y_new[0] = x_old[0], y_old[0]

    I = 0
    tol = 1.0E-06    # tolerance for interpolation
    for i in xrange(n-1):
      x_start, y_start = x_new[i], y_new[i]
      x_end, y_end = x_old[I+1], y_old[I+1]
      distance = self.get_distance(x_start, y_start, x_end, y_end)
      if ds-distance <= tol:
        # interpolation method
        x_new[i+1], y_new[i+1] = self.interpolation(x_start, y_start, 
                                                    x_end, y_end, 
                                                    ds)
      else:
        # projection method
        while I < x_old.size-2 and ds-distance > tol:
          I += 1
          x_tmp, y_tmp = x_end, y_end
          x_end, y_end = x_old[I+1], y_old[I+1]
          distance = self.get_distance(x_start, y_start, x_end, y_end)
        x_new[i+1], y_new[i+1] = self.projection(x_start, y_start, 
                                                 x_tmp, y_tmp, 
                                                 x_end, y_end, 
                                                 ds)
    # store the new discretization
    self.x, self.y = x_new.copy(), y_new.copy()

  def interpolation(self, x_start, y_start, x_end, y_end, ds):
    """Computes the coordinates of a point 
    by interpolation between two given points given a distance.
    
    Arguments
    ---------
    x_start, y_start -- coordinates of the starting point
    x_end, y_end -- coordinates of the ending point
    ds -- length between the starting point and the interpolated one

    Returns
    -------
    x_target, y_target -- coordinates of the interpolated point
    """
    length = self.get_distance(x_start, y_start, x_end, y_end)
    x_target = x_start + ds/length*(x_end-x_start)
    y_target = y_start + ds/length*(y_end-y_start)
    return x_target, y_target

  def projection(self, x_start, y_start, x_tmp, y_tmp, x_end, y_end, ds):
    """Computes the coordinates of a point
    by projection onto the segment [(x_tmp, y_tmp), (x_end, y_end)]
    such that the distance between (x_start, y_start) and the new point
    is the given cahracteristic length.

    Arguments
    ---------
    x_start, y_start -- coordinates of the starting point
    x_tmp, y_tmp -- coordinates of the intermediate point
    x_end, y_end -- coordinates of the ending point
    ds -- characteristic length

    Returns
    -------
    x_target, y_target -- coordinates of the projected point
    """
    tol = 1.0E-06
    if abs(y_end-y_tmp) >= tol:
      # solve for y
      # coefficients of the second-order polynomial
      a = (x_end-x_tmp)**2 + (y_end-y_tmp)**2
      b = 2.0*( (x_end-x_tmp)*( y_tmp*(x_start-x_end) + y_end*(x_tmp-x_start) ) 
               - y_start*(y_end-y_tmp)**2 )
      c = (y_start**2-ds**2)*(y_end-y_tmp)**2 \
          + (y_tmp*(x_start-x_end) + y_end*(x_tmp-x_start))**2
      # solve the second-order polynomial: ay^2 + by + c = 0
      y = numpy.roots([a, b, c])
      # test if the point belongs to the segment
      test = (y_tmp <= y[0] <= y_end or y_end <= y[0] <= y_tmp)
      y_target = (y[0] if test else y[1])
      x_target = x_tmp + (x_end-x_tmp)/(y_end-y_tmp)*(y_target-y_tmp)
    else:
      # solve for x
      # coefficients of the second-order polynomial
      a = (x_end-x_tmp)**2 + (y_end-y_tmp)**2
      b = 2.0*( (x_end-x_tmp)*(y_tmp-y_start)*(y_end-y_tmp) 
                - x_start*(x_end-x_tmp)**2 
                - x_tmp*(x_end-x_tmp)**2 )
      c = (x_end-x_tmp)**2*((y_tmp-y_start)**2+x_start**2-ds**2) \
          + x_tmp**2*(y_end-y_tmp)**2 \
          - 2*x_tmp*(x_end-x_tmp)*(y_tmp-y_start)*(y_end-y_tmp)
      # solve the second-order polynomial: ax^2 + bx + c = 0
      x = numpy.roots([a, b, c])
      # test if the point belongs to the segment
      test = (x_tmp <= x[0] <= x_end or x_end <= x[0] <= x_tmp)
      x_target = (x[0] if test else x[1])
      y_target = y_tmp + (y_end-y_tmp)/(x_end-x_tmp)*(x_target-x_tmp)
    return x_target, y_target

  def plot(self):
    """Plots the geometry."""
    pyplot.figure()
    pyplot.grid(True)
    pyplot.xlabel(r'$x$', fontsize=18)
    pyplot.ylabel(r'$y$', fontsize=18)
    pyplot.plot(numpy.append(self.x_old, self.x_old[0]), 
                numpy.append(self.y_old, self.y_old[0]),
                label='initial', 
                color='k', ls='-', lw=2, marker='o', markersize=6)
    pyplot.plot(numpy.append(self.x, self.x[0]), 
                numpy.append(self.y, self.y[0]),
                label='current', 
                color='r', ls='-', lw=2, marker='o', markersize=6)
    pyplot.axis('equal')
    pyplot.legend(loc='best', prop={'size': 16})
    pyplot.show()


class Points(Geometry):
  """Contains info about a body read from file."""
  def __init__(self, file_path):
    """Reads and stores the cross-section.

    Arguments
    ---------
    file_path -- path of the file with coordinates
    """
    self.read(file_path)

  def read(self, file_path):
    """Reads and stores coordinates.

    Arguments
    ---------
    file_path -- path of the file with 2D coordinates
    """
    with open(file_path, 'r') as infile:
      self.x, self.y = numpy.loadtxt(infile, dtype=float, delimiter='\t', 
                                     unpack=True, skiprows=1)
    self.x_old, self.y_old = self.x.copy(), self.y.copy()


class Circle(Geometry):
  """Contains info about a circle."""
  def __init__(self, parameters=[0.5, 0.0, 0.0], n=100, ds=None):
    """Creates the circular cross-section.

    Arguments
    ---------
    parameters -- radius and center of the circle (default 0.5, 0.0, 0.0)
    n -- number of segments on the circle (default 100)
    ds -- target segment-length (default None)
    """
    self.radius = parameters[0]
    self.xc, self.yc = parameters[1:]
    self.n = n
    self.ds = ds
    self.create()

  def create(self):
    """Creates coordinates of the circle."""
    if self.ds and not self.n:
      self.n = int(math.ceil(2.*math.pi*self.radius/self.ds))
    elif not (self.ds or self.n):
      print 'Error: missing number of segments or length of segment'
      sys.exit(0)
    theta = numpy.linspace(0.0, 2.0*math.pi, self.n+1)[:-1]
    self.x, self.y = self.radius*numpy.cos(theta), self.radius*numpy.sin(theta)


class Line(Geometry):
  """Contains info about a line."""
  def __init__(self, parameters=[1.0, 0.0, 0.0], n=100, ds=None):
    """Create the line.

    Arguments
    ---------
    parameters -- length and starting-point of the line (default 1.0, 0.0, 0.0)
    n -- number of segments on the line (default 100)
    ds -- target segment-legnth (default None)
    """
    self.length = parameters[0]
    self.x_start, self.y_start = parameters[1:]
    self.n = n
    self.ds = ds
    self.create()

  def create(self):
    """Creates coordinates of the line."""
    if self.ds and not self.n:
      self.n = int(math.ceil(self.length/self.ds))
    elif self.n and not self.ds:
      self.ds = self.length/self.n
    elif not (self.ds or self.n):
      print 'Error: missing number of segments or length of segment'
      sys.exit(0)
    self.x = numpy.linspace(self.x_start, self.x_start+self.length, self.n+1)
    self.y = numpy.zeros_like(self.x)
    self.x_old, self.y_old = self.x.copy(), self.y.copy()


class Sphere(Geometry):
  """Contains info about a spherical body."""
  def __init__(self, parameters=[0.5, 0.0, 0.0, 0.0], n=100, ds=None):
    """Creates the sperical body.

    Arguments
    ---------
    parameters -- radius and center of the sphere (default 0.5, 0.0, 0.0, 0.0)
    n -- number of segments on the great-circle (default 100)
    ds -- target segment-length (default None)
    """
    self.radius = parameters[0]
    self.xc, self.yc, self.zc = parameters[1:]
    self.n = n
    self.ds = ds
    self.create()

  def create(self):
    """Creates coordinates of the sphere."""
    if self.n and not self.ds:
      self.ds = 2.0*math.pi*self.radius/self.n
    elif not (self.ds or self.n):
      print 'Error: missing number of segments or length of segment'
      sys.exit(0)
    n_phi = int(math.ceil(math.pi*self.radius/self.ds))+1
    phi = numpy.linspace(0.0, math.pi, n_phi)[1:-1]
    # north pole
    self.x, self.y, self.z = self.xc, self.yc, self.radius+self.zc
    for phi in phi:
      n_theta = int(math.ceil(2.0*math.pi*self.radius*math.sin(phi)/self.ds))+1
      theta = numpy.linspace(0.0, 2.0*math.pi, n_theta)[:-1]
      self.x = numpy.append(self.x, 
                            self.xc + self.radius*math.sin(phi)*numpy.cos(theta))
      self.y = numpy.append(self.y, 
                            self.yc + self.radius*math.sin(phi)*numpy.sin(theta))
      self.z = numpy.append(self.z, 
                            self.zc + self.radius*math.cos(phi)*numpy.ones(theta.size))
    # south pole
    self.x = numpy.append(self.x, self.xc)
    self.y = numpy.append(self.y, self.yc)
    self.z = numpy.append(self.z, -self.radius+self.zc)


def main():
  """Generates a file containing the coordinates of a body."""
  # parse the command-line
  args = read_inputs()

  output_path = '%s/%s.%s' % (args.save_directory, args.save_name, 
                              args.extension)

  dimensions = (3 if args.extrusion else 2)

  # generate the geometry
  if args.file_path:
    print 'Reading points from file ...'
    body = Points(args.file_path)
    body.scale(args.scale)
    body.rotation(args.rotation)
    body.translation(args.translation)
    body.discretization(n=args.n, ds=args.ds)
    if dimensions == 3:
      print 'Creating cylinder ...'
      body.extrusion(args.extrusion, ds=args.ds)
    if args.show:
      body.plot()
    body.write(output_path, dimensions=dimensions)
  elif args.circle:
    print 'Creating circle ...'
    body = Circle(args.circle, args.n, args.ds)
    if dimensions == 3:
      print 'Creating cylinder ...'
      body.extrusion(args.extrusion, ds=args.ds)
    body.write(output_path, dimensions=dimensions)
  elif args.line:
    print 'Creating line ...'
    body = Line(args.line, args.n, args.ds)
    body.rotation(args.rotation)
    if dimensions == 3:
      print 'Creating flat-plate ...'
      body.extrusion(args.extrusion, ds=args.ds)
    if args.show:
      body.plot()
    body.write(output_path, dimensions=dimensions)
  elif args.sphere:
    print 'Creating sphere ...'
    body = Sphere(args.sphere, args.n, args.ds)
    body.write(output_path, dimensions=3, cylinder=False)

  print 'DONE'


if __name__ == '__main__':
  main()