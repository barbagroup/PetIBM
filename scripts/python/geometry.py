#!/usr/bin/env/ python

# file: geometry.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: To be written.


import os
import math
import copy
import argparse

import numpy
from matplotlib import pyplot
from mayavi import mlab


class Point(object):
  """Contains information about a point."""
  def __init__(self, x, y, z=None):
    """Initializes the position of the point.
    
    Arguments
    ---------
    x, y, z -- coordinates of the point (default: z=0.0)
    """
    self.dimensions = (2 if z == None else 3)
    self.x, self.y, self.z = x, y, z
      
  def as_array(self):
    """Returns the coordinates as a Numpy array."""
    if self.dimensions == 2:
      return numpy.array([self.x, self.y])
    else:
      return numpy.array([self.x, self.y, self.z])

  def distance(self, point=None):
    if not point:
      raise ValueError('a point should be given to compute the distance')
    elif point.dimensions != self.dimensions:
      raise ValueError('the two points should have the same dimension')
    return numpy.linalg.norm(self.as_array() - point.as_array())
  
  def rotation(self, center=None, 
               roll=0.0, yaw=0.0, pitch=0.0, mode='deg'):
    """Performs an intrinsic rotation of the point.
    
    Arguments
    ---------
    center -- center of rotation (default None)
    roll, yaw, pitch -- angles of rotation (default 0.0, 0.0, 0.0)
    mode -- anlges in degrees (deg) or radians (rad) (default deg)
    """
    if center:
      center = center.as_array()
    else:
      return
    if mode == 'deg':
      roll *= math.pi/180.0
      yaw *= math.pi/180.0
      pitch *= math.pi/180.0
    point = self.as_array()
    if self.dimensions == 2:
      Rz = numpy.array([[math.cos(pitch), -math.sin(pitch)],
                        [math.sin(pitch), math.cos(pitch)]])
      self.x, self.y = Rz.dot(point-center) + center
    else:
      Rz = numpy.array([[math.cos(pitch), -math.sin(pitch), 0.0],
                        [math.sin(pitch), math.cos(pitch), 0.0],
                        [0.0, 0.0, 1.0]])
      Ry = numpy.array([[math.cos(yaw), 0.0, math.sin(yaw)],
                        [0.0, 1.0, 0.0],
                        [-math.sin(yaw), 0.0, math.cos(yaw)]])
      Rx = numpy.array([[1.0, 0.0, 0.0],
                        [0.0, math.cos(roll), -math.sin(roll)],
                        [0.0, math.sin(roll), math.cos(roll)]])
      self.x, self.y, self.z = Rz.dot(Ry.dot(Rx.dot(point-center))) + center
        
  def translation(self, displacement=[0.0, 0.0, 0.0]):
    """Translates the point.
    
    Arguments
    ---------
    displacement -- displacement in each direction 
                    (default [0.0, 0.0, 0.0])
    """
    self.x += displacement[0]
    self.y += displacement[1]
    if self.dimensions == 3:
      self.z += displacement[2]


class Geometry(object):
  dimensions = None

  """Contains information about a geometry."""
  def __init__(self, points=None, file_path=None):
    """Initializes the geometry with points.
    
    Arguments
    ---------
    points -- list of points (default None)
    file_path -- path of the file with coordinates (default None)
    """
    if points:
      self.points = points
      self.points_initial = copy.deepcopy(points)
    if file_path:
      self.read_from_file(file_path)
    if self.dimensions or self.points:
      self.get_mass_center()
          
  def read_from_file(self, file_path):
    """Reads the coordinates of the geometry from a file.
    
    Arguments
    ---------
    file_path -- path of the coordinates file
    """
    if not self.dimensions:
      print('\nGet the dimension of the geometry ...')
      with open(file_path, 'r') as infile:
        dim =  len(infile.readlines()[1].strip().split())
        if dim == 2:
          self.__class__ = Geometry2d
        elif dim == 3:
          self.__class__ = Geometry3d
    print('\nRead coordinates from file ...')
    with open(file_path, 'r') as infile:
      coords = numpy.loadtxt(infile, dtype=float, skiprows=1)
    self.points = [Point(*coord) for coord in coords]
    

    self.points_initial = copy.deepcopy(self.points)
          
  def gather_coordinate(self, component, position='current'):
    """Gathers a given component of all points into a Numpy array.
    
    Arguments
    ---------
    component -- component of the point to look for
    position -- position of the body (default current)
    """
    if position == 'current':
      return numpy.array([getattr(point, component) 
                          for point in self.points])
    elif position == 'initial':
      return numpy.array([getattr(point, component) 
                          for point in self.points_initial])
      
  def broadcast_coordinate(self, array, component):
    """Broadcasts Numpy array elemens to point coordinate.
    
    Arguments
    ---------
    array -- Numpy array to broadcast
    component -- point's component to be filled
    """
    for i, point in enumerate(self.points):
      setattr(self.points[i], component, array[i])
      
  def get_mass_center(self):
    """Computes the center of mass of the geometry."""
    x_mass = self.gather_coordinate('x').mean()
    y_mass = self.gather_coordinate('y').mean()
    z_mass = (self.gather_coordinate('z').mean() if self.dimensions == 3 else None)
    self.mass_center = Point(x_mass, y_mass, z_mass)
      
  def translation(self, displacement=[0.0, 0.0, 0.0]):
    """Translates the geometry.
    
    Arguments
    ---------
    displacement -- displacement in each direction 
                    (default [0.0, 0.0, 0.0])
    """
    if not any(displacement):
      return
    print('\nTranslate the geometry ...')
    for i, point in enumerate(self.points):
      self.points[i].translation(displacement)
    self.get_mass_center()
      
  def rotation(self, center=None, 
               roll=0.0, yaw=0.0, pitch=0.0, mode='deg'):
    """Rotates the geometry.
    
    Arguments
    ---------
    center -- center of rotation (default None)
    roll, yaw, pitch -- angles of rotation (default 0.0, 0.0, 0.0)
    mode -- angles in degrees (deg) or radians (rad) (default deg)
    """
    if not any([roll, yaw, pitch]):
      return
    print('\nRotate the geometry ...')
    if not center:
      self.get_mass_center()
    for i, point in enumerate(self.points):
      self.points[i].rotation(self.mass_center, roll, yaw, pitch, mode)
    self.get_mass_center()
      
  def scale(self, ratio=1.0):
    """Scales the geometry.
    
    Arguments
    ---------
    ratio -- scaling ratio (default 1.0)
    """
    if ratio == 1.0:
      return
    print('\nScale the geometry ...')
    components = (['x', 'y', 'z'] if self.dimensions == 3 else ['x', 'y'])
    for component in components:
      center = getattr(self.mass_center, component)
      s = center + ratio*(self.gather_coordinate(component)-center)
      self.broadcast_coordinate(s, component)
  
  def write(self, file_path='{}/new_body'.format(os.getcwd())):
    """Writes the coordinates into a file.
    
    Arguments
    ---------
    file_path -- path of the ouput file (default ./new_body)
    """
    print('\nWrite coordinates into file: {}'.format(file_path))
    x = self.gather_coordinate('x')
    y = self.gather_coordinate('y')
    if self.dimensions == 2:
      coords = numpy.c_[x, y]
    elif self.dimensions == 3:
      z = self.gather_coordinate('z')
      coords = numpy.c_[x, y, z]
    with open(file_path, 'w') as outfile:
      outfile.write('{}\n'.format(x.size))
      numpy.savetxt(outfile, coords, fmt='%.6f', delimiter='\t')


class Geometry2d(Geometry):
  """Contains information about a two-dimensional geometry."""
  dimensions = 2
  
  def __init__(self, points=None, file_path=None):
    """Initializes the geometry.
    
    Arguments
    ---------
    points -- list of points (default None)
    file_path -- path of the file with coordinates (default None)
    """
    Geometry.__init__(self, points, file_path)

  def perimeter(self):
    """Returns the perimeter of the closed geometry."""
    x, y = self.gather_coordinate('x'), self.gather_coordinate('y')
    x, y = numpy.append(x, x[0]), numpy.append(y, y[0])
    return numpy.sum(numpy.sqrt((x[1:] - x[:-1])**2 + (y[1:] - y[:-1])**2))

  def extrusion(self, limits=[-0.5, 0.5], n=None, ds=None):
    """Extrudes the two-dimensional geometry in the z-direction.
    
    Arguments
    ---------
    limits -- limits of the extrusion (default [-0.5, 0.5])
    n -- number of divisions in the z-direction (default None)
    ds -- target segment-legnth (default None)
    """
    print('\nExtrude the geometry in the z-direction...')
    if not (ds or n):
      raise ValueError('both ds and n are set to None')
    elif abs(limits[0]-limits[1]) < 1.0E-06:
      raise ValueError('limits are too close from each other')
    z_start, z_end = limits
    if not n:
      n = int(math.ceil(abs(z_start-z_end)/ds))
    s = math.copysign(1.0, z_end-z_start)
    z = numpy.linspace(z_start+s*0.5*ds, z_end-s*0.5*ds, n+1)
    points = sum(([Point(point.x, point.y, z) for point in self.points] for z in z), [])
    return Geometry3d(points)
  
  def discretization(self, n=None, ds=None):
    """Discretizes the geometry.
  
    Arguments
    ---------
    n -- number of divisions (default None)
    ds -- target segment-length (default None)
    """
    print('\nDiscretize the geometry ...')
    if not (n or ds):
      return
    elif not ds:
      ds = self.perimeter()/n
    else:
      n = int(math.ceil(self.perimeter()/ds))
      ds = self.perimeter()/n

    points_old = copy.deepcopy(self.points) + [self.points[0]]
    self.points = [points_old[0]]

    I = 0
    tolerance = 1.0E-06
    for i in xrange(n-1):
      start = self.points[-1]
      end = points_old[I+1]
      distance = start.distance(end)
      if abs(ds-distance) <= tolerance:
        self.points.append(end)
        I += 1
      elif ds < distance:
        # interpolation method
        self.points.append(self.interpolation(start, end, ds))
      else:
        # projection method
        while I < len(points_old)-2 and ds-distance > tolerance:
          I += 1
          tmp, end = end, points_old[I+1]
          distance = start.distance(end)
          self.points.append(self.projection(start, tmp, end, ds))

  def interpolation(self, start, end, distance):
    """Returns a point interpolated between two given points given the distance.

    Arguments
    ---------
    start, end -- ending-points
    distance -- distance between the starting point and the point to return
    """
    length = start.distance(end)
    start, end = start.as_array(), end.as_array()
    new = start + distance/length*(end-start)
    return Point(*new)

  def projection(self, start, tmp, end, distance):
    """Returns a point projected on the segment [tmp, end] 
    at a given distance from start.

    Arguments
    ---------
    start -- starting point
    tmp -- intermediate point
    end -- end point
    distance -- the distance
    """
    tolerance = 1.0E-06
    if abs(end.y- tmp.y) >= tolerance:
      # solve for y-component
      # coefficients of the second-order polynomial
      a = end.distance(tmp)**2
      b = 2.0*((end.x-tmp.x)*(tmp.y*(start.x-end.x) + end.y*(tmp.x-start.x))
               - start.y*(end.y-tmp.y)**2)
      c = ((start.y**2-distance**2)*(end.y-tmp.y)**2
           + (tmp.y*(start.x-end.x) + end.y*(tmp.x-start.x))**2)
      # solve second-order polynomial: ay^2 + by + c = 0
      y = numpy.roots([a, b, c])
      # test if the first solution belongs to the segment
      test = (tmp.y <= y[0] <= end.y or end.y <= y[0] <= tmp.y)
      y_target = (y[0] if test else y[1])
      x_target = tmp.x + (end.x-tmp.x)/(end.y-tmp.y)*(y_target-tmp.y)
    else:
      # solve for x-component
      # coefficients of the second-order polynomial
      a = end.distance(tmp)**2
      b = 2.0*((end.x-tmp.x)*(tmp.y-start.y)*(end.y-tmp.y) 
               - start.x*(end.x-tmp.x)**2 - tmp.x*(end.x-tmp.x)**2)
      c = ((end.x-tmp.x)**2*((tmp.y-start.y)**2+start.x**2-distance**2)
           + tmp.x**2*(end.y-tmp.y)**2
           - 2.0*tmp.x*(end.x-tmp.x)*(tmp.y-start.y)*(end.y-tmp.y))
      # solve second-order polynomial: ay^2 + by + c = 0
      x = numpy.roots([a, b, c])
      # test if the first solution belongs to the segment
      test = (tmp.x <= x[0] <= end.x or end.x <= x[0] <= tmp.x)
      x_target = (x[0] if test else x[1])
      y_target = tmp.y + (end.y-tmp.y)/(end.x-tmp.x)*(x_target-tmp.x)
    return Point(x_target, y_target)

  def plot(self):
    """Plots the two-dimensional geometry using Matplotlib."""
    print('\nPlot the two-dimensional geometry ...')
    pyplot.style.use('{}/scripts/python/style/'
                     'style_PetIBM.mplstyle'.format(os.environ['PETIBM_DIR']))
    pyplot.grid(True)
    pyplot.xlabel('x')
    pyplot.ylabel('y')
    x = self.gather_coordinate('x')
    y = self.gather_coordinate('y')
    x_init = self.gather_coordinate('x', position='initial')
    y_init = self.gather_coordinate('y', position='initial')
    if len(self.points) == len(self.points_initial):
      same = (numpy.allclose(x, x_init, rtol=1.0E-06) and
              numpy.allclose(y, y_init, rtol=1.0E-06))
    else:
      same = False
    if not same:
      pyplot.plot(x_init, y_init, label='initial', 
                  lw=0, marker='o')
    pyplot.plot(x, y, label='current', lw=0, marker='o')
    pyplot.legend()
    pyplot.axis('equal')
    pyplot.show()


class Line(Geometry2d):
  """Contains information about a line."""
  def __init__(self, start=Point(0.0, 0.0), length=1.0, 
               n=None, ds=None):
    """Initializes the line.
    
    Arguments
    ---------
    start -- starting point of the line (default Point(0.0, 0.0))
    length -- length of the line (default 1.0)
    n -- number of divisions on the line (default None)
    ds -- target segment length (default None)
    """
    self.start = start
    self.length = length
    self.n, self.ds = n, ds
    self.create()
    self.get_mass_center()
  
  def create(self):
    """Creates the points on the line."""
    print('\nCreate a line ...')
    if not (self.n or self.ds):
      raise ValueError('both ds and n are set to None')
    elif not self.n:
      self.n = int(math.ceil(self.length/self.ds))
    x = self.start.x + numpy.linspace(0.0, self.length, self.n+1)
    y = self.start.y * numpy.ones(self.n+1)
    self.points = [Point(x[i], y[i]) for i in xrange(self.n+1)]
    self.points_initial = copy.deepcopy(self.points)


class Circle(Geometry2d):
  """Contains information about a circular geometry."""
  def __init__(self, center=Point(0.0, 0.0), radius=0.5, 
               n=None, ds=None):
    """Initializes the circle.
    
    Arguments
    ---------
    center -- center of the circle (default Point(0.0, 0.0))
    radius -- radius of the circle (default 0.5)
    n -- number of divisions on the circle (default None)
    ds -- target segment-length (default None)
    """
    self.center = center
    self.radius = radius
    self.n, self.ds = n, ds
    self.create()
    self.get_mass_center()
      
  def create(self):
    """Creates the points on the circle."""
    print('\nCreate a circle ...')
    if not (self.n or self.ds):
      raise ValueError('both ds and n are set to None')
    elif not self.n:
      self.n = int(math.ceil(2.0*math.pi*self.radius/self.ds))
    theta = numpy.linspace(0.0, 2.0*math.pi, self.n+1)
    x = self.radius*numpy.cos(theta)[:-1]
    y = self.radius*numpy.sin(theta)[:-1]
    self.points = [Point(x[i], y[i]) for i in xrange(x.size)]
    self.points_initial = copy.deepcopy(self.points)


class Rectangle(Geometry2d):
  """Contains information about a rectangle."""
  def __init__(self, 
               bottom_left=Point(0.0, 0.0), top_right=Point(1.0, 1.0),
               nx=None, ny=None, ds=None):
    """Initializes the rectangle.
    
    Arguments
    ---------
    bottom_left -- bottom-left point (default Point(0.0, 0.0))
    top_right -- top-right point (default Point(1.0, 1.0))
    nx -- number of divisions in the x-direction (default None)
    ny -- number of divisions in the y-direction (default None)
    ds -- target segment-length (default None)
    """
    self.bottom_left, self.top_right = bottom_left, top_right
    self.nx, self.ny, self.ds = nx, ny, ds
    self.create()
    self.get_mass_center()
      
  def create(self):
    """creates the points on the rectangle."""
    print('\nCreate a rectangle ...')
    if self.ds:
      self.nx = int(math.ceil((self.top_right.x-self.bottom_left.x)/self.ds))
      self.ny = int(math.ceil((self.top_right.y-self.bottom_left.y)/self.ds))
    elif not (self.ds or (self.nx and self.ny)):
      raise ValueError('ds is set to None '
                       'while nx and/or ny are set to None')
    self.points = []
    # bottom
    x_bottom = numpy.linspace(self.bottom_left.x, self.top_right.x, self.nx+1)[:-1]
    self.points += [Point(x, self.bottom_left.y) for x in x_bottom]
    # right
    y_right = numpy.linspace(self.bottom_left.y, self.top_right.y, self.ny+1)[:-1]
    self.points += [Point(self.top_right.x, y) for y in y_right]
    # top
    x_top = numpy.linspace(self.top_right.x, self.bottom_left.x, self.nx+1)[:-1]
    self.points += [Point(x, self.top_right.y) for x in x_top]
    # left
    y_left = numpy.linspace(self.top_right.y, self.bottom_left.y, self.ny+1)[:-1]
    self.points += [Point(self.bottom_left.x, y) for y in y_left]
    self.points_initial = copy.deepcopy(self.points)


class Geometry3d(Geometry):
  """Contains information about a three-dimensional geometry."""
  dimensions = 3
  
  def __init__(self, points=None, file_path=None):
    """Initializes the geometry with points.
    
    Arguments
    ---------
    points -- list of points (default None)
    file_path -- path of the file with coordinates (default None)
    """
    Geometry.__init__(self, points, file_path)
  
  def plot(self):
    """Plots the geometry using the package Mayavi."""
    print('\nPlot the three-dimensional geometry ...')
    x_init = self.gather_coordinate('x', position='initial')
    y_init = self.gather_coordinate('y', position='initial')
    z_init = self.gather_coordinate('z', position='initial')
    x = self.gather_coordinate('x')
    y = self.gather_coordinate('y')
    z = self.gather_coordinate('z')
    figure = mlab.figure('body', size=(600, 600))
    figure.scene.disable_render = False
    same = (numpy.allclose(x, x_init, rtol=1.0E-06) and
            numpy.allclose(y, y_init, rtol=1.0E-06) and
            numpy.allclose(z, z_init, rtol=1.0E-06))
    if not same:
      mlab.points3d(x_init, y_init, z_init, name='initial',
                    scale_factor=0.01, color=(0, 0, 1))
    mlab.points3d(x, y, z, name='current',
                  scale_factor=0.01, color=(1, 0, 0))
    mlab.axes()
    mlab.orientation_axes()
    mlab.outline()
    figure.scene.disable_render = True
    mlab.show()


class Sphere(Geometry3d):
  """Contains information about a spherical geometry."""
  def __init__(self, center=Point(0.0, 0.0, 0.0), radius=0.5, 
               n=None, ds=None):
    """Initializes the sphere.
    
    Arguments
    ---------
    center -- center of the sphere (default Point(0.0, 0.0, 0.0))
    radius -- radius of the sphere (default 0.5)
    n -- number of divisions ont eh great-circle (default 10)
    ds -- target segment length (default None)
    """
    self.center = center
    self.radius = radius
    self.n, self.ds = n, ds
    self.create()
    self.get_mass_center()
      
  def create(self):
    """Creates the points on the sphere."""
    print('\nCreate a sphere ...')
    if not (self.n or self.ds):
      raise ValueError('both ds and n are set to None')
    elif self.n and not self.ds:
      self.ds = 2.0*math.pi*self.radius/self.n
    n_phi = int(math.ceil(math.pi*self.radius/self.ds))
    phi = numpy.linspace(0.0, math.pi, n_phi)[1:-1]
    # north pole
    x = self.center.x
    y = self.center.y
    z = self.center.z + self.radius
    # between poles
    for phi in phi:
      rsinphi = self.radius*math.sin(phi)
      rcosphi = self.radius*math.cos(phi)
      n_theta = int(math.ceil(2.0*math.pi*rsinphi/self.ds))
      theta = numpy.linspace(0.0, 2.0*math.pi, n_theta+1)[:-1]
      x = numpy.append(x, self.center.x + rsinphi*numpy.cos(theta))
      y = numpy.append(y, self.center.y + rsinphi*numpy.sin(theta))
      z = numpy.append(z, self.center.z + rcosphi*numpy.ones(n_theta))
    # south pole
    x = numpy.append(x, self.center.x)
    y = numpy.append(y, self.center.y)
    z = numpy.append(z, self.center.z - self.radius)
    # create points
    self.points = [Point(x[i], y[i], z[i]) for i in xrange(x.size)]
    self.points_initial = copy.deepcopy(self.points) 


def main():
  pass


if __name__ == '__main__':
  main()