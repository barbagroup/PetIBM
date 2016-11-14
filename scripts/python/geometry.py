"""
Collection of classes and functions to generate the coordinates of a geometry.
"""

import os
import math
import copy

import numpy
from matplotlib import pyplot


class Point(object):
  """
  Contains information about a point.
  """
  def __init__(self, x, y, z=None):
    """
    Initializes the position of the point.

    Parameters
    ----------
    x, y, z: floats
      Coordinates of the point; z is optional.
    """
    self.dimensions = (2 if not z else 3)
    self.x, self.y, self.z = x, y, z

  def as_array(self):
    """
    Returns the coordinates as a Numpy array.
    """
    if self.dimensions == 2:
      return numpy.array([self.x, self.y])
    else:
      return numpy.array([self.x, self.y, self.z])

  def distance(self, point=None):
    """
    Computes the distance between the point and another one.

    Parameters
    ----------
    point: Point object, optional
      The other point;
      default: None.

    Returns
    -------
    distance: float
      The distance between the two points.
    """
    if not point:
      raise ValueError('a point should be given to compute the distance')
    elif point.dimensions != self.dimensions:
      raise ValueError('the two points should have the same dimension')
    return numpy.linalg.norm(self.as_array() - point.as_array())

  def rotation(self, center=None,
               roll=0.0, yaw=0.0, pitch=0.0, mode='deg'):
    """
    Performs an intrinsic rotation of the point.

    Parameters
    ----------
    center: Point object, optional
      Center of rotation;
      default: None.
    roll, yaw, pitch: floats, optional
      Angles of rotation in radians or degrees;
      default: 0.0, 0.0, 0.0.
    mode: string, optional
      Angles in degrees ('deg') or radians ('rad')?;
      default: 'deg'.
    """
    if center:
      center = center.as_array()
    else:
      return
    if mode == 'deg':
      roll *= math.pi / 180.0
      yaw *= math.pi / 180.0
      pitch *= math.pi / 180.0
    point = self.as_array()
    if self.dimensions == 2:
      Rz = numpy.array([[math.cos(pitch), -math.sin(pitch)],
                        [math.sin(pitch), math.cos(pitch)]])
      self.x, self.y = Rz.dot(point - center) + center
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
      self.x, self.y, self.z = Rz.dot(Ry.dot(Rx.dot(point - center))) + center

  def translation(self, displacement=[0.0, 0.0, 0.0]):
    """
    Translates the point.

    Parameters
    ----------
    displacement: list of floats, optional
      Displacement in each direction;
      default: [0.0, 0.0, 0.0].
    """
    self.x += displacement[0]
    self.y += displacement[1]
    if self.dimensions == 3:
      self.z += displacement[2]


class Geometry(object):
  """
  Contains information about a geometry.
  """
  dimensions = None

  def __init__(self, points=None, file_path=None):
    """
    Initializes the geometry with points.

    Parameters
    ----------
    points: list of Point objects, optional
      List of points that defines the geometry;
      default: None.
    file_path: string, optional
      Path of the file with coordinates;
      default: None.
    """
    if points:
      self.points = points
      self.points_initial = copy.deepcopy(points)
    if file_path:
      self.read_from_file(file_path)
    if self.dimensions or self.points:
      self.get_mass_center()

  def read_from_file(self, file_path):
    """
    Reads the coordinates of the geometry from a file.

    Parameters
    ----------
    file_path: string
      Path of the file that contains list of coordinates.
    """
    if not self.dimensions:
      print('\nGet the dimension of the geometry ...')
      with open(file_path, 'r') as infile:
        dim = len(infile.readlines()[1].strip().split())
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
    """
    Gathers a given component of all points into a Numpy array.

    Parameters
    ----------
    component: string
      Component of the point to look for ('x', 'y' or 'z').
    position: string, optional
      Position of the body ('current' or 'initial');
      default: 'current'.

    Returns
    -------
    array: 1D array of floats
      Array with the appropriate component of all points defining the geometry.
    """
    if position == 'current':
      return numpy.array([getattr(point, component)
                          for point in self.points])
    elif position == 'initial':
      return numpy.array([getattr(point, component)
                          for point in self.points_initial])

  def broadcast_coordinate(self, array, component):
    """
    Broadcasts array elements to point coordinate.

    Parameters
    ----------
    array: 1D array of floats
      Array to broadcast to a component of all points.
    component: string
      Point's component to be filled.
    """
    for i, point in enumerate(self.points):
      setattr(self.points[i], component, array[i])

  def get_mass_center(self):
    """
    Computes the center of mass of the geometry.
    """
    x_mass = self.gather_coordinate('x').mean()
    y_mass = self.gather_coordinate('y').mean()
    z_mass = (self.gather_coordinate('z').mean()
              if self.dimensions == 3 else None)
    self.mass_center = Point(x_mass, y_mass, z_mass)

  def translation(self, displacement=[0.0, 0.0, 0.0]):
    """
    Translates the geometry.

    Parameters
    ----------
    displacement: list of floats, optional
      Displacement in each direction;
      default: [0.0, 0.0, 0.0].
    """
    if not any(displacement):
      return
    print('\nTranslate the geometry ...')
    for i, point in enumerate(self.points):
      self.points[i].translation(displacement)
    self.get_mass_center()

  def rotation(self, center=None,
               roll=0.0, yaw=0.0, pitch=0.0, mode='deg'):
    """
    Rotates the geometry.

    Parameters
    ----------
    center: Point object, optional
      Center of rotation;
      default: None.
    roll, yaw, pitch: floats, optional
      Angles of rotation;
      default: 0.0, 0.0, 0.0.
    mode: string, optional
      Angles in degrees ('deg') or radians ('rad');
      default: 'deg'.
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
    """
    Scales the geometry.

    Parameters
    ----------
    ratio: float, optional
      Scaling ratio;
      default: 1.0.
    """
    if ratio == 1.0:
      return
    print('\nScale the geometry ...')
    components = (['x', 'y', 'z'] if self.dimensions == 3 else ['x', 'y'])
    for component in components:
      center = getattr(self.mass_center, component)
      s = center + ratio * (self.gather_coordinate(component) - center)
      self.broadcast_coordinate(s, component)

  def write(self, file_path=os.path.join(os.getcwd(), 'new_body')):
    """
    Writes the coordinates into a file.

    Parameters
    ----------
    file_path: string, optional
      Path of the output file;
      default: ./new_body.
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
  """
  Contains information about a two-dimensional geometry.
  """
  dimensions = 2

  def __init__(self, points=None, file_path=None):
    """
    Initializes the geometry.

    Parameters
    ----------
    points: list of Point object, optional
      List of points that defines the geometry;
      default: None.
    file_path: string, optional
      Path of the file with coordinates;
      default: None.
    """
    Geometry.__init__(self, points, file_path)

  def perimeter(self):
    """
    Returns the perimeter of the closed geometry.

    Returns
    -------
    perimeter: float
      Perimeter of the closed geometry.
    """
    x, y = self.gather_coordinate('x'), self.gather_coordinate('y')
    x, y = numpy.append(x, x[0]), numpy.append(y, y[0])
    return numpy.sum(numpy.sqrt((x[1:] - x[:-1])**2 + (y[1:] - y[:-1])**2))

  def extrusion(self, limits=[-0.5, 0.5], n=None, ds=None, force=False):
    """
    Extrudes the two-dimensional geometry in the z-direction.

    Parameters
    ----------
    limits: 2-list of float, optional
      Limits of the extrusion;
      default: [-0.5, 0.5].
    n: integer, optional
      Number of divisions in the z-direction;
      default: None.
    ds: float, optional
      Desired segment-length;
      default: None.
    force: boolean, optional
      Forces the extrusion to the limits prescribed;
      default: False.

    Returns
    -------
    geometry3d: Geometry3d object
      An instance of a three dimensional geometry.
    """
    print('\nExtrude the geometry in the z-direction...')
    if not (ds or n):
      raise ValueError('both ds and n are set to None')
    elif abs(limits[0] - limits[1]) < 1.0E-06:
      raise ValueError('limits are too close from each other')
    z_start, z_end = limits
    if not n:
      n = int(math.ceil(abs(z_start - z_end) / ds))
    ds = abs(z_start - z_end) / n
    s = math.copysign(1.0, z_end - z_start)
    if force:
      z = numpy.linspace(z_start, z_end, n + 1)
    else:
      z = numpy.linspace(z_start + s * 0.5 * ds, z_end - s * 0.5 * ds, n)
    points = sum(([Point(point.x, point.y, z)
                   for point in self.points] for z in z), [])
    return Geometry3d(points)

  def discretization(self, n=None, ds=None):
    """
    Discretizes the geometry.

    Parameters
    ----------
    n: integer, optional
      Number of divisions;
      default: None.
    ds: float, optional
      Desired segment-length;
      default: None.
    """
    if not (n or ds):
      return
    print('\nDiscretize the geometry ...')
    if not n:
      n = int(math.ceil(self.perimeter() / ds))
    ds = self.perimeter() / n
    # store initial points
    points_old = copy.deepcopy(self.points) + [self.points[0]]
    last = len(points_old) - 1
    # initialize new list of points
    self.points = [points_old[0]]
    # compute new points
    next = 1
    tolerance = 1.0E-06
    for i in xrange(1, n):
      start = self.points[-1]
      end = points_old[next]
      distance = start.distance(end)
      # copy
      if abs(ds - distance) <= tolerance:
        self.points.append(end)
        next += 1
      # interpolation
      elif ds < distance:
        length = start.distance(end)
        start, end = start.as_array(), end.as_array()
        new = start + ds / length * (end - start)
        self.points.append(Point(*new))
      # projection
      else:
        # get segment index
        while distance < ds and next < last:
          next += 1
          distance = start.distance(points_old[next])
        # gather coordinates as array
        previous = points_old[next - 1].as_array()
        end = points_old[next].as_array()
        # interpolate on segment
        precision = 1
        coeff = 0.0
        while abs(ds - distance) > tolerance:
          new = previous + coeff * (end - previous)
          distance = start.distance(Point(*new))
          if distance > ds:
            coeff -= 0.1**precision
            precision += 1
          coeff += 0.1**precision
        # check point not too close from first point before adding
        if self.points[0].distance(Point(*new)) > 0.5 * ds:
          self.points.append(Point(*new))

  def plot(self):
    """
    Plots the two-dimensional geometry using Matplotlib.
    """
    print('\nPlot the two-dimensional geometry ...')
    pyplot.style.use(os.path.join(os.environ['PETIBM_DIR'],
                                  'scripts',
                                  'python',
                                  'style',
                                  'style_PetIBM.mplstyle'))
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
  """
  Contains information about a line.
  """
  def __init__(self, start=Point(0.0, 0.0), length=1.0,
               n=None, ds=None):
    """
    Creates the line.

    Parameters
    ----------
    start: Point object, optional
      Starting point of the line;
      default: Point(0.0, 0.0).
    length: float, optional
      Length of the line;
      default: 1.0.
    n: integer, optional
      Number of divisions on the line;
      default: None.
    ds: float, optional
      Desired segment-length;
      default: None.
    """
    self.start = start
    self.length = length
    self.n, self.ds = n, ds
    self.create()
    self.get_mass_center()

  def create(self):
    """
    Creates the points on the line.
    """
    print('\nCreate a line ...')
    if not (self.n or self.ds):
      raise ValueError('both ds and n are set to None')
    elif not self.n:
      self.n = int(math.ceil(self.length / self.ds))
    x = self.start.x + numpy.linspace(0.0, self.length, self.n + 1)
    y = self.start.y * numpy.ones(self.n + 1)
    self.points = [Point(x[i], y[i]) for i in xrange(self.n + 1)]
    self.points_initial = copy.deepcopy(self.points)


class Circle(Geometry2d):
  """
  Contains information about a circular geometry.
  """
  def __init__(self, center=Point(0.0, 0.0), radius=0.5,
               n=None, ds=None):
    """Creates the circle.

    Parameters
    ----------
    center: Point object, optional
      Center of the circle;
      default: Point(0.0, 0.0).
    radius: float, optional
      Radius of the circle;
      default: 0.5.
    n: integer, optional
      Number of divisions;
      default: None.
    ds: float, optional
      Desired segment-length;
      default: None.
    """
    self.center = center
    self.radius = radius
    self.n, self.ds = n, ds
    self.create()
    self.get_mass_center()

  def create(self):
    """
    Creates the points on the circle.
    """
    print('\nCreate a circle ...')
    if not (self.n or self.ds):
      raise ValueError('both ds and n are set to None')
    elif not self.n:
      self.n = int(math.ceil(2.0 * math.pi * self.radius / self.ds))
    theta = numpy.linspace(0.0, 2.0 * math.pi, self.n + 1)
    x = self.center.x + self.radius * numpy.cos(theta)[:-1]
    y = self.center.y + self.radius * numpy.sin(theta)[:-1]
    self.points = [Point(x[i], y[i]) for i in xrange(x.size)]
    self.points_initial = copy.deepcopy(self.points)


class Rectangle(Geometry2d):
  """
  Contains information about a rectangle.
  """
  def __init__(self,
               bottom_left=Point(0.0, 0.0), top_right=Point(1.0, 1.0),
               nx=None, ny=None, ds=None):
    """
    Creates the rectangle.

    Parameters
    ----------
    bottom_left: Point object, optional
      Bottom-left point;
      default: Point(0.0, 0.0).
    top_right: Point object, optional
      Top-right point;
      default: Point(1.0, 1.0).
    nx: integer, optional
      Number of divisions in the x-direction;
      default: None.
    ny: integer, optional
      Number of divisions in the y-direction;
      default: None.
    ds: float, optional
      Desired segment-length;
      default: None.
    """
    self.bottom_left, self.top_right = bottom_left, top_right
    self.nx, self.ny, self.ds = nx, ny, ds
    self.create()
    self.get_mass_center()

  def create(self):
    """
    Creates the points on the rectangle.
    """
    print('\nCreate a rectangle ...')
    if self.ds:
      self.nx = int(math.ceil((self.top_right.x - self.bottom_left.x)
                              / self.ds))
      self.ny = int(math.ceil((self.top_right.y - self.bottom_left.y)
                              / self.ds))
    elif not (self.ds or (self.nx and self.ny)):
      raise ValueError('ds is set to None '
                       'while nx and/or ny are set to None')
    self.points = []
    # bottom
    x_bottom = numpy.linspace(self.bottom_left.x,
                              self.top_right.x,
                              self.nx + 1)[:-1]
    self.points += [Point(x, self.bottom_left.y) for x in x_bottom]
    # right
    y_right = numpy.linspace(self.bottom_left.y,
                             self.top_right.y,
                             self.ny + 1)[:-1]
    self.points += [Point(self.top_right.x, y) for y in y_right]
    # top
    x_top = numpy.linspace(self.top_right.x,
                           self.bottom_left.x,
                           self.nx + 1)[:-1]
    self.points += [Point(x, self.top_right.y) for x in x_top]
    # left
    y_left = numpy.linspace(self.top_right.y,
                            self.bottom_left.y,
                            self.ny + 1)[:-1]
    self.points += [Point(self.bottom_left.x, y) for y in y_left]
    self.points_initial = copy.deepcopy(self.points)


class Geometry3d(Geometry):
  """
  Contains information about a three-dimensional geometry.
  """
  dimensions = 3

  def __init__(self, points=None, file_path=None):
    """
    Initializes the three-dimensional geometry with points.

    Parameters
    ----------
    points: list of Point objects, optional
      List of points that defines the geometry;
      default: None.
    file_path: string, optional
      Path of the file with coordinates;
      default: None.
    """
    Geometry.__init__(self, points, file_path)

  def plot(self):
    """
    Plots the geometry using the package Mayavi.
    """
    print('\nPlot the three-dimensional geometry ...')
    from mayavi import mlab
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
  """
  Contains information about a spherical geometry.
  """
  def __init__(self, center=Point(0.0, 0.0, 0.0), radius=0.5,
               n=None, ds=None):
    """
    Creates the sphere.

    Parameters
    ----------
    center: Point object, optional
      Center of the sphere;
      default: Point(0.0, 0.0, 0.0).
    radius: float, optional
      Radius of the sphere;
      default: 0.5.
    n: integer, optional
      Number of divisions on the great-circle of the sphere;
      default: None.
    ds: float, optional
      Desired segment length;
      default: None.
    """
    self.center = center
    self.radius = radius
    self.n, self.ds = n, ds
    self.create()
    self.get_mass_center()

  def create(self):
    """
    Creates the points on the sphere.
    """
    print('\nCreate a sphere ...')
    if not (self.n or self.ds):
      raise ValueError('both ds and n are set to None')
    elif self.n and not self.ds:
      self.ds = 2.0 * math.pi * self.radius / self.n
    n_phi = int(math.ceil(math.pi * self.radius / self.ds))
    phi = numpy.linspace(0.0, math.pi, n_phi)[1:-1]
    # north pole
    x = self.center.x
    y = self.center.y
    z = self.center.z + self.radius
    # between poles
    for phi in phi:
      rsinphi = self.radius * math.sin(phi)
      rcosphi = self.radius * math.cos(phi)
      n_theta = int(math.ceil(2.0 * math.pi * rsinphi / self.ds))
      theta = numpy.linspace(0.0, 2.0 * math.pi, n_theta + 1)[:-1]
      x = numpy.append(x, self.center.x + rsinphi * numpy.cos(theta))
      y = numpy.append(y, self.center.y + rsinphi * numpy.sin(theta))
      z = numpy.append(z, self.center.z + rcosphi * numpy.ones(n_theta))
    # south pole
    x = numpy.append(x, self.center.x)
    y = numpy.append(y, self.center.y)
    z = numpy.append(z, self.center.z - self.radius)
    # create points
    self.points = [Point(x[i], y[i], z[i]) for i in xrange(x.size)]
    self.points_initial = copy.deepcopy(self.points)
