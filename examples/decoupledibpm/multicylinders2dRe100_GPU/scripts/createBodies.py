"""
Create two circles and write coordinates into files.
"""

import pathlib
import math
import numpy


def create_circle(R=0.5, center=(0.0, 0.0), ds=0.01):
    """
    Generate coordinates of circle.

    Parameters
    ----------
    R: float, optional
        radius;
        default: 0.5.
    center: 2-tuple of floats, optional
        center's coordinates;
        default: (0.0, 0.0).
    ds: float, optional
        mesh spacing;
        default: 0.01.

    Returns
    -------
    x: 1D Numpy array of floats
        x coordinates.
    y: 1D Numpy array of floats
        y coordinates.
    """
    xc, yc = center
    n = math.ceil(2 * numpy.pi * R / ds)
    theta = numpy.linspace(0.0, 2.0 * numpy.pi, num=n, endpoint=False)
    x, y = xc + R * numpy.cos(theta), yc + R * numpy.sin(theta)
    return x, y


def write_coordinates(coords, filepath):
    """
    Write coordinates into file.
    The first line contains the number of points.
    The coordinates are written into columns.

    Parameters
    ----------
    coords: 2-tuple of 1D Numpy arrays of floats
        The coordinates.
    filepath: string
        Path of the output file.
    """
    x, y = coords
    n = x.size
    with open(filepath, 'w') as outfile:
        outfile.write('{}\n'.format(n))
    with open(filepath, 'ab') as outfile:
        numpy.savetxt(outfile, numpy.c_[x, y])


simu_dir = pathlib.Path(__file__).absolute().parents[1]

# Generate coordinates of first body and write into file.
circle1 = create_circle(R=0.5, center=(0.0, -2.5), ds=0.02)
filepath = simu_dir / 'circle1.body'
write_coordinates(circle1, filepath)

# Generate coordinates of first body and write into file.
circle2 = create_circle(R=0.5, center=(0.0, 2.5), ds=0.02)
filepath = simu_dir / 'circle2.body'
write_coordinates(circle2, filepath)
