#!/usr/bin/env python

# file: generateBody.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Generates a body file.


import os
import sys
import argparse

sys.path.append('{}/scripts/python'.format(os.environ['PETIBM_DIR']))
import geometry


def read_inputs():
  """Parses the command-line."""
  # create parser
  parser = argparse.ArgumentParser(description='Geometry discretization',
                        formatter_class= argparse.ArgumentDefaultsHelpFormatter)
  # fill parser with arguments
  # geometry arguments
  parser.add_argument('--type', dest='body_type', type=str,
                      help='type of body '
                           '(file, circle, line, rectangle, sphere)')
  parser.add_argument('--file', '-f', dest='file_path', type=str,
                      help='path of the coordinates file')
  parser.add_argument('--circle', dest='circle', type=float, nargs='+',
                      default=[0.5, 0.0, 0.0],
                      help='radius and center-coordinates of the circle')
  parser.add_argument('--line', '-l', dest='line', type=float, nargs='+',
                      default=[1.0, 0.0, 0.0],
                      help='length and starting-point of the line')
  parser.add_argument('--rectangle', dest='rectangle', type=float, nargs='+',
                      default=[0.0, 0.0, 1.0, 1.0],
                      help='bottom-left and top-right coordinates')
  parser.add_argument('--sphere', dest='sphere', type=float, nargs='+',
                      default=[0.5, 0.0, 0.0, 0.0],
                      help='radius and center-coordinates of the sphere')
  # discretization arguments
  parser.add_argument('--n', '-n', dest='n', type=int, default=None,
                      help='number of divisions')
  parser.add_argument('--ds', '-ds', dest='ds', type=float, default=None,
                      help='target segment-length')
  # geometry modification arguments
  parser.add_argument('--rotation', '-r', dest='rotation', type=float, 
                      nargs='+', default=None,
                      help='center of rotation')
  parser.add_argument('--roll', dest='roll', type=float, default=0.0,
                      help='roll angle')
  parser.add_argument('--yaw', dest='yaw', type=float, default=0.0,
                      help='yaw angle')
  parser.add_argument('--pitch', dest='pitch', type=float, default=0.0,
                      help='pitch angle')
  parser.add_argument('--mode', dest='mode', type=str, default='deg',
                      help='angles in degrees(deg) or radians (rad)')
  parser.add_argument('--translation', '-t', dest='translation', type=float,
                      nargs='+', default=[0.0, 0.0, 0.0],
                      help='displacement in the x-, y- and z- directions')
  parser.add_argument('--scale', '-s', dest='scale', type=float, default=1.0,
                      help='scaling factor for 2D geometry')
  parser.add_argument('--extrusion', '-e', dest='extrusion', type=float, 
                      nargs='+',
                      help='limits of the cylinder in the third direction')
  # output arguments
  parser.add_argument('--save-name', dest='save_name', type=str, 
                      default='new_body', 
                      help='name of the new body file')
  parser.add_argument('--extension', dest='extension', type=str, default='body',
                      help='extension of the output file')
  parser.add_argument('--save-dir', dest='save_directory', type=str, 
                      default=os.getcwd(),
                      help='directory where body file will be saved')
  parser.add_argument('--no-save', dest='save', action='store_false',
                      help='does not save the geometry into a file')
  parser.add_argument('--show', dest='show', action='store_true',
                      help='displays the geometry')
  parser.set_defaults(save=True)
  return parser.parse_args()


def main():
  """Generates a file containing the coordinates of a body."""
  # parse command-line
  args = read_inputs()

  # generate the geometry
  if args.body_type == 'file':
    body = geometry.Geometry(file_path=args.file_path)
  elif args.body_type == 'circle':
    body = geometry.Circle(radius=args.circle[0], 
                           center=geometry.Point(args.circle[1], args.circle[2]),
                           n=args.n, ds=args.ds)
  elif args.body_type == 'line':
    body = geometry.Line(length=args.line[0],
                         start=geometry.Point(args.line[1], args.line[2]),
                         n=args.n, ds=args.ds)
  elif args.body_type == 'rectangle':
    body = geometry.Rectangle(bottom_left=geometry.Point(args.rectangle[0],
                                                         args.rectangle[1]),
                              top_right=geometry.Point(args.rectangle[2],
                                                       args.rectangle[3]),
                              nx=args.n, ny=args.n, ds=args.ds)
  elif args.body_type == 'sphere':
    body = geometry.Sphere(radius=args.sphere[0],
                           center=geometry.Point(args.sphere[1], 
                                                 args.sphere[2],
                                                 args.sphere[3]),
                           n=args.n, ds=args.ds)
  body.scale(ratio=args.scale)
  body.rotation(center=args.rotation, 
                roll=args.roll, yaw=args.yaw, pitch=args.pitch, mode=args.mode)
  body.translation(displacement=args.translation)
  if body.dimensions == 2 and args.body_type == 'file':
    body.discretization(n=args.n, ds=args.ds)
  if body.dimensions == 2 and args.extrusion:
    body = body.extrusion(limits=args.extrusion, n=args.n, ds=args.ds)
  if args.save:
      output_path = '{}/{}.{}'.format(args.save_directory, 
                                      args.save_name, 
                                      args.extension)
      body.write(file_path=output_path)
  if args.show:
      body.plot()

  print('\n[{}] DONE'.format(os.path.basename(__file__)))


if __name__ == '__main__':
  main()