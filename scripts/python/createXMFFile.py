"""
Creates a single XMF file that includes info for each time-step saved.
"""

import os
import argparse
from lxml import etree as ET

import numpy


def parse_command_line():
  """Parses the command-line."""
  print('\nParsing command-line ...'),
  # create parser
  formatter_class = argparse.ArgumentDefaultsHelpFormatter
  parser = argparse.ArgumentParser(description='Plots the instantaneous '
                                               'forces',
                                   formatter_class=formatter_class)
  # fill parser with arguments
  parser.add_argument('--directory', dest='directory',
                      type=str,
                      default=os.getcwd(),
                      metavar='<directory>',
                      help='directory of the simulation')
  parser.add_argument('--grid-file', dest='grid_file',
                      type=str,
                      metavar='<fileName>',
                      help='name of the grid file')
  parser.add_argument('--grid-size', dest='grid_size',
                      type=int,
                      nargs='+',
                      help='grid size')
  parser.add_argument('--variables', dest='variables',
                      type=str,
                      nargs='+',
                      default=('var1', 'var2'),
                      metavar=('<var1>', '<var2>'),
                      help='list of variables to include')
  parser.add_argument('--times-file', dest='times_file',
                      type=str,
                      metavar='<fileName>',
                      help='name of the file containing the map '
                           '(index, time-unit)')
  parser.add_argument('--outfile', dest='outfile',
                      type=str,
                      default='example.xmf',
                      metavar='<fileName>',
                      help='name of output .xmf file')
  parser.add_argument('--options',
                      type=open, action=ReadOptionsFromFile,
                      metavar='<filePath>',
                      help='path of the file with options to parse')
  print('done')
  # return namespace
  return parser.parse_args()


class ReadOptionsFromFile(argparse.Action):
  """
  Container to read parameters from file.
  """

  def __call__(self, parser, namespace, values, option_string=None):
    """
    Fills the namespace with parameters read in file.
    """
    with values as infile:
      lines = [element for line in infile.readlines()
               for element in line.strip().split()
               if not line.startswith('#')]
      lines = [os.path.expandvars(line) if '$' in line else line
               for line in lines[:]]
    parser.parse_args(lines, namespace)


def main():
  """
  Creates an XMF file with support for multiple variables
  (that share the same grid).
  """
  # get the user configuration
  args = parse_command_line()
  print('User-configuration:\n{}'.format(args))
  # read file containing index and corresponding
  # time at which solution was saved
  with open(args.times_file, 'r') as infile:
    time_values = numpy.genfromtxt(infile, dtype=None)
  # start xmf tree
  xdmf = ET.Element('Xdmf',
                    Version='2.2')
  info = ET.SubElement(xdmf, 'Information',
                       Name='MetaData',
                       Value='ID-23454')
  domain = ET.SubElement(xdmf, 'Domain')
  grid_time_series = ET.SubElement(domain, 'Grid',
                                   Name='TimeSeries',
                                   GridType='Collection',
                                   CollectionType='Temporal')
  # define type of field
  if len(args.grid_size) == 2:
    topology_type = '2DRectMesh'
    geometry_type = 'VXVY'
    components = ('x', 'y')
  elif len(args.grid_size) == 3:
    topology_type = '3DRectMesh'
    geometry_type = 'VXVYVZ'
    components = ('x', 'y', 'z')
  number_of_elements = ' '.join(str(n) for n in args.grid_size)
  # create an xmf block for each time-step saved
  for it, time_value in time_values:
    grid = ET.SubElement(grid_time_series, 'Grid',
                         Name='Grid',
                         GridType='Uniform')
    time = ET.SubElement(grid, 'Time',
                         Value=str(time_value))
    topology = ET.SubElement(grid, 'Topology',
                             TopologyType=topology_type,
                             NumberOfElements=number_of_elements)
    geometry = ET.SubElement(grid, 'Geometry',
                             GeometryType=geometry_type)
    # loop over the 3 directions (for code-reuse purpose)
    for d, n in zip(components, args.grid_size):
      dataitem = ET.SubElement(geometry, 'DataItem',
                               Dimensions=str(n),
                               NumberType='Float',
                               Precision='4',
                               Format='HDF')
      dataitem.text = '{}:/{}'.format(os.path.join(args.directory,
                                                   args.grid_file), d)
    # create a block for each variable to insert
    for variable in args.variables:
      attribute = ET.SubElement(grid, 'Attribute',
                                Name=variable,
                                AttributeType='Scalar',
                                Center='Node')
      dataitem = ET.SubElement(attribute, 'DataItem',
                               Dimensions=number_of_elements,
                               NumberType='Float',
                               Precision='4',
                               Format='HDF')
      variable_file_path = os.path.join(args.directory,
                                        '{:0>7}'.format(it),
                                        '{}.h5'.format(variable))
      dataitem.text = '{}:/{}'.format(variable_file_path, variable)
  # write the xmf file
  print('\nWriting XMF file: {} ...'.format(args.outfile)),
  tree = ET.ElementTree(xdmf)
  tree.write(args.outfile, pretty_print=True, xml_declaration=True)
  print('done')


if __name__ == '__main__':
  print('\n[{}] START\n'.format(os.path.basename(__file__)))
  main()
  print('\n[{}] END\n'.format(os.path.basename(__file__)))
