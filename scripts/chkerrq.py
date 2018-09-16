"""Merge CHKERRQ(ierr) with the previous statement on a single line."""

import os
import argparse
import pathlib


def parse_command_line():
    """Parse the command-line options."""
    formatter_class = argparse.ArgumentDefaultsHelpFormatter
    description = 'Clang-format: Allow CHKERRQ to be on same line.'
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=formatter_class)
    parser.add_argument('--version', '-V',
                        action='version',
                        version='%(prog)s (version 0.1)')
    parser.add_argument('--file', '-f', dest='filepath',
                        type=str,
                        required=True,
                        help='Path of the file.')
    parser.add_argument('--output', '-o', dest='outpath',
                        type=str,
                        default=None,
                        help='Path of the output file.'
                             'Default adds "-new" to the input file name.')
    parser.add_argument('--inplace', dest='inplace',
                        action='store_true',
                        default=False,
                        help='Use flag to modify inplace.')
    args = parser.parse_args()
    if args.outpath is None and not args.inplace:
        name, ext = os.path.splitext(args.filepath)
        args.outpath = name + '-new' + ext
    if args.inplace:
        args.outpath = args.filepath
    return args


def get_key_in_line(line, keys):
    for key in keys:
        if key in line:
            break
        key = None
    return key


def main(args):
    """Merge CHKERRQ(ierr) with the previous statement on a single line."""
    filepath = pathlib.Path(args.filepath)
    with open(filepath, 'r') as infile:
        lines = infile.readlines()
    newlines = []
    keys = ('  CHKERRQ(ierr);', '  CHKERRV(ierr);')
    for line in lines:
        key = get_key_in_line(line, keys)
        if key is None:
            newlines.append(line)
        else:
            candidate = newlines[-1].rstrip() + key[1:] + '\n'
            if len(candidate) > 80:
                newlines.append(line)
            else:
                newlines[-1] = candidate
    filepath = pathlib.Path(args.outpath)
    with open(filepath, 'w') as outfile:
        for line in newlines:
            outfile.write(line)


if __name__ == '__main__':
    args = parse_command_line()
    main(args)
