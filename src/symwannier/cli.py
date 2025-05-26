#!/usr/bin/env python
"""
symwannier - Symmetry-adapted Wannier tools

Usage:
  symwannier expand [<args>...]
  symwannier wannierize [<args>...]
  symwannier (-h | --help)

Options:
  -h --help     Show this help message and exit.

Run 'symwannier <command> --help' for more information on a command.
"""

from docopt import docopt
import sys
from symwannier import expand_wannier_inputs, wannierize

def main():
    args = docopt(__doc__, argv=sys.argv[1:], options_first=True)
    cmd_args = sys.argv[2:]

    if args["expand"]:
        expand_wannier_inputs.main(cmd_args, for_cli=True)
    elif args["wannierize"]:
        wannierize.main(cmd_args, for_cli=True)
    else:
        print(__doc__)
