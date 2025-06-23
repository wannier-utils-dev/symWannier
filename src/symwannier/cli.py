#!/usr/bin/env python
"""
symwannier - Symmetry-adapted Wannier tools
"""

import argparse
import sys
from symwannier import expand_wannier_inputs, wannierize

def main():
    parser = argparse.ArgumentParser(
        prog="symwannier",
        description="symwannier - Symmetry-adapted Wannier tools",
    )

    subparsers = parser.add_subparsers(dest="command", help="Subcommands")

    expand_parser = subparsers.add_parser("expand", help="Expand input files using symmetry information", add_help=False)

    wannier_parser = subparsers.add_parser("wannierize", help="Perform Wannierization using symmetry inputs", add_help=False)

    args, remainder = parser.parse_known_args()

    if args.command == "expand":
        expand_wannier_inputs.main(remainder, for_cli=True)
    elif args.command == "wannierize":
        wannierize.main(remainder, for_cli=True)
    else:
        parser.print_help()
