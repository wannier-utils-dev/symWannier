#!/usr/bin/env python
"""
Expand Wannier input files using symmetry information:
prefix.immn, prefix.iamn, prefix.ieig, prefix.isym => prefix.mmn, prefix.amn, prefix.eig
"""

import sys
import argparse
import logging

from symwannier.nnkp import Nnkp
from symwannier.sym import Sym
from symwannier.amn import Amn
from symwannier.mmn import Mmn
from symwannier.eig import Eig

def main(argv=None, for_cli=False):
    """Expand symmetry-reduced Wannier input files to full k-point sets.

    Parameters
    ----------
    argv : list, optional
        Command-line arguments.
    for_cli : bool, optional
        Whether called from CLI.
    """
    if argv is None:
        argv = sys.argv[1:]
    progname = "symmwanier expand" if for_cli else "python expand_wannier_inputs.py"

    logging.basicConfig(level=logging.INFO, format="%(message)s")
    log = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(
        prog=progname,
        description="Expand Wannier input files using symmetry information."
    )

    parser.add_argument(
        "prefix",
        help="Prefix name of input/output files"
    )

    prefix = parser.parse_args(argv).prefix

    nnkp = Nnkp(file_nnkp=prefix+".nnkp", log=log)
    sym = Sym(file_sym=prefix+".isym", nnkp=nnkp, log=log)

    # Eig
    log.info(f"{prefix}.ieig => {prefix}.eig")
    eig = Eig(file_eig=prefix+".ieig", sym=sym, log=log)
    eig.write_eig(prefix+".eig")

    # Amn
    log.info(f"{prefix}.iamn => {prefix}.amn")
    amn = Amn(file_amn=prefix+".iamn", sym=sym, nnkp=nnkp, log=log)
    amn.write_amn(prefix+".amn")

    # Mmn
    log.info(f"{prefix}.immn => {prefix}.mmn")
    mmn = Mmn(file_mmn=prefix+".immn", nnkp=nnkp, sym=sym, log=log)
    mmn.write_mmn(prefix+".mmn")

if __name__ == '__main__':
    main()
