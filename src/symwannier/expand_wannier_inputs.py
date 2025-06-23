#!/usr/bin/env python
"""
 Expand Wannier input files using symmetry information:
 prefix.immn, prefix.iamn, prefix.ieig, prefix.isym => prefix.mmn, prefix.amn, prefix.eig
"""

import sys
import os
import argparse

from symwannier.nnkp import Nnkp
from symwannier.sym import Sym
from symwannier.amn import Amn
from symwannier.mmn import Mmn
from symwannier.eig import Eig

def main(argv=None, for_cli=False):
    if argv is None:
        argv = sys.argv[1:]
    progname = "symmwanier expand" if for_cli else "python expand_wannier_inputs.py"

    parser = argparse.ArgumentParser(
        prog=progname,
        description="Expand Wannier input files using symmetry information."
    )

    parser.add_argument(
        "prefix",
        help="Prefix name of input/output files"
    )

    prefix = parser.parse_args(argv).prefix

    nnkp = Nnkp(file_nnkp=prefix+".nnkp")
    sym = Sym(file_sym=prefix+".isym", nnkp=nnkp)

    # Eig
    print(" {0:s}.ieig => {0:s}.eig".format(prefix))
    eig = Eig(file_eig=prefix+".ieig", sym=sym)
    eig.write_eig(prefix+".eig")

    # Amn
    print(" {0:s}.iamn => {0:s}.amn".format(prefix))
    amn = Amn(file_amn=prefix+".iamn", sym=sym, nnkp=nnkp)
    amn.write_amn(prefix+".amn")

    # Mmn
    print(" {0:s}.immn => {0:s}.mmn".format(prefix))
    mmn = Mmn(file_mmn=prefix+".immn", nnkp=nnkp, sym=sym)
    mmn.write_mmn(prefix+".mmn")

if __name__ == '__main__':
    main()
