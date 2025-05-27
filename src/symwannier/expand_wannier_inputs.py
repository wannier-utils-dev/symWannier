#!/usr/bin/env python
"""
 Expand Wannier input files using symmetry information:
 prefix.immn, prefix.iamn, prefix.ieig, prefix.isym => prefix.mmn, prefix.amn, prefix.eig
"""

DESCRIPTION = __doc__.strip()
USAGE_CLI = """
Usage:
  expand <prefix>
"""

USAGE_DIRECT = """
Usage:
  expand_wannier_inputs.py <prefix>
"""

import sys
import os
from docopt import docopt
from symwannier.nnkp import Nnkp
from symwannier.sym import Sym
from symwannier.amn import Amn
from symwannier.mmn import Mmn
from symwannier.eig import Eig

def main(argv=None, for_cli=False):
    usage = USAGE_CLI if for_cli else USAGE_DIRECT
    doc = f"{DESCRIPTION}\n{usage}"

    args = docopt(doc, argv=argv)

    prefix = args["<prefix>"]

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
