#!/usr/bin/env python

import sys
import os
from symwannier.nnkp import Nnkp
from symwannier.sym import Sym
from symwannier.amn import Amn
from symwannier.mmn import Mmn
from symwannier.eig import Eig

prefix = sys.argv[1]

sym = Sym(file_sym=prefix+"_sym.dat")
nnkp = Nnkp(file_nnkp=prefix+".nnkp")

# Eig
print(" {0:s}_ibz.eig => {0:s}.eig".format(prefix))
eig = Eig(file_eig=prefix+"_ibz.eig", sym=sym)
eig.write_eig(prefix+".eig")

# Amn
print(" {0:s}_ibz.amn => {0:s}.amn".format(prefix))
amn = Amn(file_amn=prefix+"_ibz.amn", sym=sym, nnkp=nnkp)
amn.write_amn(prefix+".amn")

# Mmn
print(" {0:s}_ibz.mmn => {0:s}.mmn".format(prefix))
mmn = Mmn(file_mmn=prefix+"_ibz.mmn", nnkp=nnkp, sym=sym)
mmn.write_mmn(prefix+".mmn")
