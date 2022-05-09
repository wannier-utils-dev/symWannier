#!/usr/bin/env python

import sys
import os
from symwannier.nnkp import Nnkp
from symwannier.sym import Sym
from symwannier.amn import Amn
from symwannier.mmn import Mmn
from symwannier.eig import Eig

prefix = sys.argv[1]

w90_dir = "w90/"

if os.path.exists(w90_dir):
    os.rename(w90_dir, "old_"+w90_dir)

os.mkdir(w90_dir)

sym = Sym(file_sym=prefix+"_sym.dat")
nnkp = Nnkp(file_nnkp=prefix+".nnkp")

# Eig
eig = Eig(file_eig=prefix+".eig", sym=sym)
eig.write_eig(w90_dir+prefix+".eig")

# Amn
amn = Amn(file_amn=prefix+".amn", sym=sym, nnkp=nnkp)
amn.write_amn(w90_dir+prefix+".amn")

# Mmn
mmn = Mmn(file_mmn=prefix+".mmn", nnkp=nnkp, sym=sym)
mmn.write_mmn(w90_dir+prefix+".mmn")
