#!/usr/bin/env python

import numpy as np
import scipy.linalg
import os
import itertools

from symwannier.nnkp import Nnkp
from symwannier.sym import Sym

class Eig():
    """
    eig file: e_n(k)
    """
    def __init__(self, file_eig, sym = None):
        self.sym = sym

        if os.path.exists(file_eig):
            dat = np.genfromtxt(file_eig)
            num_bands = int(dat[-1,0])
            nk = int(dat[-1,1])
            dat = dat.reshape(nk, num_bands, 3)
            eig = dat[:,:,2]
        else:
            raise Exception("failed to read eig file: " + file_eig)

        if self.sym is None:
            self.nk = nk
            self.num_bands = num_bands
            self.eig = eig

        else:
            self.nk = self.sym.nkf
            self.num_bands = num_bands
            self.eig = np.zeros([self.nk, self.num_bands])
            for ik, k in enumerate(self.sym.full_kpoints):
                iks = self.sym.equiv[ik]
                self.eig[ik,:] = eig[iks,:]

    def write_eig(self, file_eig):
        with open(file_eig, "w") as fp:
            for ik, n in itertools.product( range(self.nk), range(self.num_bands) ):
                fp.write("{:5d}{:5d}{:18.12f}\n".format(n+1, ik+1, self.eig[ik,n]))
