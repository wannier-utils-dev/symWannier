#!/usr/bin/env python

import numpy as np
import os
import itertools
import logging

class Eig():
    """Reader for eig files (eigenvalues e_n(k))."""
    def __init__(self, file_eig, sym = None, log=None):
        """Load eigenvalue data from file.

        Parameters
        ----------
        file_eig : str
            Path to eig file.
        sym : Sym, optional
            Symmetry data for IBZ expansion.
        log : logging.Logger, optional
            Logger instance.
        """
        self.log = log or logging.getLogger(__name__)
        if not self.log.handlers:
            logging.basicConfig(level=logging.INFO, format="%(message)s")
        self.sym = sym

        if os.path.exists(file_eig):
            dat = np.genfromtxt(file_eig)
            num_bands = int(dat[-1,0])
            nk = int(dat[-1,1])
            dat = dat.reshape(nk, num_bands, 3)
            eig = dat[:,:,2]
        else:
            raise FileNotFoundError(f"eig file not found: {file_eig}")

        if self.sym is None:
            self.nk = nk
            self.num_bands = num_bands
            self.eig = eig

        else:
            self.nk = self.sym.nkf
            self.num_bands = num_bands
            self.eig = eig[self.sym.equiv, :]

    def write_eig(self, file_eig):
        """Write eigenvalues to file in wannier90 format.

        Parameters
        ----------
        file_eig : str
            Output file path.
        """
        with open(file_eig, "w") as fp:
            for ik, n in itertools.product( range(self.nk), range(self.num_bands) ):
                fp.write("{:5d}{:5d}{:18.12f}\n".format(n+1, ik+1, self.eig[ik,n]))
