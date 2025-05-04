#!/usr/bin/env python

import numpy as np
import os
import gzip
import itertools

from symwannier.nnkp import Nnkp
from symwannier.sym import Sym

class Mmn:
    def __init__(self, file_mmn, nnkp, sym=None):
        self.nnkp = nnkp
        self.sym = sym

        if os.path.exists(file_mmn):
            with open(file_mmn) as fp:
                self._read_mmn(fp)
        elif os.path.exists(file_mmn + ".gz"):
            with gzip.open(file_mmn + ".gz", 'rt') as fp:
                self._read_mmn(fp)
        else:
            raise Exception("failed to read mmn file: " + file_mmn)

        self._mmn_full_klist()

    def write_mmn(self, file_mmn):
        with open(file_mmn, "w") as fp:
            fp.write("Mmn created by mmn.py\n")
            fp.write("{} {} {}\n".format(self.num_bands, self.nk, self.nb))
            for ik, ib in itertools.product(range(self.nk), range(self.nb)):
                k = self.sym.full_kpoints[ik]
                b = self.nnkp.bvec_crys[ib]
                ikb = self.sym.search_ik_full(k+b)
                g = k + b - self.sym.full_kpoints[ikb]
                fp.write("{0}  {1}  {2[0]}  {2[1]}  {2[2]}\n".format(ik+1, ikb+1, np.round(g).astype("int")))
                for m, n in itertools.product( range(self.num_bands), repeat=2 ):
                    fp.write("{0.real:18.12f}  {0.imag:18.12f}\n".format(self.mmn[ik,ib,n,m]))

    def _read_mmn(self, fp):
        """
        read mmn and set variables
        num_bands:  number of bands
        nks:        number of irreducible k-points
        nb:         number of b-vectors
        mmn:        Mmn
        kb2k:       index of k+b
        kpb_info:   information about k+b
        """
        ibz = "IBZ" in fp.readline()  # first line starts with "IBZ" when prefix_ibz.mmn
        if self.sym is not None and not ibz:
            raise Exception("Mmn is not for IBZ")
        if self.sym is None and ibz:
            raise Exception("IBZ Mmn but no symmetry information.")

        if ibz:
            print("  Reading IBZ mmn file")
        else:
            print("  Reading mmn file")

        self.num_bands, self.nks, self.nb = [ int(x) for x in fp.readline().split() ]

        self.mmn = np.zeros([self.nks, self.nb, self.num_bands, self.num_bands], dtype=complex)
        self.kpb_info = np.zeros([self.nks, self.nb, 5], dtype=int)

        for ik, ib in itertools.product(range(self.nks), range(self.nb)):
            d = [ int(x) for x in fp.readline().split() ]
            assert ik == d[0]-1, "{} {}".format(ik, d[0])
            self.kpb_info[ik,ib,:] = d
            for m, n in itertools.product( range(self.num_bands), repeat=2 ):
                dat = [ float(x) for x in fp.readline().split() ]
                self.mmn[ik,ib,n,m] = dat[0] + 1j*dat[1]

    def _mmn_full_klist(self):
        if self.sym is None:
            assert self.nks == self.nnkp.nk
        else:
            assert self.sym.nkf == self.nnkp.nk
        self.nk = self.nnkp.nk   # nk for full k list
        assert self.nb == self.nnkp.nb

        self.kb2k = - np.ones([self.nk, self.nb], dtype=int)  # -1 becomes non-negative when defined
        ####### simple case (without symmetry) #######
        if self.sym is None:
            for ik, ib in itertools.product(range(self.nk), range(self.nb)):
                assert ik == self.kpb_info[ik,ib,0]-1
                bvec = self.nnkp.calc_bvec(self.kpb_info[ik,ib,:])
                ibt = self.nnkp.bvec_num(bvec)
                assert ibt == ib
                self.kb2k[ik,ib] = self.kpb_info[ik,ib,1] - 1

        ####### symmetrized case #######
        else:
            mmn = np.zeros([self.nk, self.nb, self.num_bands, self.num_bands], dtype=complex)

            bvec_equiv = self.sym.kpoint_equiv_info(self.nnkp.bvec_crys)
            for ikf, kf in enumerate(self.sym.full_kpoints):
                # s[isym] . irr_k[iki] = full_k[ikf]
                iki = self.sym.equiv[ikf]
                isym1 = self.sym.equiv_sym[ikf]
                ki = self.sym.irr_kpoints[iki]
                for ibf, bf in enumerate(self.nnkp.bvec_crys):
                    ## s[isym] . bvec[ibi] ~= bvec[ bvec_equiv[ibf,isym1] ]  i.e. S.bi = bf
                    ibi = bvec_equiv[ibf,isym1]
                    bi = self.nnkp.bvec_crys[ibi]
                    ikbi = self.sym.search_ik_full(ki + bi)
                    ikbf = self.sym.search_ik_full(kf + bf)
                    isym2 = self.sym.equiv_sym[ikbi]
                    isym3 = self.sym.equiv_sym[ikbf]
                    # g[isym2]^-1 g[isym1]^-1 g[isym3] = s[isym,:,:], ft[isym,:,:] + tdiff
                    isym, factor, tdiff = self.sym.search_symop( [[isym2, -1], [isym1, -1], [isym3, 1]] )

                    ikbi_eq = self.sym.equiv[ ikbi ]

                    if np.sum(np.abs(self.sym.repmat[ikbi_eq,isym,:,:])) < 1e-5:
                        print(ikbi)
                        print(ikbi_eq)

                    self.kb2k[ikf,ibf] = ikbf
                    if self.sym.t_rev[isym2] == 1:
                        mmn[ikf,ibf,:,:] = np.einsum("mn,np->mp", self.mmn[iki,ibi,:,:], np.conj(self.sym.repmat[ikbi_eq,isym,:,:]), optimize=True) * factor
                    else:
                        mmn[ikf,ibf,:,:] = np.einsum("mn,np->mp", self.mmn[iki,ibi,:,:], self.sym.repmat[ikbi_eq,isym,:,:], optimize=True) * factor
                    if self.sym.t_rev[isym1] == 1:
                        mmn[ikf,ibf,:,:] = np.conj(mmn[ikf,ibf,:,:])

                    # e^{-i b_i T} 
                    #   from g0^-1 e^{-i b_f r} => e^{i b_i r} g0^-1 x e^{-i b_i T}
                    # e^{i (k_i+b_i) Tdiff}
                    #   from h = g0^-1(ki+bi) g0^-1(kf) g0(kf+bf)
                    kbi_eq = self.sym.irr_kpoints[ikbi_eq]
                    arg1 = - np.dot(bi, self.sym.ft[isym1,:])
                    arg2 = - np.dot(kbi_eq, tdiff)
                    if self.sym.t_rev[isym2] == 1:
                        arg2 *= -1
                    if self.sym.t_rev[isym1] == 1:
                        arg1 *= -1
                        arg2 *= -1
                    if self.sym.t_rev[isym] == 1:
                        arg2 *= -1
                    phase = np.exp( 2j*np.pi* (arg1 + arg2) )
                    mmn[ikf,ibf,:,:] *= phase

            self.mmn = mmn

        assert np.all( self.kb2k[:,:] >= 0 )  # check all kb2k is defined

