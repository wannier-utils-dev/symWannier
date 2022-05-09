#!/usr/bin/env python

import numpy as np
import scipy.linalg
import os
import itertools
import gzip

from symwannier.nnkp import Nnkp
from symwannier.sym import Sym

class Amn():
    """
    amn file
    Amn(k) = <psi_mk|g_n>
    """
    def __init__(self, file_amn, nnkp, sym = None):

        self.nnkp = nnkp
        self.sym = sym

        if os.path.exists(file_amn):
            with open(file_amn) as fp:
                self._read_amn(fp)
        elif os.path.exists(file_amn + ".gz"):
            with gzip.open(file_amn + ".gz", 'rt') as fp:
                self._read_amn(fp)
        else:
            raise Exception("failed to read amn file: " + file_amn)

    def _read_amn(self, fp):
        print("  Reading amn file")
        lines = fp.readlines()
        num_bands, nk, num_wann = [ int(x) for x in lines[1].split() ]
        dat = np.genfromtxt(lines[2:]).reshape(nk, num_wann, num_bands, 5)
        amn = np.transpose(dat[:,:,:,3] + 1j*dat[:,:,:,4], axes=(0,2,1))

        self.num_bands = num_bands
        self.num_wann = num_wann

        ####### simple case (without symmetry) #######
        if self.sym is None:
            self.nk = nk
            self.amn = amn

        ####### symmetrized case #######
        else:
            self.nk = self.sym.nkf
            amn = self.symmetrize_Gk(amn)
            self.amn = self.symmetrize_expand(amn)

    def symmetrize_Gk(self, amn, thr=0):
        amn_sym = self._symmetrize_Gk_internal(amn)
        diff = np.sum(np.abs(amn_sym - amn))/self.sym.nks
        print("  symmetrize Gk diff1 = {:15.5e}".format(diff))
        
        if thr > 0:
            for i in range(20):
                amn = amn_sym
                amn_sym = self._symmetrize_Gk_internal(amn)
                diff = np.sum(np.abs(amn_sym - amn))/self.sym.nks
                if diff < thr:
                    break
            print("  symmetrize Gk diff2 = {:15.5e}".format(diff))
        return amn_sym

    def _symmetrize_Gk_internal(self, amn):
        """
        symmetrize Amn/U using G_k
        Amn = <psi_m k | g_n> = 1/N(h) sum_h <psi_m k | h h^-1 | g_n>
            = 1/N(h) sum_h <psi_m k | h | psi_l > < psi_l | h^-1 | g_n>

        | w_n k> = U_mn(k) | psi_m k>
        U_mn(k) = < psi_m k | w_n k >
        """
        amn_sym = np.zeros_like(amn)
        Rmat, Rshift, pos = self.projection_sym_mat()
        for ik, k in enumerate(self.sym.irr_kpoints):
            # for all h (h in G_k)
            nh = 0
            for isym, s in enumerate(self.sym.s):
                sk = np.dot(s, k)
                if self.sym.t_rev[isym] == 1 : sk = -sk
                kdiff = k - sk
                if not np.allclose(kdiff, np.round(kdiff)): continue

                nh += 1
                # Rmat[isym,m,n] = <g_m| S^-1 |g_n>
                phase2 = np.einsum("a,na->n", k, Rshift[isym, :, :], optimize=True)
                phase = np.exp(-1j * 2*np.pi * phase2)
                amn_ik = np.einsum("ml,ln,n->mn", amn[ik,:,:], Rmat[isym,:,:], phase, optimize=True)
                if self.sym.t_rev[isym] == 1:
                    amn_ik = np.conj(amn_ik)
                amn_ik = np.einsum("ml,ln->mn", self.sym.repmat[ik,isym,:,:], amn_ik, optimize=True)
                amn_sym[ik,:,:] += amn_ik
            amn_sym[ik,:,:] /= nh

        return amn_sym

    def symmetrize_expand(self, Umat_irk):
        """
        Generate Amn/Umat(full_k) from Amn/Umat(irr_k)
        Umat(full_k) = Umat(irr_k) (Rmat+phase)

        Amn(k_full) = <psi_m k_f | g_n>
                    = <psi_m k_i | g_0^-1(kf) | g_n>
                    = <psi_m k_i | Rmat(n,l) | g_l(r-R_l)>   ! R_l: lattice vector
                    = exp(i k_i R_l) <psi_m k_i | Rmat(n,l) | g_l(r)>

        irr_k[iks] x s[isym] = full_k[ik]
        """
        _, dim1, dim2 = Umat_irk.shape
        Umat = np.zeros([self.sym.nkf, dim1, dim2], dtype=complex)
        Rmat, Rshift, pos = self.projection_sym_mat()
        for ik, k in enumerate(self.sym.full_kpoints):
            iks = self.sym.equiv[ik]
            ks = self.sym.irr_kpoints[iks,:]
            isym = self.sym.equiv_sym[ik]
            phase2 = np.einsum("a,na->n", ks, Rshift[isym, :, :], optimize=True)
            phase = np.exp(-1j * 2*np.pi * phase2)
            Umat[ik,:,:] = np.einsum("ml,ln,n->mn", Umat_irk[iks,:,:], Rmat[isym,:,:], phase, optimize=True)
            if self.sym.t_rev[isym] == 1:
                Umat[ik,:,:] = np.conj(Umat[ik,:,:])
        return Umat

    def Umat_symmetrize(self, Umat, Umat_opt = None):
        Umat_irk = Umat[self.sym.iks2ik[:]]
        if Umat_opt is not None:
            Umat_opt_irk = Umat_opt[self.sym.iks2ik[:]]
            Umat_full_irk = np.einsum("klm,kmn->kln", Umat_opt_irk, Umat_irk, optimize=True)
            Umat_full_irk = self.symmetrize_Gk(Umat_full_irk)
            Umat_full = self.symmetrize_expand(Umat_full_irk)
            Umat_new = np.einsum("kml,kmn->kln", np.conj(Umat_opt), Umat_full, optimize=True)
        else:
            Umat_irk = self.symmetrize_Gk(Umat_irk)
            Umat_new = self.symmetrize_expand(Umat_irk)
        print("  symmetrize expand diff = {:15.5e}".format(np.sum(np.abs(Umat_new - Umat))/self.nk))
        return Umat_new

    def write_amn(self, file_amn):
        with open(file_amn, "w") as fp:
            fp.write("Amn created by amn.py\n")
            fp.write("{} {} {}\n".format(self.num_bands, self.nk, self.num_wann))
            for ik, m, n in itertools.product( range(self.nk), range(self.num_wann), range(self.num_bands) ):
                # num_bands, num_wann, nk
                fp.write("{0} {1} {2}  {3.real:18.12f}  {3.imag:18.12f}\n".format(n+1, m+1, ik+1, self.amn[ik,n,m]))

    def Umat(self, index_win=None, check=True, Umat_opt=None):
        """
        calculate initial Umat

        amn[:num_bands, :num_wann, :nk]
        m[:num_bands, :num_wann] => SVD => u[:num_bands, :num_bands], v[:num_wann, :num_wann]
        Umat[:num_bands, :num_wann] = u[:num_bands, :num_wann] . v[:num_wann, :num_wann]
        """
        if Umat_opt is None:   # Umat simply from Amn
            Umat = np.zeros_like(self.amn)
            for k in range(self.nk):
                m = self.amn[k,:,:]
                if index_win is not None:
                    m[~index_win[k,:],:] = 0   # m = 0 for outside of index window
                    u, s, v = scipy.linalg.svd(m)
                    Umat[k,index_win[k,:],:] = np.matmul(u[index_win[k,:],:self.num_wann], v)
                else:
                    u, s, v = scipy.linalg.svd(m)
                    Umat[k,:,:] = np.matmul(u[:,:self.num_wann], v)

        else:  # Umat after disentanglement
            Umat = np.zeros([self.nk, self.num_wann, self.num_wann], dtype=complex)
            for k in range(self.nk):
                m = np.einsum("lm,ln->mn", np.conj(Umat_opt[k,:,:]), self.amn[k,:,:], optimize=True)
                u, s, v = scipy.linalg.svd(m)
                Umat[k,:,:] = np.matmul(u, v)
        if self.sym is not None:
            Umat = self.Umat_symmetrize( Umat, Umat_opt )

        if check and self.sym is None:
            # check Umat^\dagger Umat = I
            for k in range(self.nk):
                assert np.allclose( np.matmul(np.transpose(np.conj(Umat[k,:,:])), Umat[k,:,:]), np.eye(self.num_wann), atol=1e-5 ), "Umat is not unitary for k = {}".format(k)

        return Umat

    def projection_sym_mat(self):
        pos = self.nnkp.nw2r     # pos[num_wann, 3]: position of each Wannier
        Rmat = self.sym.rotmat   # Rotation matrix

        # calculate Rshift (Rotated Wannier center R - corresponding Wannier center)
        Rshift = np.zeros([self.sym.nsym, self.num_wann, 3])
        for isym in range(self.sym.nsym):
            #assert np.allclose( np.einsum("ab,cb->ac", Rmat[isym], Rmat[isym]), np.eye(Rmat.shape[1]) )  # Rmat is orthogonal
            for iw in range(self.num_wann):
                rpos = pos[ np.logical_not(Rmat[isym,:,iw] == 0) ]
                r0 = rpos[0]
                for r1 in rpos:
                    assert np.all(r0 == r1)
                Rshift[isym, iw, :] = self.sym.apply_r(isym, pos[iw]) - r0

        assert np.all(np.abs(Rshift - np.round(Rshift)) < 1e-3), "Rshift is not a lattice vector"

        return Rmat, Rshift, pos

