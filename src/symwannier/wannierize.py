#!/usr/bin/env python
"""

 Python version of wannier90
 Maximally localized Wannier functions code

 Usage: wannierize.py [-s] [-S] <prefix>

 Options:
   -s    use symmetry ( prefix_sym.dat is required. )
   -S    use symmetry + site symmetry ( prefix_sym.dat is required. )

 Reference: 
 N. Marzari and D. Vanderbilt PRB 56 12847 (1997):  Referred as R1
 I. Souza, N. Marzari and D. Vanderbilt PRB 65 035109 (2001): R2


 copyright (C) 2021-2022  Takashi Koretsune

"""

import numpy as np
import scipy.linalg
import itertools
import textwrap
import datetime
from docopt import docopt

from symwannier.win import Win
from symwannier.mmn import Mmn
from symwannier.amn import Amn
from symwannier.sym import Sym
from symwannier.eig import Eig
from symwannier.nnkp import Nnkp
from symwannier.timedata import TimeData

class Wannierize:
    def __init__(self, prefix, lsym=False, lsite_sym=False):
        self.time = TimeData()
        self.win = Win(prefix)
        self.nnkp = Nnkp(prefix + ".nnkp")
        self.lsym = lsym
        self.lsite_sym = lsite_sym
        if self.lsite_sym:
            self.lsym = True

        if self.lsym:
            self.time.start_clock("sym read")
            self.sym = Sym(file_sym=prefix + ".isym", nnkp=self.nnkp)
            self.time.stop_clock("sym read")
            ext = 'i'
        else:
            self.sym = None
            ext = ''

        self.time.start_clock("mmn read")
        mmn = Mmn(file_mmn = prefix+"." + ext + "mmn", nnkp=self.nnkp, sym=self.sym)
        self.time.stop_clock("mmn read")

        self.mmn0 = mmn.mmn
        self.num_bands = mmn.num_bands
        self.num_wann = self.nnkp.num_wann
        self.nk = self.nnkp.nk
        self.kpts = self.nnkp.kpoints
        self.nb = self.nnkp.nb
        self.bvec = self.nnkp.bvec
        self.wb = self.nnkp.wb
        self.kb2k = mmn.kb2k   # index of k+b

        self.time.start_clock("amn read")
        self.amn = Amn(prefix+"." + ext + "amn", nnkp=self.nnkp, sym=self.sym)
        self.time.stop_clock("amn read")
        self.eig = Eig(prefix+"." + ext + "eig", sym=self.sym)
        self.Umat_opt = None

        self.file_hr_dat = prefix + "_py_hr.dat"
        self.file_tb_dat = prefix + "_py_tb.dat"

        self.disentangle = (self.num_bands > self.num_wann)

        self.omega_thr = 1e-8    # threshold for wannierization
        self.womega_thr = 1e-10  # threshold for disentanglement

    def run(self):
        """
        main subroutine
        """
        # disentanglement. minimize omega_I
        if self.disentangle:
            self.dis_window()
            self.dis_project()
            self.dis_extract()

        # minimize omega_tot
        self.wannierization()

        # show final spreads and write hr.dat
        self.post_process()

    def wannierization(self):
        """
        wannierization procedure
        """
        self.init_Umat_and_Mmn()

        omega = self.calc_omega(self.mmn)
        print("! Initial:  omega_tot = {:15.8f}    time: {:12.6f}".format(omega, self.time.get_time()))
        self.show_spreads()
        converged = False
        for i in range(self.win.num_iter):
            self.time.start_clock("wannierization")
            omega_prev = omega
            dw = self.calc_dw()
            omega = self.update_mmn(dw, omega)
            print("! {:5d}-th  omega_tot = {:15.8f}    time: {:12.6f}".format(i+1, omega, self.time.get_time()))
            self.show_spreads()
            self.calc_omega_detail(self.mmn)

            converged = np.abs(omega - omega_prev) < self.omega_thr
            if converged:
                print("  convergence has been achieved")
                break
            self.time.stop_clock("wannierization")
        print("! Final:    omega_tot = {:15.8f}    time: {:12.6f}".format(omega, self.time.get_time()))

    def post_process(self):
        if self.disentangle:
            self.Umat = np.einsum("kmn,knp->kmp", self.Umat_opt, self.Umat, optimize=True)

        if self.lsite_sym:
            Umat_irk = self.Umat[self.sym.iks2ik[:]]
            Umat_irk = self.amn.symmetrize_Gk(Umat_irk)
            self.Umat = self.amn.symmetrize_expand(Umat_irk)

        self.mmn = self.update_Mmn_by_Umat(self.mmn0, self.Umat)
        self.calc_omega_detail(self.mmn)
        self.show_spreads()
        self.write_hr()
        self.write_tb()
        self.time.show_all()

    def show_spreads(self):
        print("{:>8s},{:>46s},{:>16s}".format("WF num", "center in cartesian coords", "spread (A^2)"))
        for nw in range(self.num_wann):
            print("  {0:6d}   ( {1[0]:12.6f}, {1[1]:12.6f}, {1[2]:12.6f} )  {2:15.8f}".format(nw+1, self.r[nw].real, self.spreads[nw].real))

    def update_Mmn_by_Umat(self, mmn, Umat):
        """
        return Mmn = U^k Mmn^(0) U^(k+b)
        """
        #for k in range(self.nk):
        #    assert np.allclose( np.matmul(Umat[:,:,k], np.transpose(np.conj(Umat[:,:,k]))), np.eye(self.num_bands) ), "Umat is not unitary for k = {}".format(k)
        return np.einsum("klm, kblp, kbpn->kbmn", np.conj(Umat), mmn, Umat[self.kb2k[:,:],:,:], optimize=True) # Eq. (61)

    def init_Umat_and_Mmn(self):
        """
        initialize self.Umat (U^k) and calculate self.mmn (Mmn = U^k Mmn^(0) U^k+b)
        Mmn0 (input) -> Mmn1 (disentangle) -> Mmn
        """
        if self.disentangle:
            print("  init Umatrix and Mmn using Umatrix_opt")
            self.Umat = self.amn.Umat(index_win=self.index_win, Umat_opt=self.Umat_opt)
            self.mmn1 = self.update_Mmn_by_Umat(self.mmn0, self.Umat_opt)
            self.mmn = self.update_Mmn_by_Umat(self.mmn1, self.Umat)
            self.calc_omega_detail(self.mmn)

        else:  # without disentangle
            print("  init Umatrix and Mmn")
            if hasattr(self, "amn"):
                self.Umat = self.amn.Umat()
                self.mmn1 = self.mmn0
                self.mmn = self.update_Mmn_by_Umat(self.mmn1, self.Umat)

            else:
                self.Umat = np.zeros([self.nk, self.num_bands, self.num_bands], dtype=complex)
                for n in range(self.num_bands):
                    self.Umat[:,n,n] = 1
                self.mmn1 = self.mmn0
                self.mmn = self.mmn1

    def calc_omega(self, mmn):
        """
        calculate Omega = sum <r^2> - <r>^2
        update self.r (<r>)
        """
        #r = 1j/self.nk * np.einsum("b,ba,nnkb->na", self.wb, self.bvec, self.mmn, optimize=True)
        #r -= 1j * np.einsum("b,ba->a", self.wb, self.bvec)
        #print("r = ", r)
        mnn = np.einsum("kbnn->kbn", mmn, optimize=True)
        imlnmnn = np.log(mnn).imag
        r = -1/self.nk * np.einsum("b,ba,kbn->na", self.wb, self.bvec, imlnmnn, optimize=True)

        r2a = np.sum(self.wb)*self.nk - np.einsum("b,kbn,kbn->n", self.wb, mnn, np.conj(mnn), optimize=True)
        r2b = np.einsum("b,kbn->n", self.wb, imlnmnn**2)
        r2 = 1/self.nk * (r2a + r2b).real
        self.r = r
        
        self.spreads = r2 - np.sum(r[:,:].real**2, axis=1)
        return np.sum(self.spreads)

    def calc_omega_detail(self, mmn):
        """
        calculate OmegaI, Omega_D, Omega_OD
        """
        mmn2 = np.einsum("kbmn,kbmn->kb", mmn, np.conj(mmn), optimize=True).real
        OmegaI = np.einsum("b, kb->", self.wb, self.num_wann - mmn2, optimize=True)/self.nk

        mnn2 = np.einsum("kbnn,kbnn->kb", mmn, np.conj(mmn), optimize=True).real
        OmegaOD = np.einsum("b, kb->", self.wb, mmn2 - mnn2, optimize=True)/self.nk

        mnn = np.einsum("kbnn->kbn", mmn, optimize=True)
        imlnmnn = np.log(mnn).imag
        r = -1/self.nk * np.einsum("b,ba,kbn->na", self.wb, self.bvec, imlnmnn, optimize=True)
        qn = imlnmnn + np.einsum("ba, na->bn", self.bvec, r, optimize=True)[np.newaxis,:,:]
        OmegaD = np.einsum("b, kbn->", self.wb, qn**2, optimize=True)/self.nk

        print("  OmegaI  = {:15.10f}".format(OmegaI))
        print("  OmegaD  = {:15.10f}".format(OmegaD))
        print("  OmegaOD = {:15.10f}".format(OmegaOD))

    def calc_dw(self):
        """
        calculate G in R1_Eq.(52)
        """
        mnn = np.einsum("kbnn->kbn", self.mmn, optimize=True)
        imlnmnn = np.log(mnn).imag
        r = -1/self.nk * np.einsum("b,ba,kbn->na", self.wb, self.bvec, imlnmnn, optimize=True)
        Rmn = np.einsum("kbmn,kbn->kbmn", self.mmn, np.conj(mnn), optimize=True)  # R1_Eq.(45)
        qn = imlnmnn + np.einsum("ba, na->bn", self.bvec, r, optimize=True)[np.newaxis,:,:] # R1_Eq.(47)
        T = np.einsum("kbmn, kbn->kbmn", self.mmn, qn/mnn, optimize=True)  # R1_Eq.(48,51)

        AR = (Rmn - np.conj(np.transpose(Rmn, axes=(0,1,3,2))))/2  # A[R] = (R-R^dagger)/2
        ST = (T + np.conj(np.transpose(T, axes=(0,1,3,2))))/(2j)  # S[T] = (T+T^dagger)/2i
        dw = 4/self.nk * np.einsum("b,kbmn->kmn", self.wb, AR - ST, optimize=True)  # R1_Eq.(52)
        return dw

    def calc_expdw(self, dw):
        """
        calculate exp(dW) by diagonalizing i*dW  (i*dW is hermite)
        """
        expdw = np.zeros_like(dw)
        for k in range(self.nk):
            e, v = scipy.linalg.eigh(1j*dw[k,:,:])
            #idw = np.einsum("ab,b,cb->ac", v, e, np.conj(v), optimize=True)
            #assert np.allclose(idw, 1j*dw[:,:,k], atol=1e-5), np.sum(np.abs(idw - 1j*dw[:,:,k]))/np.sum(np.abs(idw))
            expdw[k,:,:] = np.einsum("ab,b,cb->ac", v, np.exp(-1j*e), np.conj(v), optimize=True)
        return expdw

    def update_mmn(self, dw, omega):
        """
        calculate optimal alpha and update self.Umat (U^(k)) and self.mmn (Mmn)
        return omega
        """
        Umat1 = np.einsum("kmn, knl->kml", self.Umat, self.calc_expdw(dw), optimize=True)  # R1_Eq.(60)
        mmn1 = self.update_Mmn_by_Umat(self.mmn1, Umat1)
        omega1 = self.calc_omega(mmn1)

        Umat2 = np.einsum("kmn, knl->kml", self.Umat, self.calc_expdw(dw/2), optimize=True)  # R1_Eq.(60)
        mmn2 = self.update_Mmn_by_Umat(self.mmn1, Umat2)
        omega2 = self.calc_omega(mmn2)

        # get optimum alpha from omega(alpha=0), omega(alpha=1/2)=omega2, omega(alpha=1)=omega1
        if 2*omega + 2*omega1 - 4*omega2 < 0:
            alpha = 0.01 if omega < omega1 else 1.0
        else:
            alpha = -(4*omega2-omega1-3*omega) /2/ (2*omega+2*omega1-4*omega2)
            alpha = min(1, alpha)
            alpha = max(0.01, alpha)
        #print("{:10.5f} {:10.5f} {:10.5f} {:10.5f}".format(omega, omega2, omega1, alpha))

        self.Umat = np.einsum("kmn, knl->kml", self.Umat, self.calc_expdw(alpha*dw), optimize=True)  # R1_Eq.(60)
        #if self.lsite_sym:
        #    self.Umat = self.amn.Umat_symmetrize(self.Umat, self.Umat_opt)
        self.mmn = self.update_Mmn_by_Umat(self.mmn1, self.Umat)
        omega = self.calc_omega(self.mmn)
        return omega

    def dis_window(self):
        # array with self.nband size
        self.index_froz = (self.eig.eig > self.win.dis_froz_min) & (self.eig.eig < self.win.dis_froz_max)
        self.index_win = (self.eig.eig > self.win.dis_win_min) & (self.eig.eig < self.win.dis_win_max)
        self.index_nfroz = self.index_win & (~self.index_froz)
        self.len_nfroz = np.array([ np.sum(self.index_nfroz[k,:]) for k in range(self.nk) ])
        self.ndimwin = np.array([ np.sum(self.index_win[k,:]) for k in range(self.nk) ])
        self.ndimfroz = np.array([ np.sum(self.index_froz[k,:]) for k in range(self.nk) ])

    def dis_project(self):
        """
        calculate self.Umat_opt and update mmn0U from amn file
        """
        assert hasattr(self, "amn"), "file_amn is not specified"
        self.Umat_opt = self.amn.Umat(index_win = self.index_win)  # Umat_opt[:num_bands, :num_wann, :nk]

        # dis_proj_froz in wannier90
        # wannier90 uses ortho-fix which is not implemented here
        for k in range(self.nk):
            # cqpq = Q_inner P_G Q_inner in R2_Eq.(27)
            Umatk = self.Umat_opt[k,self.index_win[k,:],:]
            lnfroz = (~self.index_froz[ k, self.index_win[k,:] ]).astype(int)
            cqpq = np.einsum("m,ml,nl,n->mn", lnfroz, Umatk, np.conj(Umatk), lnfroz, optimize=True)
            # check hermiticity
            #assert np.sum(np.abs(cqpq - np.transpose(np.conj(cqpq)))) < 1e-5  # check hermiticity

            e, v = scipy.linalg.eigh(-cqpq)  # "-" sign for descending order
            self.Umat_opt[k,:,:] = 0
            self.Umat_opt[k, self.index_froz[k,:], :self.ndimfroz[k]] = np.identity(self.ndimfroz[k])
            self.Umat_opt[k, self.index_win[k,:], self.ndimfroz[k]:self.num_wann] = v[:, :self.num_wann - self.ndimfroz[k]]

            # check unitary
            #assert np.allclose( np.matmul(np.transpose(np.conj(self.Umat_opt[:,:,k])), self.Umat_opt[:,:,k]), np.eye(self.num_wann) )

        self.update_mmn0U()

    def update_mmn0U(self):
        """
        update mmn0U = mmn0[num_bands, num_win] x Umat_opt[num_win, num_wann]
        """
        self.mmn0U = np.einsum("kbml, kbln->kbmn", self.mmn0, self.Umat_opt[self.kb2k[:,:],:,:], optimize=True) # Eq. (61)

    def init_zmatrix(self):
        """
        calculate Zmatrix (zmat)
        Zmatrix = sum_b wb <u^(0)_mk | P^(i)_k+b | u^0_nk>
                = sum_{bl} wb <u^(0)_mk | u^(i)_lk+b > < u^(i)_lk+b | u^0_nk >
        mu : Mmn_k U^(i)_k+b = <u^(0)_mk | u^(i)_lk+b >
        """
        # mu = mmn0[num_bands, num_win] x Umat_opt[num_win, num_wann]
        # zmat = mu x mu^\dagger
        zmat_froz = np.zeros([self.nk, self.num_bands, self.num_bands], dtype=complex)
        zmat_nfroz = np.zeros([self.nk, self.num_bands, self.num_bands], dtype=complex)
        #for k in range(self.nk):
        #    mmn0U_froz = np.compress(self.index_froz[k,:], self.mmn0U[k,:,:,:], axis=1)
        #    zmat_froz[k, :self.ndimfroz[k], :self.ndimfroz[k]] = np.einsum("b, bml, bnl->mn", self.wb, mmn0U_froz, np.conj(mmn0U_froz), optimize=True)
        #    mmn0U_nfroz = np.compress(self.index_nfroz[k,:], self.mmn0U[k,:,:,:], axis=1)
        #    zmat_nfroz[k, :self.len_nfroz[k], :self.len_nfroz[k]] = np.einsum("b, bml, bnl->mn", self.wb, mmn0U_nfroz, np.conj(mmn0U_nfroz), optimize=True)
        zmat = np.einsum("b, kbml, kbnl->kmn", self.wb, self.mmn0U, np.conj(self.mmn0U), optimize=True)
        for k in range(self.nk):
            zmat_froz[k,:self.ndimfroz[k],:self.ndimfroz[k]] = zmat[k, self.index_froz[k,:], :][:, self.index_froz[k,:]]
            zmat_nfroz[k,:self.len_nfroz[k],:self.len_nfroz[k]] = zmat[k, self.index_nfroz[k,:], :][:, self.index_nfroz[k,:]]

        return zmat_froz, zmat_nfroz

    def calc_womegaI(self):
        mmn = np.einsum("kml,kbmn->kbln", np.conj(self.Umat_opt), self.mmn0U, optimize=True)
        return self.num_wann * np.sum(self.wb) - np.einsum("b,kbmn,kbmn->", self.wb, mmn, np.conj(mmn), optimize=True).real/self.nk

    def dis_extract(self):
        for i in range(self.win.dis_num_iter):
            zmat_froz, zmat_nfroz = self.init_zmatrix()
            if not i == 0:
                # update zmatrix
                zmat_nfroz = self.win.dis_mix_ratio * zmat_nfroz + (1-self.win.dis_mix_ratio) * zmat_nfroz_old

            wkomegaI1 = self.num_wann * np.sum(self.wb) * np.ones([self.nk])

            for k in range(self.nk):
                wkomegaI1[k] -= np.trace(zmat_froz[k, :self.ndimfroz[k], :self.ndimfroz[k]]).real

            for k in range(self.nk):
                # Here, we consider -Z because we need eigvals of Z in descending order
                e, v = scipy.linalg.eigh(-zmat_nfroz[k, :self.len_nfroz[k], :self.len_nfroz[k]])
                e = -e
                wkomegaI1[k] -= np.sum( e[:(self.num_wann-self.ndimfroz[k])] )

                self.Umat_opt[k, self.index_win[k,:], self.ndimfroz[k]:self.num_wann] = 0
                self.Umat_opt[k, self.index_nfroz[k,:], self.ndimfroz[k]:self.num_wann] = v[:,:self.num_wann-self.ndimfroz[k]]

                # check unitary
                #assert np.allclose( np.matmul(np.transpose(np.conj(self.Umat_opt[:,:,k])), self.Umat_opt[:,:,k]), np.eye(self.num_wann) )
            self.update_mmn0U()
            
            womegaI1 = np.sum(wkomegaI1)/self.nk

            womegaI = self.calc_womegaI()

            print("{:5d} {:15.8f} {:15.8f} {:15.5e}  {:12.6f}".format(i, womegaI1, womegaI, womegaI1/womegaI - 1, self.time.get_time()))
            if np.abs(womegaI1/womegaI - 1) < self.womega_thr:
                break

            zmat_nfroz_old = zmat_nfroz

    def calc_irvec(self):
        real_metric = np.einsum("ij,lj->il", self.nnkp.a, self.nnkp.a)
        mp_grid = self.win.mp_grid
        ilist = np.array(list(itertools.product( range(-2,3), repeat=3 )))
        mp_list = np.array(list(
                  itertools.product( range(-mp_grid[0], mp_grid[0]+1),
                                     range(-mp_grid[1], mp_grid[1]+1),
                                     range(-mp_grid[2], mp_grid[2]+1) )))
        # distance from Wannier on superlattice (LxMxN). size R_list[npos,LMN,3], dist[npos,LMN]
        R_list = mp_list[:,None,:] - np.einsum("na,a->na", ilist, mp_grid, optimize=True)[None,:,:]
        dist = np.einsum("ija,ab,ijb->ij", R_list, real_metric, R_list, optimize=True)
        # dist[:,62] = dist[:, i1=i2=i3=0]
        dist0 = dist[:,62]
        dist_min = np.min(dist, axis=1)
        ind = np.where(dist0 - dist_min < 1e-5)
        ndeg = np.sum( np.abs(dist - dist0[:,None]) < 1e-5, axis=1 )
        self.irvec = mp_list[ind]
        self.ndegen = ndeg[ind]
        # check sum rule
        assert np.sum(1/self.ndegen) - np.prod(mp_grid) < 1e-8, "error in finding Wigner-Seitz points"

    def write_hr(self):
        self.calc_irvec()

        ham_k = np.einsum("kni,kn,knj->kij", np.conj(self.Umat), self.eig.eig, self.Umat, optimize=True)
        kr = np.einsum("ka,ra->kr", self.kpts, self.irvec, optimize=True)
        ham_r = np.einsum("kij,kr->ijr", ham_k, np.exp(-2 * np.pi * 1j * kr), optimize=True)/self.nk

        t = datetime.datetime.now()
        with open(self.file_hr_dat, "w") as fp:
            fp.write(" written {}\n".format(t.strftime('on %d%b%Y at %H:%M:%S')))
            fp.write("{:12d}\n{:12d}\n".format(self.num_wann, len(self.ndegen)))
            fp.write(textwrap.fill("".join(["{:5d}".format(x) for x in self.ndegen]), 75, drop_whitespace=False))
            fp.write("\n")
            for irpts in range(len(self.ndegen)):
                for i, j in itertools.product( range(self.num_wann), repeat=2 ):
                    line = "".join( "{:5d}".format(x) for x in self.irvec[irpts,:] )
                    line += "{:5d}{:5d}".format(j+1, i+1)
                    line += "{:12.6f}{:12.6f}".format(ham_r[j,i,irpts].real, ham_r[j,i,irpts].imag)
                    #line += "{:15.10f}{:15.10f}".format(ham_r[i,j,irpts].real, ham_r[i,j,irpts].imag)
                    fp.write(line + "\n")

    def write_tb(self):
        self.calc_irvec()

        ham_k = np.einsum("kni,kn,knj->kij", np.conj(self.Umat), self.eig.eig, self.Umat, optimize=True)
        kr = np.einsum("ka,ra->kr", self.kpts, self.irvec, optimize=True)
        fac = np.exp(-2 * np.pi * 1j * kr)
        ham_r = np.einsum("kij,kr->ijr", ham_k, fac, optimize=True)/self.nk

        mnn = np.einsum("kbnn->kbn", self.mmn, optimize=True)
        imlnmnn = np.log(mnn).imag
        rk = - np.einsum("b,ba,kbn->ank", self.wb, self.bvec, imlnmnn, optimize=True)
        rr = np.einsum("ank,kr->anr", rk, fac, optimize=True)/self.nk
        rr2 = 1j * np.einsum("b,ba,kbmn,kr->amnr", self.wb, self.bvec, self.mmn, fac, optimize=True)/self.nk

        t = datetime.datetime.now()
        with open(self.file_tb_dat, "w") as fp:
            fp.write(" written {}\n".format(t.strftime('on %d%b%Y at %H:%M:%S')))
            fp.write(" {0[0]:15.8f} {0[1]:15.8f} {0[2]:15.8f}\n".format(self.nnkp.a[0,:]))
            fp.write(" {0[0]:15.8f} {0[1]:15.8f} {0[2]:15.8f}\n".format(self.nnkp.a[1,:]))
            fp.write(" {0[0]:15.8f} {0[1]:15.8f} {0[2]:15.8f}\n".format(self.nnkp.a[2,:]))
            fp.write("{:12d}\n{:12d}\n".format(self.num_wann, len(self.ndegen)))
            fp.write(textwrap.fill("".join(["{:5d}".format(x) for x in self.ndegen]), 75, drop_whitespace=False))
            fp.write("\n")
            #
            # write <0n|H|Rm>
            #
            for irpts in range(len(self.ndegen)):
                fp.write("\n")
                fp.write("".join( "{:5d}".format(x) for x in self.irvec[irpts,:] ))
                fp.write("\n")
                for i, j in itertools.product( range(self.num_wann), repeat=2 ):
                    line = "{:5d}{:5d}  ".format(j+1, i+1)
                    line += " {:15.8e} {:15.8e}".format(ham_r[j,i,irpts].real, ham_r[j,i,irpts].imag)
                    fp.write(line + "\n")
            #
            # write <0n|r|Rm>
            #
            for irpts in range(len(self.ndegen)):
                fp.write("\n")
                fp.write("".join( "{:5d}".format(x) for x in self.irvec[irpts,:] ))
                fp.write("\n")
                for i, j in itertools.product( range(self.num_wann), repeat=2 ):
                    line = "{:5d}{:5d}  ".format(j+1, i+1)
                    if i == j:
                        line += "".join([" {:15.8e} {:15.8e}".format(x.real, x.imag) for x in rr[:,i,irpts]])
                    else:
                        line += "".join([" {:15.8e} {:15.8e}".format(x.real, x.imag) for x in rr2[:,j,i,irpts]])
                    fp.write(line + "\n")


if __name__ == '__main__':
    args = docopt(__doc__)

    prefix = args["<prefix>"]

    wann = Wannierize(prefix=prefix, lsym=args["-s"], lsite_sym=args["-S"])
    wann.run()
