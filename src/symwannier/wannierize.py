#!/usr/bin/env python
"""
 Compute maximally localized Wannier functions using symmetry information.

 Reference: 
   R1. N. Marzari and D. Vanderbilt PRB 56 12847 (1997)
   R2. I. Souza, N. Marzari and D. Vanderbilt PRB 65 035109 (2001)
   R3. T. Koretsune Comp. Phys. Comm. 285 108645 (2023)
"""

import numpy as np
import scipy.linalg
import scipy.optimize
import itertools
import textwrap
import datetime
import argparse
import sys
import logging

from symwannier.win import Win
from symwannier.mmn import Mmn
from symwannier.amn import Amn
from symwannier.sym import Sym
from symwannier.eig import Eig
from symwannier.nnkp import Nnkp
from symwannier.timedata import TimeData

class Wannierize:
    """Compute maximally localized (and symmetry-adapted) Wannier functions.

    Parameters
    ----------
    prefix : str
        Prefix for input/output files (e.g., ``prefix.win``, ``prefix.mmn``).
    lsym : bool, optional
        Use symmetry-aware inputs (``*.isym``, ``*.immn``, ``*.iamn``, ``*.ieig``).
    lsite_sym : bool, optional
        Apply site symmetry after Wannierization (implies ``lsym=True``).
    prec : bool, optional
        If True, write high-precision outputs for hr.dat / tb.dat.
    optimize_memory_usage : bool, optional
        If True, prioritize memory efficiency over speed (default: False).
        If False (default), prioritize speed (may use more memory). 
    projectability_disentangle : bool, optional
        If True, build disentanglement windows from AMN projectability instead of energy windows.
    log_level : int or str, optional
        Logging level (e.g., logging.DEBUG, logging.INFO, or 'DEBUG', 'INFO').
    log : logging.Logger, optional
        Logger to use; if not provided a module logger is created.
    """

    def __init__(self, prefix, lsym=False, lsite_sym=False, prec=False, optimize_memory_usage=False, projectability_disentangle=False, log_level=None, log=None, snap_kp=True):
        self.log = log or logging.getLogger(__name__)
        if not self.log.handlers:
            # Determine log level
            if log_level is None:
                level = logging.INFO
            elif isinstance(log_level, str):
                level = getattr(logging, log_level.upper(), logging.INFO)
            else:
                level = log_level
            logging.basicConfig(level=level, format="%(message)s")
        self.log.debug(f"Initializing Wannierize with prefix={prefix}, lsym={lsym}, lsite_sym={lsite_sym}")
        self.prefix = prefix
        self.time = TimeData(log=self.log)
        self.win = Win(prefix, log=self.log)
        self.log.debug(f"Loaded win: num_wann={self.win.num_wann}")
        self.nnkp = Nnkp(prefix + ".nnkp", log=self.log, snap_kp=snap_kp)
        self.log.debug(f"Loaded nnkp: nk={self.nnkp.nk}, nb={self.nnkp.nb}, num_wann={self.nnkp.num_wann}")
        self.lsym = lsym
        self.lsite_sym = lsite_sym
        self.prec = prec
        self.optimize_memory_usage = optimize_memory_usage
        self.projectability_disentangle = projectability_disentangle
        if self.lsite_sym:
            self.lsym = True

        if self.lsym:
            self.time.start_clock("sym read")
            self.sym = Sym(file_sym=prefix + ".isym", nnkp=self.nnkp, log=self.log)
            self.time.stop_clock("sym read")
            self.log.debug(f"Loaded symmetry: nsym={len(self.sym.s)}")
            ext = 'i'
        else:
            self.sym = None
            ext = ''

        self.time.start_clock("mmn read")
        mmn = Mmn(file_mmn = prefix+"." + ext + "mmn", nnkp=self.nnkp, sym=self.sym, log=self.log)
        self.time.stop_clock("mmn read")
        self.log.debug(f"Loaded Mmn: shape={mmn.mmn.shape}")

        self.mmn0 = mmn.mmn
        self.num_bands = mmn.num_bands
        self.num_wann = self.win.num_wann
        self.nk = self.nnkp.nk
        self.kpts = self.nnkp.kpoints
        self.nb = self.nnkp.nb
        self.bvec = self.nnkp.bvec
        self.wb = self.nnkp.wb
        self.kb2k = mmn.kb2k   # index of k+b

        self.time.start_clock("amn read")
        self.amn = Amn(prefix+"." + ext + "amn", nnkp=self.nnkp, sym=self.sym, log=self.log)
        self.time.stop_clock("amn read")
        self.eig = Eig(prefix+"." + ext + "eig", sym=self.sym, log=self.log)
        self.log.debug(f"Loaded Amn: shape={self.amn.amn.shape}, Eig: shape={self.eig.eig.shape}")
        # Projectability p_mk = sum_n |<psi_mk|g_n>|^2 for optional disentanglement mode
        self.projectability = np.einsum(
            "knm,knm->kn", self.amn.amn, np.conj(self.amn.amn), optimize=True
        ).real
        is_finite = np.isfinite(self.projectability)
        num_nan = np.isnan(self.projectability).sum()
        num_inf = np.isinf(self.projectability).sum()
        num_neg = (self.projectability < 0).sum()
        self.log.info(
            "projectability stats: "
            f"min={np.min(self.projectability):.6e}, max={np.max(self.projectability):.6e}, "
            f"nan={num_nan}, inf={num_inf}, neg={num_neg}"
        )
        if not is_finite.all():
            raise ValueError("projectability contains NaN or inf")
        if num_neg != 0:
            raise ValueError("projectability contains negative values")
        self._validate_inputs()
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
        self.log.info("Starting wannierization for prefix '%s' (lsym=%s, lsite_sym=%s)", self.prefix, self.lsym, self.lsite_sym)
        # disentanglement. minimize omega_I
        if self.disentangle:
            if self.projectability_disentangle:
                self.dis_window_projectability()
            else:
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
        self.log.debug("Starting wannierization loop")
        self.init_Umat_and_Mmn()

        omega = self.calc_omega(self.mmn)
        self.log.debug(f"Initial omega_tot = {omega:.8f}")
        print("! Initial:  omega_tot = {:15.8f}    time: {:12.6f}".format(omega, self.time.get_time()))
        self.show_spreads()
        converged = False
        self.time.start_clock("wannierization")
        for i in range(self.win.num_iter):
            self.log.debug(f"Wannierization iteration {i+1}/{self.win.num_iter}")
            omega_prev = omega
            dw = self.calc_dw()
            omega = self.update_mmn(dw, omega)
            self.log.debug(f"Iteration {i+1}: omega_tot = {omega:.8f}, delta_omega = {abs(omega - omega_prev):.2e}")
            print("! {:5d}-th  omega_tot = {:15.8f}    time: {:12.6f}".format(i+1, omega, self.time.get_time()))
            self.show_spreads()
            self.calc_omega_detail(self.mmn)

            converged = np.abs(omega - omega_prev) < self.omega_thr
            if converged:
                self.log.debug("Wannierization converged")
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

    def _validate_inputs(self):
        """Validate basic consistency of loaded inputs to fail fast with clear errors."""
        self.log.debug(f"Validating inputs: num_wann={self.num_wann}, num_bands={self.num_bands}, nk={self.nnkp.nk}, nb={self.nnkp.nb}")
        if self.num_wann > self.num_bands:
            raise ValueError("num_wann ({}) must not exceed num_bands ({})".format(self.num_wann, self.num_bands))

        if self.mmn0.shape[0] != self.nnkp.nk or self.mmn0.shape[1] != self.nnkp.nb:
            raise ValueError("mmn shape {} inconsistent with nk={} or nb={}".format(self.mmn0.shape, self.nnkp.nk, self.nnkp.nb))

        if getattr(self.amn, "nk", None) and self.amn.nk != self.nnkp.nk:
            raise ValueError("amn nk={} inconsistent with nnkp nk={}".format(self.amn.nk, self.nnkp.nk))

        if getattr(self.eig, "nk", None) and self.eig.nk != self.nnkp.nk:
            raise ValueError("eig nk={} inconsistent with nnkp nk={}".format(self.eig.nk, self.nnkp.nk))

    def show_spreads(self):
        """Print Wannier function centers and spreads in Cartesian coordinates."""
        print("{:>8s},{:>46s},{:>16s}".format("WF num", "center in cartesian coords", "spread (A^2)"))
        for nw in range(self.num_wann):
            print("  {0:6d}   ( {1[0]:12.6f}, {1[1]:12.6f}, {1[2]:12.6f} )  {2:15.8f}".format(nw+1, self.r[nw].real, self.spreads[nw].real))

    def update_Mmn_by_Umat(self, mmn, Umat):
        """Transform overlap matrix Mmn by unitary rotation.
        
        Computes Mmn = U^k Mmn^(0) U^(k+b) following R1_Eq.(61).
        Uses either global pre-compute (fast, memory-heavy) or batched k-loop (balanced).
        
        Parameters
        ----------
        mmn : ndarray, shape (nk, nb, num_bands, num_bands)
            Overlap matrices between neighboring k-points.
        Umat : ndarray, shape (nk, num_bands, num_wann)
            Unitary rotation matrices at each k-point.
        
        Returns
        -------
        ndarray, shape (nk, nb, num_wann, num_wann)
            Rotated overlap matrices.
        """
        if not self.optimize_memory_usage:
            # Pre-compute indexed Umat for all k+b pairs to avoid repeated indexing
            Umat_kpb = Umat[self.kb2k.ravel()].reshape(self.nk, self.nb, Umat.shape[1], Umat.shape[2])
            return np.einsum("klm, kblp, kbpn->kbmn", np.conj(Umat), mmn, Umat_kpb, optimize=True)
        else:
            # Memory-efficient: re-index on each call (no pre-computation)
            return np.einsum("klm, kblp, kbpn->kbmn", np.conj(Umat), mmn, Umat[self.kb2k[:,:],:,:], optimize=True)

    def init_Umat_and_Mmn(self):
        """
        initialize self.Umat (U^k) and calculate self.mmn (Mmn = U^k Mmn^(0) U^k+b)
        Mmn0 (input) -> Mmn1 (disentangle) -> Mmn
        """
        self.log.debug(f"init_Umat_and_Mmn: disentangle={self.disentangle}")
        if self.disentangle:
            self.log.info("init Umatrix and Mmn using Umatrix_opt")
            self.Umat = self.amn.Umat(index_win=self.index_win, Umat_opt=self.Umat_opt)
            self.mmn1 = self.update_Mmn_by_Umat(self.mmn0, self.Umat_opt)
            self.mmn = self.update_Mmn_by_Umat(self.mmn1, self.Umat)
            self.calc_omega_detail(self.mmn)

        else:  # without disentangle
            self.log.info("init Umatrix and Mmn")
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
        """Calculate total spread functional Omega.
        
        Computes Omega = sum_n [<r^2>_n - <r>_n^2] following R1_Eqs.(11,31,32).
        Updates self.r (Wannier function centers) and self.spreads as side effects.
        Caches mnn and imlnmnn for reuse in calc_omega_detail.
        
        Parameters
        ----------
        mmn : ndarray, shape (nk, nb, num_wann, num_wann)
            Overlap matrices M^(k,b)_mn.
        
        Returns
        -------
        float
            Total spread Omega (sum of individual spreads).
        """
        self.log.debug(f"calc_omega: mmn.shape={mmn.shape}")
        # Cache for reuse in calc_omega_detail and calc_dw
        self._cached_mnn = np.einsum("kbnn->kbn", mmn, optimize=True)
        self._assert_nonzero_mnn(self._cached_mnn, "calc_omega")
        self._cached_imlnmnn = np.log(self._cached_mnn).imag
        self._cached_r = -1/self.nk * np.einsum("b,ba,kbn->na", self.wb, self.bvec, self._cached_imlnmnn, optimize=True)
        
        r2a = np.sum(self.wb)*self.nk - np.einsum("b,kbn,kbn->n", self.wb, self._cached_mnn, np.conj(self._cached_mnn), optimize=True)
        r2b = np.einsum("b,kbn->n", self.wb, self._cached_imlnmnn**2)
        r2 = 1/self.nk * (r2a + r2b).real
        self.r = self._cached_r
        
        self.spreads = r2 - np.sum(self.r[:,:].real**2, axis=1)
        return np.sum(self.spreads)

    def calc_omega_detail(self, mmn):
        """Decompose spread into gauge-invariant (I), diagonal (D), and off-diagonal (OD) parts.
        
        Following R1_Eq.(13,18), Omega = OmegaI + OmegaD + OmegaOD.
        Prints the three components to stdout.
        Reuses cached mnn and imlnmnn from calc_omega if available.
        
        Parameters
        ----------
        mmn : ndarray, shape (nk, nb, num_wann, num_wann)
            Overlap matrices M^(k,b)_mn.
        """
        mmn2 = np.einsum("kbmn,kbmn->kb", mmn, np.conj(mmn), optimize=True).real
        OmegaI = np.einsum("b, kb->", self.wb, self.num_wann - mmn2, optimize=True)/self.nk

        mnn2 = np.einsum("kbnn,kbnn->kb", mmn, np.conj(mmn), optimize=True).real
        OmegaOD = np.einsum("b, kb->", self.wb, mmn2 - mnn2, optimize=True)/self.nk

        # Reuse cached values if available (set by calc_omega)
        if hasattr(self, '_cached_mnn') and hasattr(self, '_cached_imlnmnn'):
            mnn = self._cached_mnn
            imlnmnn = self._cached_imlnmnn
        else:
            mnn = np.einsum("kbnn->kbn", mmn, optimize=True)
            self._assert_nonzero_mnn(mnn, "calc_omega_detail")
            imlnmnn = np.log(mnn).imag
        
        r = -1/self.nk * np.einsum("b,ba,kbn->na", self.wb, self.bvec, imlnmnn, optimize=True)
        qn = imlnmnn + np.einsum("ba, na->bn", self.bvec, r, optimize=True)[np.newaxis,:,:]
        OmegaD = np.einsum("b, kbn->", self.wb, qn**2, optimize=True)/self.nk

        print("  OmegaI  = {:15.10f}".format(OmegaI))
        print("  OmegaD  = {:15.10f}".format(OmegaD))
        print("  OmegaOD = {:15.10f}".format(OmegaOD))

    def calc_dw(self):
        """Compute the gradient matrix G for steepest-descent update.
        
        Calculates the anti-Hermitian matrix dW(k) = G(k) following R1_Eq.(52),
        used to update the unitary rotations in the direction of decreasing spread.
        Reuses cached mnn, imlnmnn, r from calc_omega if available.
        
        Returns
        -------
        ndarray, shape (nk, num_wann, num_wann)
            Gradient matrices dW(k), anti-Hermitian at each k-point.
        """
        self.log.debug("calc_dw: computing gradient matrices")
        # Reuse cached values if available (set by calc_omega)
        if hasattr(self, '_cached_mnn') and hasattr(self, '_cached_imlnmnn') and hasattr(self, '_cached_r'):
            mnn = self._cached_mnn
            imlnmnn = self._cached_imlnmnn
            r = self._cached_r
        else:
            mnn = np.einsum("kbnn->kbn", self.mmn, optimize=True)
            self._assert_nonzero_mnn(mnn, "calc_dw")
            imlnmnn = np.log(mnn).imag
            r = -1/self.nk * np.einsum("b,ba,kbn->na", self.wb, self.bvec, imlnmnn, optimize=True)
        
        Rmn = np.einsum("kbmn,kbn->kbmn", self.mmn, np.conj(mnn), optimize=True)  # R1_Eq.(45)
        qn = imlnmnn + np.einsum("ba, na->bn", self.bvec, r, optimize=True)[np.newaxis,:,:] # R1_Eq.(47)
        T = np.einsum("kbmn, kbn->kbmn", self.mmn, qn/mnn, optimize=True)  # R1_Eq.(48,51)

        AR = (Rmn - np.conj(np.transpose(Rmn, axes=(0,1,3,2))))/2  # A[R] = (R-R^dagger)/2
        ST = (T + np.conj(np.transpose(T, axes=(0,1,3,2))))/(2j)  # S[T] = (T+T^dagger)/2i
        dw = 4/self.nk * np.einsum("b,kbmn->kmn", self.wb, AR - ST, optimize=True)  # R1_Eq.(52)
        return dw

    def _assert_nonzero_mnn(self, mnn, context):
        """Validate that diagonal overlaps are nonzero before logarithm.
        
        Parameters
        ----------
        mnn : ndarray
            Diagonal overlap elements M^(k,b)_nn.
        context : str
            Caller name for error messages.
        
        Raises
        ------
        ValueError
            If any element of mnn is exactly zero.
        """
        self.log.debug(f"_assert_nonzero_mnn: checking {context}, mnn.shape={mnn.shape}")
        if np.any(mnn == 0):
            count = np.count_nonzero(mnn == 0)
            self.log.error(f"{context}: mnn contains {count} zero elements")
            raise ValueError("{}: mnn contains {} zero elements; log undefined".format(context, count))

    def calc_expdw(self, dw):
        """Compute matrix exponential exp(dW) via eigendecomposition.
        
        Since dW is anti-Hermitian, i*dW is Hermitian and can be diagonalized with real eigenvalues.
        Then exp(dW) = V exp(-i*lambda) V^dagger.
        
        Parameters
        ----------
        dw : ndarray, shape (nk, num_wann, num_wann)
            Anti-Hermitian gradient matrices.
        
        Returns
        -------
        ndarray, shape (nk, num_wann, num_wann)
            Unitary matrices exp(dW(k)).
        """
        # Vectorized eigendecomposition for all k-points
        e = np.empty((self.nk, self.num_wann), dtype=float)
        v = np.empty((self.nk, self.num_wann, self.num_wann), dtype=complex)
        for k in range(self.nk):
            e[k], v[k] = scipy.linalg.eigh(1j*dw[k,:,:])
        # Vectorized matrix multiplication: V @ diag(exp(-i*lambda)) @ V^H
        exp_e = np.exp(-1j * e)  # (nk, num_wann)
        expdw = np.einsum("kab,kb,kcb->kac", v, exp_e, np.conj(v), optimize=True)
        return expdw

    def update_mmn(self, dw, omega):
        """Update unitary rotations and overlaps using line search for optimal step size.
        
        Uses Brent's method for robust 1D minimization to find optimal mixing parameter alpha.
        Evaluates Omega(alpha) and finds the minimum with adaptive sampling.
        
        Parameters
        ----------
        dw : ndarray, shape (nk, num_wann, num_wann)
            Gradient matrices dW(k).
        omega : float
            Current total spread (at alpha=0).
        
        Returns
        -------
        float
            New total spread after optimal rotation.
        """
        self.log.debug("update_mmn: computing line search with Brent's method")
        
        # Define objective function for line search
        def omega_at_alpha(alpha):
            """Compute Omega for given alpha step size."""
            expdw_alpha = self.calc_expdw(alpha * dw)
            Umat_alpha = np.einsum("kmn, knl->kml", self.Umat, expdw_alpha, optimize=True)
            mmn_alpha = self.update_Mmn_by_Umat(self.mmn1, Umat_alpha)
            return self.calc_omega(mmn_alpha)
        
        # Use Brent's method for 1D minimization in range [0, 1]
        # bracket=(0, 1) ensures search within [0, 1]
        result = scipy.optimize.minimize_scalar(
            omega_at_alpha,
            bounds=(0.01, 1.0),
            method='bounded',
            options={'xatol': 1e-4}  # Tolerance for alpha convergence
        )
        
        alpha = result.x
        omega_new = result.fun
        self.log.debug(f"update_mmn: optimal alpha={alpha:.4f}, omega={omega_new:.8f} (nfev={result.nfev})")
        
        # Update with optimal alpha
        expdw_opt = self.calc_expdw(alpha * dw)
        self.Umat = np.einsum("kmn, knl->kml", self.Umat, expdw_opt, optimize=True)
        self.mmn = self.update_Mmn_by_Umat(self.mmn1, self.Umat)
        omega = self.calc_omega(self.mmn)
        
        return omega

    def dis_window(self):
        """Define energy windows for disentanglement."""
        # array with self.nband size
        self.index_froz = (self.eig.eig > self.win.dis_froz_min) & (self.eig.eig < self.win.dis_froz_max)
        self.index_win = (self.eig.eig > self.win.dis_win_min) & (self.eig.eig < self.win.dis_win_max)
        self.index_nfroz = self.index_win & (~self.index_froz)
        self.len_nfroz = np.array([ np.sum(self.index_nfroz[k,:]) for k in range(self.nk) ])
        self.ndimwin = np.array([ np.sum(self.index_win[k,:]) for k in range(self.nk) ])
        self.ndimfroz = np.array([ np.sum(self.index_froz[k,:]) for k in range(self.nk) ])
        self.log.debug(f"dis_window: ndimwin min={np.min(self.ndimwin)}, max={np.max(self.ndimwin)}, ndimfroz min={np.min(self.ndimfroz)}, max={np.max(self.ndimfroz)}")

    def dis_window_projectability(self, dis_proj_max = None, dis_proj_min = None):
        """Define disentanglement windows using projectability from Amn overlaps.

        The projectability of each Bloch state is p_mk = sum_n |<psi_mk|g_n>|^2.
        Bands with high projectability are frozen; the outer window collects the
        best-projected bands until a reasonable buffer over num_wann is reached.
        
        If dis_proj_max is not specified, it is automatically determined by:
        - Sorting all bands by projectability across all k-points
        - Freezing approximately 70-80% of num_wann bands with highest projectability
        - Leaving room for disentanglement to optimize the remaining bands
        """
        proj = self.projectability
        if dis_proj_max is not None and 0 <= dis_proj_max <= 1:
            self.log.info(f"Using user-specified dis_proj_max = {dis_proj_max:.4f}")
        else:
            if hasattr(self.win, "dis_proj_max"):
                dis_proj_max = self.win.dis_proj_max
                self.log.info(f"Using win-specified dis_proj_max = {dis_proj_max:.4f}")
            else:
                # Smart determination: freeze ~75% of num_wann to leave room for disentanglement
                # This balances stability (frozen bands) with flexibility (optimizable bands)
                proj_flat = proj.flatten()
                proj_sorted = np.sort(proj_flat)[::-1]  # descending order
                num_frozen_target = int(0.75 * self.num_wann * self.nk)
                idx_cutoff = min(num_frozen_target, len(proj_sorted) - 1)
                dis_proj_max = proj_sorted[idx_cutoff]
                # Add small margin to avoid floating point issues
                dis_proj_max = max(0.90, dis_proj_max - 0.01)
                self.log.info(f"Auto-determined dis_proj_max = {dis_proj_max:.4f} (targeting ~75% of num_wann as frozen)")
        
        if dis_proj_min is not None and 0 <= dis_proj_min <= 1:
            self.log.info(f"Using user-specified dis_proj_min = {dis_proj_min:.4f}")
        else:
            if hasattr(self.win, "dis_proj_min"):
                dis_proj_min = self.win.dis_proj_min
                self.log.info(f"Using win-specified dis_proj_min = {dis_proj_min:.4f}")
            else:
                # Smart determination: include enough bands for disentanglement
                # Typically want num_wann + buffer bands in the window
                buffer = min(self.num_wann // 2, self.num_bands - self.num_wann)
                proj_flat = proj.flatten()
                proj_sorted = np.sort(proj_flat)[::-1]  # descending order
                idx_cutoff = min(self.num_wann + buffer * self.nk, len(proj_sorted) - 1)
                dis_proj_min = proj_sorted[idx_cutoff]
                # Ensure minimum threshold is reasonable
                dis_proj_min = max(0.05, min(0.25, dis_proj_min - 0.01))
                self.log.info(f"Auto-determined dis_proj_min = {dis_proj_min:.4f} (buffer={buffer} bands per k-point)")

        self.index_froz = dis_proj_max <= self.projectability
        self.index_win = dis_proj_min <= self.projectability
        self.index_nfroz = self.index_win & (~self.index_froz)
        self.len_nfroz = np.array([np.sum(self.index_nfroz[k, :]) for k in range(self.nk)])
        self.ndimwin = np.array([np.sum(self.index_win[k, :]) for k in range(self.nk)])
        self.ndimfroz = np.array([np.sum(self.index_froz[k, :]) for k in range(self.nk)])
        self.log.info(f"num_inner (min, max, mean): {np.min(self.ndimfroz)}, {np.max(self.ndimfroz)}, {np.mean(self.ndimfroz)}")
        self.log.info(f"num_outer (min, max, mean): {np.min(self.ndimwin)}, {np.max(self.ndimwin)}, {np.mean(self.ndimwin)}")

    def dis_project(self):
        """Calculate self.Umat_opt and update mmn0 from amn file."""
        self.log.debug("dis_project: starting projection")
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
                    if self.prec:
                        line += "{:20.15f}{:20.15f}".format(ham_r[j,i,irpts].real, ham_r[j,i,irpts].imag)
                    else:
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
                    if self.prec:
                        line += " {:22.15e} {:22.15e}".format(ham_r[j,i,irpts].real, ham_r[j,i,irpts].imag)
                    else:
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


def main(argv=None, for_cli=False):
    if argv is None:
        argv = sys.argv[1:]
    progname = "symmwanier wannierize" if for_cli else "python wannierize.py"

    parser = argparse.ArgumentParser(
        prog=progname,
        description="Compute maximally localized Wannier functions using symmetry information."
    )

    parser.add_argument(
        "-s", "--symmetry",
        action="store_true",
        help="Generate Wannier functions in the conventional manner using prefix.isym, prefix.immn, prefix.iamn, and prefix.ieig"
    )

    parser.add_argument(
        "-S", "--site-symmetry",
        action="store_true",
        help="Generate symmetry-adapted Wannier functions using prefix.isym, prefix.immn, prefix.iamn, and prefix.ieig"
    )

    parser.add_argument(
        "-H", "--high-precision",
        action="store_true",
        help="Enable high-precision output for hr.dat and tb.dat"
    )

    parser.add_argument(
        "--log-level",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set logging level (default: INFO)"
    )

    parser.add_argument(
        "--optimize-memory",
        action="store_true",
        dest="optimize_memory_usage",
        help="Prioritize memory efficiency over speed (default: optimize for speed)"
    )

    def _str2bool(s):
        return str(s).lower() in ("1", "true", "yes", "on", "t")

    parser.add_argument(
        "--snap-kp",
        type=_str2bool, default=True, metavar="BOOL",
        help="Snap .nnkp k-points to the exact i/mp_grid rationals of the "
             "Monkhorst-Pack grid (default: True). Keep the finite-digit "
             "values with: --snap-kp false"
    )

    parser.add_argument(
        "-P", "--projectability-disentangle",
        action="store_true",
        help="Use AMN projectability to define disentanglement windows instead of energy windows"
    )

    parser.add_argument(
        "prefix",
        help="Prefix name of input/output files"
    )

    args = parser.parse_args(argv)

    wann = Wannierize(prefix=args.prefix, lsym=args.symmetry, lsite_sym=args.site_symmetry, prec=args.high_precision, optimize_memory_usage=args.optimize_memory_usage, projectability_disentangle=args.projectability_disentangle, log_level=args.log_level, snap_kp=args.snap_kp)
    wann.run()


if __name__ == '__main__':
    main()
