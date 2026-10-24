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
import itertools
import textwrap
import datetime
import argparse
import os
import sys
import logging

from symwannier.win import Win
from symwannier.mmn import Mmn
from symwannier.amn import Amn
from symwannier.sym import Sym
from symwannier.eig import Eig
from symwannier.nnkp import Nnkp
from symwannier.timedata import TimeData
from symwannier.kdist import KDist, world_comm

class _TangentSymmetrizer:
    """Projection of a tangent direction onto those that keep the rotations symmetric.

    A direction d holds an anti-Hermitian num_wann x num_wann matrix per k-point,
    in the local frame of Umat(k): the rotations move to Umat(k) exp(a d(k)). The
    full-BZ states are defined as the images of the irreducible ones with the band
    index untouched, and Umat as the irreducible rotations times the projector
    rotation W(g,k) = Rmat[g] diag(exp(-2 pi i k.Rshift[g])) (Amn.symmetrize_expand).
    Keeping that relation through the step needs

        d(gk) = W(g,k)^dagger d(k) W(g,k)        (conjugated when g is antiunitary)

    across a star and the same with d(k) on both sides within the little group
    G_k. The projection takes every direction of a star back to its irreducible
    point, averages the star and then G_k, and maps the result out again. It is
    orthogonal in <A,B> = sum_k Re Tr(A_k^dagger B_k), so a projected descent
    direction stays a descent direction.
    """

    def __init__(self, sym, amn, kd):
        Rmat, Rshift, _ = amn.projection_sym_mat()

        def rot(isym, k):
            return Rmat[isym] * np.exp(-2j*np.pi * np.einsum("a,na->n", k, Rshift[isym]))[None, :]

        self.kd = kd
        self.nks = sym.nks
        ks = range(kd.k0, kd.k1)
        self.iks = sym.equiv[kd.k0:kd.k1]
        self.trev = np.array([sym.t_rev[sym.equiv_sym[k]] == 1 for k in ks], dtype=bool)
        self.W = np.array([rot(sym.equiv_sym[k], sym.irr_kpoints[sym.equiv[k]]) for k in ks])
        self.count = np.bincount(sym.equiv, minlength=sym.nks).astype(float)
        # every (irreducible k-point, little-group operation) pair, flat, since the
        # little groups differ in size and most k-points have a small one
        pairs = [ (i, h) for i in range(sym.nks) for h in sym.little_group(i) ]
        self.pair_k = np.array([ i for i, _ in pairs ], dtype=int)
        self.pair_W = np.array([ rot(h, sym.irr_kpoints[i]) for i, h in pairs ])
        self.pair_trev = np.array([ sym.t_rev[h] == 1 for _, h in pairs ], dtype=bool)
        self.order = np.bincount(self.pair_k, minlength=sym.nks).astype(float)

    @staticmethod
    def _forward(W, X, trev):
        """W^dagger X W, conjugated where trev"""
        Y = np.einsum("kml,kmn,knj->klj", np.conj(W), X, W, optimize=True)
        return np.where(trev[:, None, None], np.conj(Y), Y)

    def __call__(self, d):
        # back to the irreducible point, and the average over the star
        dc = np.where(self.trev[:, None, None], np.conj(d), d)
        back = np.einsum("kml,kln,kjn->kmj", self.W, dc, np.conj(self.W), optimize=True)
        irr = np.zeros((self.nks,) + d.shape[1:], dtype=complex)
        np.add.at(irr, self.iks, back)
        irr = self.kd.sum(irr) / self.count[:, None, None]
        # the average over the little group
        avg = np.zeros_like(irr)
        np.add.at(avg, self.pair_k, self._forward(self.pair_W, irr[self.pair_k], self.pair_trev))
        avg /= self.order[:, None, None]
        # out to the star
        return self._forward(self.W, avg[self.iks], self.trev)


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
    minimizer : {"cg", "sd", "lbfgs"}, optional
        How the spread is minimized: conjugate gradients with wannier90's parabolic
        line search, as wannier90 does (default), steepest descent with the same
        line search, or L-BFGS, which needs the fewest iterations and can reach
        lower minima than wannier90 does. See update_mmn.
    comm : mpi4py communicator, optional
        Distribute the k-points of the full mesh over the ranks of comm. Mmn and
        everything derived from it, and the rotations Umat and Umat_opt, then hold
        only this rank's block of k-points, while Amn, the eigenvalues and the
        windows are kept whole on every rank. The spreads, centers and outputs are
        the same on every rank, and only rank 0 writes files. Two things are left
        to the caller, as main() does them: silencing the other ranks, which print
        too, and calling comm.Abort() on an exception. An error found at one
        k-point is raised only on the rank that holds it, and the others would wait
        for that rank in their next exchange for as long as the job runs.
    """

    def __init__(self, prefix, lsym=False, lsite_sym=False, prec=False, optimize_memory_usage=False, projectability_disentangle=False, log_level=None, log=None, snap_kp=True, comm=None, minimizer="cg"):
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
        self.kd = KDist(self.nnkp.nk, comm)
        if not self.kd.serial:
            self.log.info("k-points distributed over %d MPI ranks, %d to %d per rank",
                          self.kd.size, min(self.kd.counts), max(self.kd.counts))
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
        mmn = Mmn(file_mmn = prefix+"." + ext + "mmn", nnkp=self.nnkp, sym=self.sym, log=self.log, kd=self.kd)
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
        self.kb2k = mmn.kb2k   # index of k+b, for every k-point of the mesh
        # U(k+b) for the local k-points, wherever k+b is held
        self._kb_plan = self.kd.plan(self.kd.local(self.kb2k))

        self.time.start_clock("amn read")
        self.amn = Amn(prefix+"." + ext + "amn", nnkp=self.nnkp, sym=self.sym, log=self.log)
        self.time.stop_clock("amn read")
        self.eig = Eig(prefix+"." + ext + "eig", sym=self.sym, log=self.log)
        self.log.debug(f"Loaded Amn: shape={self.amn.amn.shape}, Eig: shape={self.eig.eig.shape}")

        if self.sym is not None:
            # the expansion cannot reproduce a band whose degenerate multiplet is cut
            # by num_bands; warn when the outer window reaches above those bands
            _, e_max = mmn.check_bands(eig=self.eig.eig)
            if e_max is not None and self.win.dis_win_max > e_max:
                self.log.warning(
                    "  Warning: dis_win_max = %g of %s.win is above %.4f eV, so the "
                    "disentanglement may pick up bands the expansion does not reproduce",
                    self.win.dis_win_max, prefix, e_max)
        self._setup_guide_centres()

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
        # p_mk = sum_n |<psi_mk|g_n>|^2 is bounded by 1 only for orthonormal trial
        # orbitals. The projections written by hand in the win file overlap, and then
        # the thresholds dis_proj_min and dis_proj_max, which live in [0,1], mean
        # nothing: everything with any weight ends up frozen.
        self.projectability_orthonormal = np.max(self.projectability) <= 1 + 1e-6
        if not self.projectability_orthonormal:
            self.log.info(
                "projectability reaches %.3f, so the projections are not orthonormal; "
                "only the energy windows are meaningful for this calculation",
                np.max(self.projectability))
        self._validate_inputs()
        self.Umat_opt = None

        self.file_hr_dat = prefix + "_py_hr.dat"
        self.file_tb_dat = prefix + "_py_tb.dat"

        self.disentangle = (self.num_bands > self.num_wann)

        if minimizer not in ("sd", "cg", "lbfgs"):
            raise ValueError("minimizer must be 'sd', 'cg' or 'lbfgs', not {!r}".format(minimizer))
        self.minimizer = minimizer
        # the minimization stops when the gradient, root mean square over the
        # k-points of dOmega/dW(k) (A^2 per radian), falls below this
        self.grad_thr = 1e-4
        self.lbfgs_memory = 10   # steps the L-BFGS minimization remembers
        self.num_cg_steps = 5    # conjugate steps between restarts, wannier90's default
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
        self._reset_minimizer()
        self.time.start_clock("wannierization")
        for i in range(self.win.num_iter):
            self.log.debug(f"Wannierization iteration {i+1}/{self.win.num_iter}")
            omega_prev = omega
            dw = self.calc_dw()
            omega = self.update_mmn(dw, omega)
            if self._converged:
                self.log.debug("Wannierization converged")
                print("  convergence has been achieved (gradient {:.2e} A^2)".format(self.grad_norm))
                break
            if self._stalled:
                print("  no step lowers the spread any more (gradient {:.2e} A^2); stopping".format(self.grad_norm))
                break
            self.log.debug(f"Iteration {i+1}: omega_tot = {omega:.8f}, delta_omega = {abs(omega - omega_prev):.2e}, gradient = {self.grad_norm:.2e}")
            print("! {:5d}-th  omega_tot = {:15.8f}    time: {:12.6f}".format(i+1, omega, self.time.get_time()))
            self.show_spreads()
            self.calc_omega_detail(self.mmn)
        self.time.stop_clock("wannierization")
        print("! Final:    omega_tot = {:15.8f}    time: {:12.6f}".format(omega, self.time.get_time()))

    def post_process(self):
        if self.disentangle:
            self.Umat = np.einsum("kmn,knp->kmp", self.Umat_opt, self.Umat, optimize=True)

        if self.lsite_sym:
            Umat_irk = self.kd.allgather(self.Umat)[self.sym.iks2ik[:]]
            Umat_irk = self.amn.symmetrize_Gk(Umat_irk)
            self.Umat = self.kd.local(self.amn.symmetrize_expand(Umat_irk))

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

        if self.mmn0.shape[0] != self.kd.nk_loc or self.mmn0.shape[1] != self.nnkp.nb:
            raise ValueError("mmn shape {} inconsistent with nk={} or nb={}".format(self.mmn0.shape, self.kd.nk_loc, self.nnkp.nb))

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
        # U^k(dagger)[k,1,m,l] Mmn[k,b,l,p] U^(k+b)[k,b,p,n]; batched matmul uses BLAS.
        # Both are this rank's block of k-points, and U^(k+b) is fetched from
        # whichever rank holds k+b
        Umat_dag = np.conj(Umat).swapaxes(1, 2)[:, None, :, :]
        return np.matmul(np.matmul(Umat_dag, mmn), self._kb_plan.fetch(Umat))

    def init_Umat_and_Mmn(self):
        """
        initialize self.Umat (U^k) and calculate self.mmn (Mmn = U^k Mmn^(0) U^k+b)
        Mmn0 (input) -> Mmn1 (disentangle) -> Mmn
        """
        self.log.debug(f"init_Umat_and_Mmn: disentangle={self.disentangle}")
        if self.disentangle:
            self.log.info("init Umatrix and Mmn using Umatrix_opt")
            # Amn.Umat symmetrizes through the irreducible k-points, so it needs
            # Umat_opt on the whole mesh
            Umat_opt = self.kd.allgather(self.Umat_opt)
            self.Umat = self.kd.local(self.amn.Umat(index_win=self.index_win, Umat_opt=Umat_opt))
            self.mmn1 = self.update_Mmn_by_Umat(self.mmn0, self.Umat_opt)
            self.mmn = self.update_Mmn_by_Umat(self.mmn1, self.Umat)
            self.calc_omega_detail(self.mmn)

        else:  # without disentangle
            self.log.info("init Umatrix and Mmn")
            if hasattr(self, "amn"):
                self.Umat = self.kd.local(self.amn.Umat())
                self.mmn1 = self.mmn0
                self.mmn = self.update_Mmn_by_Umat(self.mmn1, self.Umat)

            else:
                self.Umat = np.zeros([self.kd.nk_loc, self.num_bands, self.num_bands], dtype=complex)
                for n in range(self.num_bands):
                    self.Umat[:,n,n] = 1
                self.mmn1 = self.mmn0
                self.mmn = self.mmn1

    def _setup_guide_centres(self):
        """Pick the centres the branch of Im ln M_nn is taken closest to.

        The Wannier center is -1/N sum_kb w_b b Im ln M_nn, and Im ln is multivalued.
        Taking the principal branch is only safe while |b.r| stays below pi, i.e. while
        the center is within about half the supercell of the origin; past that the
        branch wraps, the centers land in the wrong cell and the minimization fails.
        wannier90 avoids it with guiding_centres, which fits guide centers to the data
        every num_guide_cycles iterations. Here the projection centers are already
        known - the symmetry expansion needs them - so they serve as the guide
        directly, and the safe window moves from the origin to them.

        This is on unless the win file says guiding_centres = .false., which restores
        the principal branch. That default is deliberately not wannier90's, whose own
        default is .false.: with the guide the results of the bundled systems move by
        at most 2e-14, while without it a projection written at a periodic image far
        from the origin silently fails to localize.
        """
        self.guide_centres = None
        if getattr(self.win, "guiding_centres", None) is False:
            self.log.info("guiding_centres = .false.: taking the principal branch of Im ln M_nn")
            return
        centers = getattr(self.amn, "projection_centers", None)
        if centers is None:
            self.log.info("no projection centers available: taking the principal branch "
                          "of Im ln M_nn")
            return
        self.guide_centres = np.array(centers) @ self.nnkp.a      # crystal -> cartesian
        self.log.info("taking the branch of Im ln M_nn closest to the projection centers")

    def _imln_mnn(self, mnn):
        """Im ln M_nn, on the branch closest to the guide centres."""
        if self.guide_centres is None:
            return np.log(mnn).imag
        # the phase a Wannier function sitting on the guide center would have
        pred = - np.einsum("ba,na->bn", self.bvec, self.guide_centres, optimize=True)
        return pred[np.newaxis,:,:] + np.log(mnn * np.exp(-1j*pred)[np.newaxis,:,:]).imag

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
        self._cached_imlnmnn = self._imln_mnn(self._cached_mnn)
        self._cached_r = -1/self.nk * self.kd.sum(np.einsum("b,ba,kbn->na", self.wb, self.bvec, self._cached_imlnmnn, optimize=True))
        
        r2a = np.sum(self.wb)*self.nk - self.kd.sum(np.einsum("b,kbn,kbn->n", self.wb, self._cached_mnn, np.conj(self._cached_mnn), optimize=True))
        r2b = self.kd.sum(np.einsum("b,kbn->n", self.wb, self._cached_imlnmnn**2))
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
        OmegaI = self.kd.sum(np.einsum("b, kb->", self.wb, self.num_wann - mmn2, optimize=True))/self.nk

        mnn2 = np.einsum("kbnn,kbnn->kb", mmn, np.conj(mmn), optimize=True).real
        OmegaOD = self.kd.sum(np.einsum("b, kb->", self.wb, mmn2 - mnn2, optimize=True))/self.nk

        # Reuse cached values if available (set by calc_omega)
        if hasattr(self, '_cached_mnn') and hasattr(self, '_cached_imlnmnn'):
            mnn = self._cached_mnn
            imlnmnn = self._cached_imlnmnn
        else:
            mnn = np.einsum("kbnn->kbn", mmn, optimize=True)
            self._assert_nonzero_mnn(mnn, "calc_omega_detail")
            imlnmnn = self._imln_mnn(mnn)
        
        r = -1/self.nk * self.kd.sum(np.einsum("b,ba,kbn->na", self.wb, self.bvec, imlnmnn, optimize=True))
        qn = imlnmnn + np.einsum("ba, na->bn", self.bvec, r, optimize=True)[np.newaxis,:,:]
        OmegaD = self.kd.sum(np.einsum("b, kbn->", self.wb, qn**2, optimize=True))/self.nk

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
            imlnmnn = self._imln_mnn(mnn)
            r = -1/self.nk * self.kd.sum(np.einsum("b,ba,kbn->na", self.wb, self.bvec, imlnmnn, optimize=True))
        
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
        # Batched eigendecomposition for all k-points
        e, v = np.linalg.eigh(1j*dw)
        # Vectorized matrix multiplication: V @ diag(exp(-i*lambda)) @ V^H
        exp_e = np.exp(-1j * e)  # (nk, num_wann)
        expdw = np.einsum("kab,kb,kcb->kac", v, exp_e, np.conj(v), optimize=True)
        return expdw

    def _ip(self, a, b):
        """<A,B> = sum_k Re Tr(A_k^dagger B_k), over the whole mesh."""
        return self.kd.sum(np.sum(np.real(np.conj(a) * b)))

    def _reset_minimizer(self):
        """Forget the L-BFGS history, at the start of the minimization."""
        self._lbfgs_s, self._lbfgs_y = [], []
        self._lbfgs_last = None          # (step, gradient) of the previous iteration
        self._trial_step = 1.0          # trial step of the parabolic line search
        self._converged = False          # the gradient fell below grad_thr
        self._stalled = False            # no step along the direction lowered Omega
        self.grad_norm = None            # rms over the k-points of dOmega/dW(k)
        self._cg_dir = None              # conjugate gradients: previous direction,
        self._cg_gg = None               # |gradient|^2 there,
        self._cg_count = 0               # and conjugate steps since the last restart
        if self.lsite_sym and not hasattr(self, "_tangent_symmetrizer"):
            self._tangent_symmetrizer = _TangentSymmetrizer(self.sym, self.amn, self.kd)

    def _symmetrize_tangent(self, d):
        """The direction kept to symmetric rotations for symmetry-adapted Wannier
        functions (-S); unchanged otherwise, -s included, where symmetry serves
        only to build the inputs."""
        if not self.lsite_sym:
            return d
        return self._tangent_symmetrizer(d)

    def _omega_along(self, d, a):
        """Omega at Umat exp(a d), with the rotations and overlaps it is evaluated at."""
        umat = np.einsum("kmn, knl->kml", self.Umat, self.calc_expdw(a * d), optimize=True)
        mmn = self.update_Mmn_by_Umat(self.mmn1, umat)
        return self.calc_omega(mmn), umat, mmn

    def _parabolic_step(self, d, omega, slope):
        """wannier90's line search: Omega at a trial step, then the minimum of the
        parabola through Omega(0), its slope and that value. Returns
        (omega, umat, mmn, step) of the best point found, or None."""
        a_t = self._trial_step
        om_t, u_t, m_t = self._omega_along(d, a_t)
        curv = 2 * (om_t - omega - slope * a_t) / a_t**2
        if curv > 0:
            a = -slope / curv
            om, u, m = self._omega_along(d, a)
            if om < min(omega, om_t):
                self._trial_step = a
                return om, u, m, a
        if om_t < omega:
            return om_t, u_t, m_t, a_t
        a = a_t
        for _ in range(30):
            a *= 0.5
            om, u, m = self._omega_along(d, a)
            if om < omega:
                self._trial_step = a
                return om, u, m, a
        return None

    def _armijo_step(self, d, omega, slope):
        """Quasi-Newton step 1, backtracking along a parabola until Omega falls
        enough. Returns (omega, umat, mmn, step), or None."""
        a = 1.0
        while a > 1e-6:
            om, u, m = self._omega_along(d, a)
            if om <= omega + 1e-4 * a * slope:
                return om, u, m, a
            a = max(0.1 * a, min(0.5 * a, -slope * a * a / (2 * (om - omega - slope * a))))
        return None

    def _lbfgs_step(self, g, omega):
        """L-BFGS: the direction from the last lbfgs_memory steps and gradients by
        the two-loop recursion, the quasi-Newton step 1 with Armijo backtracking,
        and wannier90's parabola on the first step and after a restart.
        Returns (direction, step found or None)."""
        if self._lbfgs_last is not None:
            s_prev, g_prev = self._lbfgs_last
            y = g - g_prev
            sy = self._ip(s_prev, y)
            if sy > 1e-12 * np.sqrt(self._ip(s_prev, s_prev) * self._ip(y, y)):
                self._lbfgs_s.append(s_prev)
                self._lbfgs_y.append(y)
                if len(self._lbfgs_s) > self.lbfgs_memory:
                    self._lbfgs_s.pop(0)
                    self._lbfgs_y.pop(0)

        # two-loop recursion
        q = g.copy()
        history = []
        for s_i, y_i in zip(reversed(self._lbfgs_s), reversed(self._lbfgs_y)):
            rho = 1 / self._ip(y_i, s_i)
            a_i = rho * self._ip(s_i, q)
            q -= a_i * y_i
            history.append((rho, a_i, s_i, y_i))
        if self._lbfgs_s:
            q *= self._ip(self._lbfgs_s[-1], self._lbfgs_y[-1]) / self._ip(self._lbfgs_y[-1], self._lbfgs_y[-1])
        for rho, a_i, s_i, y_i in reversed(history):
            q += (a_i - rho * self._ip(y_i, q)) * s_i
        d = self._symmetrize_tangent(-q)
        slope = self._ip(g, d)
        if slope >= 0:
            # the history no longer gives a descent direction
            self._lbfgs_s, self._lbfgs_y = [], []
            d = -g
            slope = self._ip(g, d)
        if slope == 0:
            return d, None

        step = self._armijo_step(d, omega, slope) if self._lbfgs_s else None
        if step is None:
            if self._lbfgs_s:
                self.log.debug("update_mmn: quasi-Newton step failed, dropping the history")
                self._lbfgs_s, self._lbfgs_y = [], []
                d = -g
                slope = self._ip(g, d)
            step = self._parabolic_step(d, omega, slope)
        return d, step

    def _cg_step(self, g, omega):
        """Conjugate gradients as wannier90 does them: Fletcher-Reeves, restarted
        to steepest descent after num_cg_steps conjugate steps, when the
        coefficient exceeds 3 and when the direction stops descending; the parabolic
        line search. Returns (direction, step found or None)."""
        gg = self._ip(g, g)
        beta = 0.0
        if self._cg_dir is not None and self._cg_count < self.num_cg_steps and self._cg_gg > 0:
            beta = gg / self._cg_gg
            if beta > 3:
                beta = 0.0
        d = -g if beta == 0 else self._symmetrize_tangent(-g + beta * self._cg_dir)
        slope = self._ip(g, d)
        if slope >= 0:
            beta, d = 0.0, -g
            slope = self._ip(g, d)
        self._cg_count = 0 if beta == 0 else self._cg_count + 1
        self._cg_dir, self._cg_gg = d, gg
        if slope == 0:
            return d, None
        return d, self._parabolic_step(d, omega, slope)

    def update_mmn(self, dw, omega):
        """One step of the spread minimization.

        dw is the descent direction at the current point, minus the gradient of
        Omega in <A,B> = sum_k Re Tr(A_k^dagger B_k): along Umat exp(a d) the
        slope at a = 0 is -<dw, d>. Two minimizers, chosen by self.minimizer:

        "cg"     conjugate gradients, as wannier90 does them (Fletcher-Reeves,
                 restarted every num_cg_steps), with wannier90's line search: Omega
                 at a trial step, then the minimum of the parabola through Omega(0),
                 its slope and that value; the trial step is the last step taken.
                 The default.
        "sd"     steepest descent along dw with the same line search; slow where
                 the spread is badly conditioned.
        "lbfgs"  L-BFGS, with the tangent vectors kept in the local frame of each
                 Umat(k) as wannier90 keeps its conjugate-gradient directions. It
                 needs far fewer iterations, and goes on to lower minima.

        For symmetry-adapted Wannier functions (-S), the gradient and the
        direction are projected onto the directions that keep the rotations
        symmetric. The spread can have a symmetric stationary point that is not a
        minimum - the atom-centred s, p, d orbitals of Fe are one - and a
        minimizer that follows round-off out of it ends in symmetry-broken
        orbitals, differently for every summation order. The projection keeps the
        minimization on the symmetric solutions. With -s the minimization is
        unconstrained, as in wannier90, and may leave such a point.

        Convergence is decided here, before the step: when the gradient, root mean
        square over the k-points of dOmega/dW(k), is below grad_thr, nothing moves
        and self._converged is set; when no step along the direction lowers Omega,
        self._stalled is.

        Returns the new Omega; Umat, mmn and the centres and spreads are updated.
        """
        g = -self._symmetrize_tangent(dw)
        # converged before stepping: the rms gradient per k-point, which unlike the
        # norm of the whole array does not grow with the mesh (dw carries 1/N)
        self.grad_norm = np.sqrt(self.nk * self._ip(g, g))
        if self.grad_norm < self.grad_thr:
            self._converged = True
            return omega
        if self.minimizer == "lbfgs":
            d, step = self._lbfgs_step(g, omega)
        elif self.minimizer == "cg":
            d, step = self._cg_step(g, omega)
        else:
            d = -g
            slope = self._ip(g, d)
            step = self._parabolic_step(d, omega, slope) if slope < 0 else None
        if step is None:
            self.log.debug("update_mmn: no step lowers Omega")
            self.calc_omega(self.mmn)
            self._stalled = True
            return omega
        omega_new, self.Umat, self.mmn, a = step
        self.log.debug(f"update_mmn ({self.minimizer}): step {a:.4g}, omega={omega_new:.8f}")
        # the caches calc_dw and calc_omega_detail read must belong to the accepted point
        omega_new = self.calc_omega(self.mmn)
        self._lbfgs_last = (a * d, g)
        return omega_new

    def dis_window(self):
        """Define energy windows for disentanglement."""
        # array with self.nband size
        self.index_froz = (self.eig.eig > self.win.dis_froz_min) & (self.eig.eig < self.win.dis_froz_max)
        self.index_win = (self.eig.eig > self.win.dis_win_min) & (self.eig.eig < self.win.dis_win_max)
        self.index_nfroz = self.index_win & (~self.index_froz)
        self.len_nfroz = np.sum(self.index_nfroz, axis=1)
        self.ndimwin = np.sum(self.index_win, axis=1)
        self.ndimfroz = np.sum(self.index_froz, axis=1)
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
        if not self.projectability_orthonormal:
            raise ValueError(
                "projectability disentanglement needs orthonormal projections, and "
                "p_mk = sum_n |<psi_mk|g_n>|^2 reaches {:.3f} here. Use the atom_proj "
                "interface of pw2wannier90, which orthonormalizes the projectors, or "
                "the energy windows instead.".format(np.max(self.projectability)))

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

        if not 0 <= dis_proj_min <= dis_proj_max <= 1:
            raise ValueError("projectability thresholds must satisfy 0 <= dis_proj_min <= dis_proj_max <= 1")

        energy_outer = ((self.eig.eig >= self.win.dis_win_min) &
                        (self.eig.eig <= self.win.dis_win_max))
        energy_frozen = np.zeros_like(energy_outer)
        if self.win.has_dis_froz_window:
            energy_frozen = ((self.eig.eig >= self.win.dis_froz_min) &
                             (self.eig.eig <= self.win.dis_froz_max))

        # Match Wannier90 dis_windows_proj: the outer energy window is always
        # applied, while the frozen subspace is the union of the energy-frozen
        # states and states above dis_proj_max. Energy-frozen states are kept
        # even if their projectability is below dis_proj_min.
        self.index_froz = energy_outer & (energy_frozen | (self.projectability >= dis_proj_max))
        self.index_win = energy_outer & (energy_frozen | (self.projectability >= dis_proj_min))
        self.index_nfroz = self.index_win & (~self.index_froz)
        self.len_nfroz = np.sum(self.index_nfroz, axis=1)
        self.ndimwin = np.sum(self.index_win, axis=1)
        self.ndimfroz = np.sum(self.index_froz, axis=1)
        self.log.info(f"num_inner (min, max, mean): {np.min(self.ndimfroz)}, {np.max(self.ndimfroz)}, {np.mean(self.ndimfroz)}")
        self.log.info(f"num_outer (min, max, mean): {np.min(self.ndimwin)}, {np.max(self.ndimwin)}, {np.mean(self.ndimwin)}")

    def dis_project(self):
        """Calculate self.Umat_opt and update mmn0 from amn file."""
        self.log.debug("dis_project: starting projection")
        assert hasattr(self, "amn"), "file_amn is not specified"
        # Umat_opt[:nk_loc, :num_bands, :num_wann], this rank's block of k-points
        self.Umat_opt = self.kd.local(self.amn.Umat(index_win = self.index_win)).copy()

        # dis_proj_froz in wannier90
        # wannier90 uses ortho-fix which is not implemented here
        for kl in range(self.kd.nk_loc):
            k = self.kd.k0 + kl   # index_win, index_froz and ndimfroz cover the whole mesh
            # cqpq = Q_inner P_G Q_inner in R2_Eq.(27)
            Umatk = self.Umat_opt[kl,self.index_win[k,:],:]
            lnfroz = (~self.index_froz[ k, self.index_win[k,:] ]).astype(int)
            cqpq = np.einsum("m,ml,nl,n->mn", lnfroz, Umatk, np.conj(Umatk), lnfroz, optimize=True)
            # check hermiticity
            #assert np.sum(np.abs(cqpq - np.transpose(np.conj(cqpq)))) < 1e-5  # check hermiticity

            if cqpq.shape[0] == 0:
                # the outer window is empty, which the eigensolver reports as an
                # opaque LAPACK error about the index range
                raise ValueError(
                    "the disentanglement window at k-point {} contains no state; "
                    "widen dis_win_min / dis_win_max".format(k+1))
            e, v = scipy.linalg.eigh(-cqpq)  # "-" sign for descending order
            self.Umat_opt[kl,:,:] = 0
            self.Umat_opt[kl, self.index_froz[k,:], :self.ndimfroz[k]] = np.identity(self.ndimfroz[k])
            self.Umat_opt[kl, self.index_win[k,:], self.ndimfroz[k]:self.num_wann] = v[:, :self.num_wann - self.ndimfroz[k]]

            # check unitary
            #assert np.allclose( np.matmul(np.transpose(np.conj(self.Umat_opt[:,:,k])), self.Umat_opt[:,:,k]), np.eye(self.num_wann) )

        self.update_mmn0U()

    def update_mmn0U(self):
        """
        update mmn0U = mmn0[num_bands, num_win] x Umat_opt[num_win, num_wann]
        """
        self.mmn0U = np.matmul(self.mmn0, self._kb_plan.fetch(self.Umat_opt))  # Eq. (61)

    def init_zmatrix(self):
        """
        calculate Zmatrix (zmat)
        Zmatrix = sum_b wb <u^(0)_mk | P^(i)_k+b | u^0_nk>
                = sum_{bl} wb <u^(0)_mk | u^(i)_lk+b > < u^(i)_lk+b | u^0_nk >
        mu : Mmn_k U^(i)_k+b = <u^(0)_mk | u^(i)_lk+b >
        """
        # mu = mmn0[num_bands, num_win] x Umat_opt[num_win, num_wann]
        # zmat = mu x mu^\dagger
        nk_loc = self.kd.nk_loc
        zmat_froz = np.zeros([nk_loc, self.num_bands, self.num_bands], dtype=complex)
        zmat_nfroz = np.zeros([nk_loc, self.num_bands, self.num_bands], dtype=complex)
        # zmat[k,m,n] = sum_{b,l} wb[b] mu[k,b,m,l] conj(mu[k,b,n,l]); contract over (b,l) with batched matmul
        mu = self.mmn0U.transpose(0, 2, 1, 3).reshape(nk_loc, self.num_bands, -1)
        mu_wb = (self.mmn0U * self.wb[None, :, None, None]).transpose(0, 2, 1, 3).reshape(nk_loc, self.num_bands, -1)
        zmat = np.matmul(mu_wb, np.conj(mu).swapaxes(1, 2))
        for kl in range(nk_loc):
            k = self.kd.k0 + kl
            zmat_froz[kl,:self.ndimfroz[k],:self.ndimfroz[k]] = zmat[kl, self.index_froz[k,:], :][:, self.index_froz[k,:]]
            zmat_nfroz[kl,:self.len_nfroz[k],:self.len_nfroz[k]] = zmat[kl, self.index_nfroz[k,:], :][:, self.index_nfroz[k,:]]

        return zmat_froz, zmat_nfroz

    def calc_womegaI(self):
        mmn = np.matmul(np.conj(self.Umat_opt).swapaxes(1, 2)[:, None, :, :], self.mmn0U)
        return self.num_wann * np.sum(self.wb) - self.kd.sum(np.einsum("b,kbmn,kbmn->", self.wb, mmn, np.conj(mmn), optimize=True)).real/self.nk

    def dis_extract(self):
        for i in range(self.win.dis_num_iter):
            zmat_froz, zmat_nfroz = self.init_zmatrix()
            if not i == 0:
                # update zmatrix
                zmat_nfroz = self.win.dis_mix_ratio * zmat_nfroz + (1-self.win.dis_mix_ratio) * zmat_nfroz_old

            # zmat_froz is zero outside the frozen block, so the full trace equals the block trace
            wkomegaI1 = self.num_wann * np.sum(self.wb) - np.trace(zmat_froz, axis1=1, axis2=2).real

            for kl in range(self.kd.nk_loc):
                k = self.kd.k0 + kl
                # Here, we consider -Z because we need eigvals of Z in descending order
                e, v = scipy.linalg.eigh(-zmat_nfroz[kl, :self.len_nfroz[k], :self.len_nfroz[k]])
                e = -e
                wkomegaI1[kl] -= np.sum( e[:(self.num_wann-self.ndimfroz[k])] )

                self.Umat_opt[kl, self.index_win[k,:], self.ndimfroz[k]:self.num_wann] = 0
                self.Umat_opt[kl, self.index_nfroz[k,:], self.ndimfroz[k]:self.num_wann] = v[:,:self.num_wann-self.ndimfroz[k]]

                # check unitary
                #assert np.allclose( np.matmul(np.transpose(np.conj(self.Umat_opt[:,:,k])), self.Umat_opt[:,:,k]), np.eye(self.num_wann) )
            self.update_mmn0U()
            
            womegaI1 = self.kd.sum(np.sum(wkomegaI1))/self.nk

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

    def calc_ham_r(self):
        """Return H(R) = 1/Nk sum_k e^{-2 pi i k.R} U^k+ diag(e_k) U^k and the phase factors e^{-2 pi i k.R}.

        Requires ``self.irvec`` (see ``calc_irvec``).
        """
        ham_k = np.einsum("kni,kn,knj->kij", np.conj(self.Umat), self.kd.local(self.eig.eig), self.Umat, optimize=True)
        kr = np.einsum("ka,ra->kr", self.kd.local(self.kpts), self.irvec, optimize=True)
        fac = np.exp(-2 * np.pi * 1j * kr)   # for this rank's k-points
        ham_r = self.kd.sum(np.einsum("kij,kr->ijr", ham_k, fac, optimize=True))/self.nk
        return ham_r, fac

    def write_hr(self):
        self.calc_irvec()
        ham_r, _ = self.calc_ham_r()
        if not self.kd.is_root:
            return

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
        ham_r, fac = self.calc_ham_r()

        mnn = np.einsum("kbnn->kbn", self.mmn, optimize=True)
        imlnmnn = self._imln_mnn(mnn)
        rk = - np.einsum("b,ba,kbn->ank", self.wb, self.bvec, imlnmnn, optimize=True)
        rr = self.kd.sum(np.einsum("ank,kr->anr", rk, fac, optimize=True))/self.nk
        rr2 = 1j * self.kd.sum(np.einsum("b,ba,kbmn,kr->amnr", self.wb, self.bvec, self.mmn, fac, optimize=True))/self.nk
        if not self.kd.is_root:
            return

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
        "--minimizer",
        choices=["cg", "sd", "lbfgs"], default="cg",
        help="Spread minimization: conjugate gradients with a parabolic line "
             "search, as wannier90 (default), steepest descent with the same line "
             "search, or L-BFGS, which needs the fewest iterations and can reach "
             "lower minima than wannier90 does"
    )

    parser.add_argument(
        "prefix",
        help="Prefix name of input/output files"
    )

    args = parser.parse_args(argv)

    # under mpirun / srun the k-points are distributed over the ranks; only rank 0
    # reports, and the others keep their errors
    comm = world_comm()
    log_level = args.log_level
    if comm is not None and comm.Get_rank() != 0:
        sys.stdout = open(os.devnull, "w")
        log_level = "ERROR"

    try:
        wann = Wannierize(prefix=args.prefix, lsym=args.symmetry, lsite_sym=args.site_symmetry, prec=args.high_precision, optimize_memory_usage=args.optimize_memory_usage, projectability_disentangle=args.projectability_disentangle, log_level=log_level, snap_kp=args.snap_kp, comm=comm, minimizer=args.minimizer)
        wann.run()
    except Exception:
        if comm is None:
            raise
        # an error found at one k-point is raised only on the rank that holds it, and
        # the others would wait for it in their next sum for as long as the job runs
        import traceback
        traceback.print_exc()
        sys.stderr.flush()
        comm.Abort(1)


if __name__ == '__main__':
    main()
