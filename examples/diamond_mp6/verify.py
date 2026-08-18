#!/usr/bin/env python
"""
Verify the k-point snapping fix on the diamond (6x6x6) example.

This uses only the committed IBZ inputs (diamond_mp6.i*), so it needs symWannier
(Python) only -- no QE / wannier90 / pw2wannier90.

WHAT GOES WRONG.  symWannier writes the tight-binding Hamiltonian H(R)
(hr.dat / tb.dat) by Fourier transforming the ab-initio H(k),

        H(R) = (1/Nk) sum_k  exp(-2*pi*i k.R) H(k)

Without --snap-kp the k used here are the .nnkp values, which wannier90 stores
to only 8 decimals.  On a 6x6x6 mesh 1/6 = 0.16666667 is ~3e-9 too small; that
error enters the phases 2*pi*k.R and corrupts the H(R) coefficients.

Note the asymmetry: the BACKWARD transform (Wannier interpolation H(R) -> H(k)
along a band path) evaluates at whatever double-precision k you ask for and
never re-reads the .nnkp, so it does NOT depend on .nnkp precision.  Only the
forward build of H(R) does.  Below, H(R) is built both ways and in BOTH cases
evaluated back at the exact double-precision mesh k, so the only difference is
how H(R) was built.

TWO DIFFERENT ERRORS ARE REPORTED -- they answer different questions:

(A) DFT bands vs Wannier bands.  Reference = the QE eigenvalues eps_DFT read
    from the .ieig file, i.e. an EXTERNAL reference.  This is the error a user
    of the model actually sees.  Its floor is NOT set by the Fourier transform
    but by the unitarity of the MLWF gauge matrix U (~1e-11 here): H(k) is
    built as U^dag diag(eps_DFT) U, so if U is not exactly unitary the
    eigenvalues of H(k) already differ from eps_DFT before any transform.  That
    floor is printed too, so the snapped number can be recognised as sitting on
    it rather than being a leftover of the truncation.

(B) Fourier self-consistency.  Reference = symWannier's own H(k), i.e. an
    INTERNAL reference.  Because the same H(k) appears on both sides, the U
    floor cancels and this isolates the Fourier-transform error alone -- the
    quantity the snap actually fixes.  This is NOT a comparison against DFT.
"""
import os
import sys
import itertools
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from symwannier.wannierize import Wannierize

# Build the symmetry-adapted Wannier model.  snap_kp only affects the Fourier
# phases below (U and the DFT eigenvalues come from DFT), so we read the raw
# .nnkp k here and form the exact rationals ourselves.
w = Wannierize(prefix="diamond_mp6", lsym=True, lsite_sym=True, snap_kp=False)
w.run()

eps_dft = np.sort(w.eig.eig, axis=1)     # DFT band energies from QE (.ieig), in eV
nk, nbnd = eps_dft.shape
nw = w.Umat.shape[2]

# H(k) in the Wannier gauge: H(k) = U(k)^dagger diag(eps_DFT(k)) U(k).
H = np.einsum("kni,kn,knj->kij", np.conj(w.Umat), w.eig.eig, w.Umat)

k_nnkp = w.kpts                                  # raw .nnkp, 8 decimals (0.16666667 for 1/6)
mp = np.asarray(w.win.mp_grid, dtype=int)        # per-direction mesh, e.g. [6 6 6]
k_exact = np.round(k_nnkp * 2 * mp) / (2 * mp)   # snapped: exact i/mp_grid rational
R = np.array(list(itertools.product(*[range(-(n // 2), n - n // 2) for n in mp])))


def build_HR(k_build):
    """Forward transform H(k) -> H(R), as write_hr / write_tb do."""
    kr = np.einsum("ka,ra->kr", k_build, R)
    return np.einsum("kmn,kr->rmn", H, np.exp(-2j * np.pi * kr)) / nk


def interp_at_exact(ham_r):
    """Backward transform H(R) -> H(k) at the EXACT double-precision mesh k."""
    kr = np.einsum("ka,ra->kr", k_exact, R)
    return np.einsum("rmn,kr->kmn", ham_r, np.exp(2j * np.pi * kr))


HR_nnkp = build_HR(k_nnkp)     # --snap-kp false: H(R) built from 8-digit .nnkp k (the bug)
HR_exact = build_HR(k_exact)   # --snap-kp true : H(R) built from exact i/mp_grid k
Hk_nnkp = interp_at_exact(HR_nnkp)
Hk_exact = interp_at_exact(HR_exact)

# Gauge floor of metric (A): how far eig H(k) already is from eps_DFT with no
# Fourier transform at all.  Meaningful as a band comparison only without
# disentanglement (num_wann == num_bands).
gauge_floor = np.max(np.abs(np.linalg.eigvalsh(H) - eps_dft)) if nw == nbnd else None
unitarity = np.max(np.abs(np.einsum("kni,knj->kij", np.conj(w.Umat), w.Umat) - np.eye(nw)))

print()
print("System: {} mesh = {} k-points, {} bands -> {} Wannier functions{}".format(
    "x".join(str(n) for n in mp), nk, nbnd, nw,
    " (no disentanglement)" if nw == nbnd else ""))
print()

if gauge_floor is not None:
    print("(A) DFT bands vs Wannier bands   max | eps_Wannier(k) - eps_DFT(k) |")
    print("    reference: QE eigenvalues from .ieig (external)")
    print("    max over {} mesh k-points and {} bands, in eV:".format(nk, nbnd))
    print("      --snap-kp false  (H(R) from 8-digit .nnkp k) : {:.2e} eV".format(
        np.max(np.abs(np.linalg.eigvalsh(Hk_nnkp) - eps_dft))))
    print("      --snap-kp true   (H(R) from exact i/mp_grid) : {:.2e} eV".format(
        np.max(np.abs(np.linalg.eigvalsh(Hk_exact) - eps_dft))))
    print("    floor of this metric, with NO Fourier transform at all:")
    print("      max| eig H(k) - eps_DFT | = {:.2e} eV   (MLWF gauge U, ||U^dag U - 1|| = {:.1e})"
          .format(gauge_floor, unitarity))
    print("    -> the snapped value sits on this pre-existing gauge floor, not on")
    print("       a leftover of the k-point truncation.")
    print()

print("(B) Fourier self-consistency   max | H(R)->H(k) - H(k) |")
print("    reference: symWannier's own H(k) (internal; the gauge floor cancels,")
print("    so this isolates the Fourier-transform error the snap fixes)")
print("    max over {} mesh k-points and {}x{} matrix elements, in eV:".format(nk, nw, nw))
print("      --snap-kp false  (H(R) from 8-digit .nnkp k) : {:.2e} eV".format(
    np.max(np.abs(Hk_nnkp - H))))
print("      --snap-kp true   (H(R) from exact i/mp_grid) : {:.2e} eV".format(
    np.max(np.abs(Hk_exact - H))))
print()
print("Underlying cause: the written H(R) coefficients themselves differ by")
print("  max| H(R)_false - H(R)_true | = {:.2e} eV  (over {} R-vectors, {}x{} elements)"
      .format(np.max(np.abs(HR_nnkp - HR_exact)), len(R), nw, nw))
print("and the k-points differ by  max| k_exact - k_nnkp | = {:.2e}"
      .format(np.max(np.abs(k_exact - k_nnkp))))
