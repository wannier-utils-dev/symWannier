#!/usr/bin/env python
"""
Verify the k-point snapping fix on the diamond (6x6x6) example.

This uses only the committed IBZ inputs (diamond_mp6.i*), so it needs symWannier
(Python) only -- no QE / wannier90 / pw2wannier90.

A Wannier tight-binding model must reproduce the ab-initio Hamiltonian H(k)
exactly at the k-points of the Wannierization mesh. We therefore Fourier
transform the model's H(k) to real space H(R) and back, and measure how well
H(k) is recovered at the mesh points.

The 6x6x6 grid has k = i/6; 1/6 = 0.16666... is not representable in the 8
decimals that wannier90 writes into the .nnkp, so the stored value is ~3e-9
below the exact double. Summed over the real-space cell in the Fourier phases
2*pi*k*R this grows to ~1e-6. Snapping k back to the exact rational i/6
(the default, --snap-kp true) restores machine precision.

Addendum -- what the numbers really mean.  The round-trip above uses the SAME k
for H(k)->H(R) and H(R)->H(k), which slightly overstates the role of the
backward transform.  Wannier interpolation (H(R)->H(k) along a band path)
always evaluates at whatever double-precision k you ask for and never re-reads
the .nnkp, so it does NOT depend on .nnkp precision.  The error lives entirely
in the forward build of H(R): the coefficients written to hr.dat/tb.dat are
themselves wrong by ~2e-7, and interpolating that H(R) even at the EXACT mesh k
still misses the true H(k) by ~1e-6.  Lines (1)/(2) below show this directly.
"""
import os
import sys
import itertools
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from symwannier.wannierize import Wannierize

# Build the symmetry-adapted Wannier model on the 6x6x6 mesh.
w = Wannierize(prefix="diamond_mp6", lsym=True, lsite_sym=True, snap_kp=True)
w.run()

# H(k) = U(k)^dagger diag(eps(k)) U(k) -- what symWannier computes and writes as H(R).
H = np.einsum("kni,kn,knj->kij", np.conj(w.Umat), w.eig.eig, w.Umat)
k_snap = w.kpts                    # exact i/mp_grid (bit-identical to QE internal k)
k_trunc = np.round(k_snap, 8)      # the finite-digit values as stored in .nnkp

N = int(w.win.mp_grid[0])
nk = len(H)
R = np.array(list(itertools.product(range(-(N // 2), N - N // 2), repeat=3)))


def roundtrip_err(k):
    kr = np.einsum("ka,ra->kr", k, R)
    ham_r = np.einsum("kmn,kr->rmn", H, np.exp(-2j * np.pi * kr)) / nk   # H(k) -> H(R)
    ham_k = np.einsum("rmn,kr->kmn", ham_r, np.exp(2j * np.pi * kr))     # H(R) -> H(k)
    return np.max(np.abs(ham_k - H))


print()
print("Fourier round-trip error  max| H(R)->H(k) - H(k) |  at the mesh points:")
print("  --snap-kp false  (finite-digit .nnkp k) : {:.2e}".format(roundtrip_err(k_trunc)))
print("  --snap-kp true   (exact i/mp_grid k)    : {:.2e}".format(roundtrip_err(k_snap)))
print()
print("max |k_snap - k_trunc| = {:.2e}".format(np.max(np.abs(k_snap - k_trunc))))

# --- Addendum: the error is in H(R) itself; the backward transform is double
# --- precision and does not depend on the .nnkp.  Build H(R) with the two
# --- forward k, and in BOTH cases interpolate back at the EXACT mesh k.


def build_HR(k_build):
    """Forward transform H(k) -> H(R), as write_hr / write_tb do."""
    kr = np.einsum("ka,ra->kr", k_build, R)
    return np.einsum("kmn,kr->rmn", H, np.exp(-2j * np.pi * kr)) / nk


def interp_at_exact(ham_r):
    """Backward transform H(R) -> H(k) at the EXACT (double-precision) mesh k."""
    kr = np.einsum("ka,ra->kr", k_snap, R)
    return np.einsum("rmn,kr->kmn", ham_r, np.exp(2j * np.pi * kr))


HR_trunc = build_HR(k_trunc)   # --snap-kp false: H(R) built from 8-digit .nnkp k (the bug)
HR_snap = build_HR(k_snap)     # --snap-kp true : H(R) built from exact i/mp_grid k

print()
print("(1) deliverable H(R) coefficients (written to hr.dat/tb.dat), no interpolation:")
print("      max| H(R)_false - H(R)_true | = {:.2e}".format(np.max(np.abs(HR_trunc - HR_snap))))
print()
print("(2) interpolate H(R) back at the EXACT double-precision mesh k, vs true H(k):")
print("      --snap-kp false  (H(R) from 8-digit .nnkp k) : {:.2e}".format(
    np.max(np.abs(interp_at_exact(HR_trunc) - H))))
print("      --snap-kp true   (H(R) from exact i/mp_grid ) : {:.2e}".format(
    np.max(np.abs(interp_at_exact(HR_snap) - H))))
