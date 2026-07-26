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
