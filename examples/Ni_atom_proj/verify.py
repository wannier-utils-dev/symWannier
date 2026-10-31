#!/usr/bin/env python
"""
Verify the unitarity fix of -S on the bundled Ni_atom_proj data.

This uses only the IBZ inputs committed under tests/inputs (Ni_atom_proj.i*), so
it needs symWannier (Python) only -- no QE / wannier90 / pw2wannier90.

WHAT GOES WRONG.  With -S, post_process averages the final U(k) over the little
group G_k of each irreducible k and expands it to the full mesh,

        U(k) <- 1/|G_k| sum_h  d(h) U(k) D(h)          (d: bands, D: projections)

An average of unitary matrices is not unitary, so U^dag U != 1 afterwards, and
H(k) = U^dag diag(eps) U -- what hr.dat / tb.dat are built from -- inherits it.
The fix projects each U(k) back onto the nearest unitary matrix (U = W V^dag
from the SVD).  Below, "without projection" switches that step off in the same
checkout, so both sides run the same code otherwise.

WHY IT IS LARGE HERE.  Ni_atom_proj has num_bands = 25, and at k = (0,0,1/4) the
degenerate multiplet of band 25 is cut by num_bands, so d(h) is not unitary
there (Sym.check_repmat warns at load time).  The average is then not a
projector, and U shrinks by ~1e-5 instead of only being symmetrized.  Keeping
the outer window below the ceiling Mmn.check_bands reports (86.12 eV here)
removes band 25 from U and with it most of the effect; both windows are run.

WHAT IS REPORTED -- each part answers a different question:

(A) Unitarity of the final U   max_k | U^dag U - 1 |.  This is what the fix
    restores.  No reference other than the identity.

(B) Point-group residual of H(k) = U^dag diag(eps) U on the full 4x4x4 mesh,
        max | D^dag H(k) D - H(gk) |         (conj for time reversal)
    over every k and every operation g mapping k to a mesh point.  Reference:
    symWannier's own projection rotation matrices D (prefix.isym), which agree
    with an independent construction (wannier_utils' analytic D(g)) to 1e-15.
    This is the symmetry a user of hr.dat / tb.dat gets.  It is built from the
    final U directly; (C) shows that writing hr.dat adds nothing to it.

(D) shows two more consequences of U^dag U != 1.  U^dag U is the overlap of the
    Wannier functions at k, and H(k) = U^dag diag(eps) U takes it to be 1, so
    the eigenvalues of H(k) change when eps is measured from another energy
    zero, and they differ from those of the generalized problem H c = e S c.
"""
import gzip
import logging
import os
import shutil
import sys
import tempfile

import numpy as np
import scipy.linalg

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..", "src"))
from symwannier.amn import Amn
from symwannier.wannierize import Wannierize

inputs = os.path.join(here, "..", "..", "tests", "inputs")
prefix = "Ni_atom_proj"

nearest_unitary = Amn.nearest_unitary
Umat = Amn.Umat


def run(project, dis_win_max=None):
    """Run -S in a scratch directory; project=False switches the fix off."""
    if project:
        Amn.nearest_unitary = nearest_unitary
        Amn.Umat = Umat
    else:
        Amn.nearest_unitary = lambda self, U: U
        # the unitarity check in Amn.Umat would stop the unprojected run
        Amn.Umat = lambda self, index_win=None, check=True, Umat_opt=None: \
            Umat(self, index_win=index_win, check=False, Umat_opt=Umat_opt)
    work = tempfile.mkdtemp()
    for ext in ["nnkp", "isym", "ieig", "iamn", "immn"]:
        src = os.path.join(inputs, prefix + "." + ext)
        if os.path.exists(src):
            shutil.copy(src, work)
        else:
            with gzip.open(src + ".gz") as fi, open(os.path.join(work, prefix + "." + ext), "wb") as fo:
                shutil.copyfileobj(fi, fo)
    with open(os.path.join(inputs, prefix + ".win")) as fi, open(os.path.join(work, prefix + ".win"), "w") as fo:
        for line in fi:
            fo.write(line)
            if dis_win_max is not None and line.strip().startswith("dis_win_min"):
                fo.write("dis_win_max = {}\n".format(dis_win_max))
    cwd = os.getcwd()
    os.chdir(work)
    try:
        w = Wannierize(prefix=prefix, lsite_sym=True, prec=True, log_level=logging.ERROR)
        sys.stdout = open(os.devnull, "w")
        try:
            w.run()
        finally:
            sys.stdout.close()
            sys.stdout = sys.__stdout__
        hr = read_hr(w)
    finally:
        os.chdir(cwd)
        shutil.rmtree(work)
    Amn.nearest_unitary, Amn.Umat = nearest_unitary, Umat
    return w, hr


def read_hr(w):
    """H(k) at the exact mesh k rebuilt from the hr.dat just written (-H precision)."""
    nw, nr = w.num_wann, len(w.ndegen)
    data = np.loadtxt(w.file_hr_dat, skiprows=3 + int(np.ceil(nr / 15)))
    irvec = data[::nw * nw, :3]
    ham_r = (data[:, 5] + 1j * data[:, 6]).reshape(nr, nw, nw).transpose(0, 2, 1)
    k = np.round(w.sym.full_kpoints * w.win.mp_grid) / w.win.mp_grid
    fac = np.exp(2j * np.pi * k @ irvec.T) / w.ndegen[None, :]
    return np.einsum("kr,rij->kij", fac, ham_r)


def unitarity(w):
    U = w.Umat
    return np.max(np.abs(np.matmul(np.conj(U).swapaxes(1, 2), U) - np.eye(U.shape[2])))


def ham(w):
    return np.einsum("kni,kn,knj->kij", np.conj(w.Umat), w.eig.eig, w.Umat)


def residual(w, H):
    """max |D^dag H(k) D - H(gk)| over the full mesh; returns (max, worst irreducible k)."""
    sym = w.sym
    Rmat, Rshift, _ = w.amn.projection_sym_mat()

    def transform(Hk, isym, k):
        D = Rmat[isym] * np.exp(-2j * np.pi * Rshift[isym] @ k)[None, :]
        Ht = np.conj(D).T @ Hk @ D
        return np.conj(Ht) if sym.t_rev[isym] == 1 else Ht

    worst, kworst = 0.0, None
    for iks, k in enumerate(sym.irr_kpoints):      # g in G_k: k -> k
        Hk = H[sym.iks2ik[iks]]
        for h in sym.little_group(iks):
            d = np.max(np.abs(transform(Hk, h, k) - Hk))
            if d > worst:
                worst, kworst = d, k
    for ik in range(sym.nkf):                      # g taking k_irr to the rest of the star
        iks = sym.equiv[ik]
        d = np.max(np.abs(H[ik] - transform(H[sym.iks2ik[iks]], sym.equiv_sym[ik], sym.irr_kpoints[iks])))
        if d > worst:
            worst, kworst = d, sym.irr_kpoints[iks]
    return worst, kworst


win_cut = 86.0
runs = {}
for window in ["default", "cut"]:
    for project in [False, True]:
        runs[window, project] = run(project, dis_win_max=win_cut if window == "cut" else None)

w = runs["default", True][0]
print()
print("System: Ni_atom_proj, {} mesh = {} k-points, {} bands -> {} Wannier functions, -S".format(
    "x".join(str(n) for n in w.win.mp_grid), w.nk, w.num_bands, w.num_wann))
print("        frozen bands per k: {} to {}; num_iter = {}".format(
    np.min(w.ndimfroz), np.max(w.ndimfroz), w.win.num_iter))
print()
label = {"default": "{:34s}".format("outer window as in the .win"),
         "cut": "{:34s}".format("dis_win_max = {:.0f} eV (below cut)".format(win_cut))}

print("(A) unitarity of the final U   max_k | U^dag U - 1 |")
for window in ["default", "cut"]:
    for project in [False, True]:
        print("      {}  {}: {:.1e}".format(label[window], "with projection   " if project else "without projection",
                                           unitarity(runs[window, project][0])))
print()

print("(B) point-group residual of H(k) = U^dag diag(eps) U   max | D^dag H(k) D - H(gk) |")
print("    reference: symWannier's rotation matrices D of the projections (prefix.isym)")
print("    max over {} mesh k-points, {} operations and {}x{} matrix elements, in eV:".format(
    w.nk, w.sym.nsym, w.num_wann, w.num_wann))
for window in ["default", "cut"]:
    for project in [False, True]:
        r, k = residual(runs[window, project][0], ham(runs[window, project][0]))
        print("      {}  {}: {:.2e} eV   (worst at k = {})".format(
            label[window], "with projection   " if project else "without projection", r, np.round(k, 4)))
print()

print("(C) the same residual of H(k) rebuilt from the written hr.dat (-H)")
for window in ["default", "cut"]:
    w, hr = runs[window, True]
    print("      {}  with projection   : {:.2e} eV   (max | H_hr.dat - H | = {:.1e})".format(
        label[window], residual(w, hr)[0], np.max(np.abs(hr - ham(w)))))
print()

# U^dag U is the overlap of the Wannier functions at each k; H(k) = U^dag diag(eps) U
# assumes it is 1.  Both checks below vanish for a unitary U.
shift = 10.0
print("(D) other effects of U^dag U != 1 (outer window as in the .win), in eV")
for project in [False, True]:
    w = runs["default", project][0]
    H = ham(w)
    H_shift = np.einsum("kni,kn,knj->kij", np.conj(w.Umat), w.eig.eig - shift, w.Umat) + shift * np.eye(w.num_wann)
    S = np.matmul(np.conj(w.Umat).swapaxes(1, 2), w.Umat)
    zero = np.max(np.abs(np.linalg.eigvalsh(H) - np.linalg.eigvalsh(H_shift)))
    overlap = max(np.max(np.abs(np.linalg.eigvalsh(H[k]) - scipy.linalg.eigh(H[k], S[k], eigvals_only=True)))
                  for k in range(w.nk))
    print("      {}: energy zero moved by {:.0f} eV {:.1e}   eig(H) vs eig(H, S) {:.1e}".format(
        "with projection   " if project else "without projection", shift, zero, overlap))
print()

# where the remainder of (B) lives
sym = runs["default", True][0].sym
for iks, k in enumerate(sym.irr_kpoints):
    d = sym.repmat[iks, sym.little_group(iks)]
    norm = np.min(np.real(np.einsum("hmn,hmn->hn", np.conj(d), d)), axis=0)
    cut = np.flatnonzero(norm < 1 - 1e-6)
    weight = {window: np.sum(np.abs(runs[window, True][0].Umat[sym.iks2ik[iks]][cut]) ** 2)
              for window in ["default", "cut"]}
    if len(cut) and weight["default"] > 0:
        print("Cut multiplet at k = {}: band(s) {} keep {} of their norm under G_k;".format(
            np.round(k, 4), [int(n) + 1 for n in cut], np.round(norm[cut], 4).tolist()))
        print("  weight of the final U on them = {:.1e} (default window), {:.1e} (dis_win_max = {:.0f} eV)".format(
            weight["default"], weight["cut"], win_cut))
