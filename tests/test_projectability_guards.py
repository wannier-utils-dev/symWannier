import os

import numpy as np
import pytest

from symwannier.nnkp import Nnkp
from symwannier.sym import Sym
from symwannier.amn import Amn
from symwannier.wannierize import Wannierize


def test_projectability_needs_orthonormal_projections(test_data_dir, tmp_path, copy_inputs):
    """p_mk = sum_n |<psi_mk|g_n>|^2 is bounded by 1 only for orthonormal projections.

    The projections written by hand in a win file overlap, and graphene reaches 5.0.
    The thresholds dis_proj_min and dis_proj_max live in [0,1], so everything with
    any weight ends up frozen, the outer window is left with no room, and the
    eigensolver used to fail with an opaque LAPACK message about its index range.
    """
    copy_inputs("graphene", tmp_path)
    cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        wann = Wannierize(prefix="graphene", lsym=True, projectability_disentangle=True,
                          log_level="ERROR")
        assert wann.projectability.max() > 1
        assert not wann.projectability_orthonormal
        with pytest.raises(ValueError, match="orthonormal"):
            wann.run()
    finally:
        os.chdir(cwd)


def test_projectability_accepts_orthonormalized_projectors(test_data_dir, tmp_path, copy_inputs):
    """The atom_proj interface orthonormalizes its projectors, and then p stays in [0,1]."""
    copy_inputs("Fe_atom_proj", tmp_path)
    cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        wann = Wannierize(prefix="Fe_atom_proj", lsym=True, projectability_disentangle=True,
                          log_level="ERROR")
        assert wann.projectability.max() <= 1 + 1e-6
        assert wann.projectability_orthonormal
        wann.dis_window_projectability()          # does not raise
    finally:
        os.chdir(cwd)


def read_amn(test_data_dir, material):
    nnkp = Nnkp(file_nnkp=str(test_data_dir / f"{material}.nnkp"))
    sym = Sym(file_sym=str(test_data_dir / f"{material}.isym"), nnkp=nnkp)
    return Amn(file_amn=str(test_data_dir / f"{material}.iamn"), nnkp=nnkp, sym=sym), sym


def test_projection_sym_mat_rejects_a_missing_rotation(test_data_dir):
    """pw2wannier90 omits a rotation matrix whose entries are all below 1e-10.

    The lattice shift of the projection centers cannot be reconstructed then, and
    the code used to raise IndexError on an empty array.
    """
    amn, sym = read_amn(test_data_dir, "diamond")
    sym.rotmat[3, :, 1] = 0.0

    with pytest.raises(ValueError, match="does not rotate projection"):
        amn.projection_sym_mat()


def test_projection_sym_mat_rejects_inconsistent_centers(test_data_dir):
    """A rotation that mixes projections sitting on different centers has no shift.

    The check used to be an exact float comparison inside an assert, which python -O
    removes.
    """
    amn, sym = read_amn(test_data_dir, "diamond")
    nonzero = np.flatnonzero(sym.rotmat[3, :, 1] != 0)
    other = next(i for i in range(amn.num_wann) if i not in nonzero)
    sym.rotmat[3, other, 1] = 1.0

    with pytest.raises(ValueError, match="different centers"):
        amn.projection_sym_mat()
