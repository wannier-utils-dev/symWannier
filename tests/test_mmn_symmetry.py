import os

import numpy as np
import pytest

from symwannier.nnkp import Nnkp
from symwannier.sym import Sym
from symwannier.mmn import Mmn


def adjoint_deviation(test_data_dir, material):
    """|M(k,b) - M(k+b,-b)^dagger|, as a matrix over the band indices.

    <u_mk|u_n,k+b> is by definition the adjoint of <u_n,k+b|u_mk>, whatever gauge
    the states are in, so the expanded Mmn has to satisfy the identity. Measured on
    the Mmn that pw2wannier90 writes for the full BZ of Fe_atom_proj, it holds to
    1e-12.
    """
    nnkp = Nnkp(file_nnkp=str(test_data_dir / f"{material}.nnkp"))
    sym = Sym(file_sym=str(test_data_dir / f"{material}.isym"), nnkp=nnkp)
    mmn = Mmn(file_mmn=str(test_data_dir / f"{material}.immn"), nnkp=nnkp, sym=sym)

    b = nnkp.bvec_crys
    opp = np.full(len(b), -1, dtype=int)
    for i in range(len(b)):
        for j in range(len(b)):
            if np.allclose(b[i], -b[j]):
                opp[i] = j
                break
    assert np.all(opp >= 0), "the b-vector shell is not symmetric"

    dev = np.zeros([sym.nbnd, sym.nbnd])
    for ik in range(sym.nkf):
        for ib in range(len(b)):
            dev = np.maximum(dev, np.abs(mmn.mmn[ik, ib]
                                         - np.conj(mmn.mmn[mmn.kb2k[ik, ib], opp[ib]]).T))
    return dev, sym.nbnd


@pytest.mark.parametrize(
    "material, n_good, tol",
    [
        ("diamond", 4, 1e-6),          # complete
        ("H", 1, 1e-6),                # complete
        ("GaAs", 4, 1e-5),             # complete; needs time reversal, and its Amn is
                                       # only symmetric to 6e-07 to begin with
        ("graphene", 15, 1e-6),        # band 16 is cut
        ("Sn", 22, 1e-6),              # bands 23, 24 are cut
        ("Ni_atom_proj", 21, 1e-6),    # check_repmat only flags band 25
        ("Fe_atom_proj", 49, 1e-6),    # check_repmat flags nothing at all
    ],
)
def test_expanded_mmn_is_a_set_of_overlaps(test_data_dir, material, n_good, tol):
    """How many of the leading bands the symmetry expansion actually reproduces.

    A band whose degenerate multiplet is cut by num_bands cannot be reproduced: the
    symmetry mixes it with a band that is not in the file. The expansion is exact
    below the first such band and wrong from there on, so the test pins where that
    boundary is. Sym.check_repmat finds the cut multiplets of the little groups of
    the irreducible k-points, which is necessary but not sufficient - it misses
    bands 22-24 of Ni_atom_proj and band 50 of Fe_atom_proj, both of which this
    identity catches.
    """
    dev, num_bands = adjoint_deviation(test_data_dir, material)

    assert dev[:n_good, :n_good].max() < tol
    if n_good < num_bands:
        # and the boundary really is where it is said to be. The step is abrupt for
        # graphene (5.5e-07 -> 4.7e-01), Sn (2.9e-07 -> 6.1e-01) and Fe_atom_proj
        # (4.1e-07 -> 2.6e-01), and gradual for Ni_atom_proj (9.1e-07 -> 1.8e-06)
        assert dev[:n_good+1, :n_good+1].max() > tol


@pytest.mark.parametrize(
    "material, n_min, n_affected, e_max",
    [
        ("diamond", 4, 0, None),
        ("Sn", 22, 12, 19.5350),
        ("Fe_atom_proj", 49, 10, 80.3933),
    ],
)
def test_check_bands(test_data_dir, material, n_min, n_affected, e_max):
    """Mmn.check_bands reports which bands the expansion reproduces, and up to which energy.

    The number of reproduced bands changes from k-point to k-point - a multiplet cut
    at one k-point need not be degenerate at another - so num_bands cannot be chosen
    to avoid it. The energy ceiling can be used instead, as the outer window.
    """
    from symwannier.nnkp import Nnkp
    from symwannier.sym import Sym
    from symwannier.mmn import Mmn
    from symwannier.eig import Eig

    nnkp = Nnkp(file_nnkp=str(test_data_dir / f"{material}.nnkp"))
    sym = Sym(file_sym=str(test_data_dir / f"{material}.isym"), nnkp=nnkp)
    mmn = Mmn(file_mmn=str(test_data_dir / f"{material}.immn"), nnkp=nnkp, sym=sym)
    eig = Eig(str(test_data_dir / f"{material}.ieig"), sym=sym)

    n_bands, found = mmn.check_bands(eig=eig.eig)

    assert n_bands.min() == n_min
    assert np.sum(n_bands < mmn.num_bands) == n_affected
    if e_max is None:
        assert found is None
    else:
        assert found == pytest.approx(e_max, abs=1e-3)

    # it agrees with the deviation computed independently above
    dev, num_bands = adjoint_deviation(test_data_dir, material)
    assert dev[:n_bands.min(), :n_bands.min()].max() < 1e-4


def test_wannierize_warns_when_the_window_is_too_high(test_data_dir, tmp_path, copy_inputs, caplog):
    """Fe_atom_proj asks for dis_win_max = 100, above the 80.4 eV the expansion covers."""
    from symwannier.wannierize import Wannierize

    copy_inputs("Fe_atom_proj", tmp_path)
    cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        with caplog.at_level("WARNING"):
            Wannierize(prefix="Fe_atom_proj", lsym=True, log_level="WARNING")
    finally:
        os.chdir(cwd)

    assert "does not reproduce the highest bands" in caplog.text
    assert "dis_win_max = 100" in caplog.text


def test_wannierize_is_quiet_when_every_band_is_covered(test_data_dir, tmp_path, copy_inputs, caplog):
    from symwannier.wannierize import Wannierize

    copy_inputs("diamond", tmp_path)
    cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        with caplog.at_level("WARNING"):
            Wannierize(prefix="diamond", lsym=True, log_level="WARNING")
    finally:
        os.chdir(cwd)

    assert "dis_win_max" not in caplog.text


def test_check_bands_when_nothing_is_reproduced(test_data_dir):
    """No energy ceiling can be given when some k-point has no reproduced band.

    Reached here by choosing the representative of every k-point among the spatial
    operations instead of taking the first matching one. For Fe_atom_proj that is a
    different choice at 9 of its 14 spatial orbits, and the expansion then fails the
    overlap identity from the first band onwards - the two choices define the states
    at those k-points through an antiunitary and a unitary operation respectively,
    which is not a gauge change.
    """
    from symwannier.nnkp import Nnkp
    from symwannier.sym import Sym
    from symwannier.mmn import Mmn
    from symwannier.eig import Eig

    nnkp = Nnkp(file_nnkp=str(test_data_dir / "Fe_atom_proj.nnkp"))
    sym = Sym(file_sym=str(test_data_dir / "Fe_atom_proj.isym"), nnkp=nnkp)
    for ik, k in enumerate(sym.full_kpoints):
        ks = sym.irr_kpoints[sym.equiv[ik]]
        for isym in range(sym.nsym):
            if sym.t_rev[isym]:
                continue
            kdiff = np.dot(sym.s[isym], ks) - k
            if np.allclose(kdiff, np.round(kdiff)):
                sym.equiv_sym[ik] = isym
                break

    mmn = Mmn(file_mmn=str(test_data_dir / "Fe_atom_proj.immn"), nnkp=nnkp, sym=sym)
    eig = Eig(str(test_data_dir / "Fe_atom_proj.ieig"), sym=sym)

    n_bands, e_max = mmn.check_bands(eig=eig.eig)

    assert np.all(n_bands == 0)
    assert e_max is None          # and it does not raise while looking for one
