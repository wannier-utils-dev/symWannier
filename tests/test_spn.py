import numpy as np
import pytest

from symwannier.nnkp import Nnkp
from symwannier.sym import Sym
from symwannier.eig import Eig
from symwannier.spn import Spn


def read_spn(test_data_dir, material, expand):
    nnkp = Nnkp(file_nnkp=str(test_data_dir / f"{material}.nnkp"))
    sym = Sym(file_sym=str(test_data_dir / f"{material}.isym"), nnkp=nnkp)
    spn = Spn(file_spn=str(test_data_dir / f"{material}.ispn"), nnkp=nnkp,
              sym=sym if expand else None)
    return spn, sym, nnkp


def test_spn_is_hermitian_and_kramers(test_data_dir):
    """The file packs the upper triangle column by column, not row by row.

    Reading it the other way round still gives a hermitian matrix, so the test that
    catches it is the Kramers one: with time reversal and spin-orbit coupling the
    bands come in degenerate pairs, and the spin operator traced over such a pair
    has to vanish. The multiplets that reach the top of the window are left out,
    since num_bands cuts one there and the pair is incomplete.
    """
    spn, sym, _ = read_spn(test_data_dir, "Sn", expand=False)
    eig = Eig(str(test_data_dir / "Sn.ieig"))

    assert spn.num_bands == 24
    assert spn.nk == sym.nks
    assert np.allclose(spn.spn, np.conj(spn.spn).transpose(0, 1, 3, 2), atol=1e-12)

    worst = 0.0
    for ik in range(spn.nk):
        e = eig.eig[ik]
        i = 0
        while i < spn.num_bands:
            j = i
            while j + 1 < spn.num_bands and abs(e[j+1] - e[i]) < 1e-4:
                j += 1
            if j < 22:
                worst = max(worst, np.max(np.abs(np.einsum("amm->a", spn.spn[ik][:, i:j+1, i:j+1]).real)))
            i = j + 1
    assert worst < 1e-6


def test_spin_rotation_is_a_rotation(test_data_dir):
    """The spinor matrix of every operation induces an orthogonal rotation of sigma."""
    spn, sym, _ = read_spn(test_data_dir, "Sn", expand=False)
    spn.sym = sym

    for isym in range(sym.nsym):
        rot = spn.spin_rotation(isym)
        assert np.allclose(rot @ rot.T, np.eye(3), atol=1e-10)
        assert np.isclose(abs(np.linalg.det(rot)), 1.0, atol=1e-10)


def test_pure_time_reversal_flips_the_spin(test_data_dir):
    """For the operation that is time reversal alone, the matrix has to be -1.

    It is the statement that time reversal flips the spin, and it pins the sign
    convention of the antiunitary case: the rotation of the spinor matrix -i sigma_y
    is diag(-1,+1,-1) on its own, and only the complex conjugation of the operator,
    which flips sigma_y, turns it into -1 on every component.
    """
    spn, sym, _ = read_spn(test_data_dir, "Sn", expand=False)
    spn.sym = sym

    pure = [isym for isym in range(sym.nsym)
            if np.array_equal(sym.s[isym], np.eye(3, dtype=int))
            and np.allclose(sym.ft[isym], 0) and sym.t_rev[isym] == 1]
    assert len(pure) == 1
    assert np.allclose(spn.spin_rotation(pure[0]), -np.eye(3), atol=1e-10)


def test_time_reversed_representative_gives_the_same_spn(test_data_dir):
    """Which operation represents a k-point is a gauge choice, and Spn has to agree.

    The bundled data reaches every k-point without time reversal, so the antiunitary
    branch is never taken on its own. Forcing it at every k-point has to leave the
    trace, which does not depend on the band basis, where it was.
    """
    nnkp = Nnkp(file_nnkp=str(test_data_dir / "Sn.nnkp"))
    sym = Sym(file_sym=str(test_data_dir / "Sn.isym"), nnkp=nnkp)
    spn = Spn(file_spn=str(test_data_dir / "Sn.ispn"), nnkp=nnkp, sym=sym)

    alt = Sym(file_sym=str(test_data_dir / "Sn.isym"), nnkp=nnkp)
    switched = 0
    for ik, k in enumerate(alt.full_kpoints):
        ks = alt.irr_kpoints[alt.equiv[ik]]
        for isym in range(alt.nsym):
            if alt.t_rev[isym] == 0:
                continue
            kdiff = -np.dot(alt.s[isym], ks) - k
            if np.allclose(kdiff, np.round(kdiff)):
                alt.equiv_sym[ik] = isym
                switched += 1
                break
    assert switched == alt.nkf
    spn_alt = Spn(file_spn=str(test_data_dir / "Sn.ispn"), nnkp=nnkp, sym=alt)

    n_gapped = 8
    tr = np.einsum("kamm->ka", spn.spn[:, :, :n_gapped, :n_gapped]).real
    tr_alt = np.einsum("kamm->ka", spn_alt.spn[:, :, :n_gapped, :n_gapped]).real
    assert np.allclose(tr, tr_alt, atol=1e-9)


def test_spn_expansion(test_data_dir):
    """Expanding to the full BZ rotates the spin operator and leaves the bands alone."""
    raw, sym, nnkp = read_spn(test_data_dir, "Sn", expand=False)
    spn, _, _ = read_spn(test_data_dir, "Sn", expand=True)

    assert spn.nk == sym.nkf > sym.nks
    assert np.allclose(spn.spn, np.conj(spn.spn).transpose(0, 1, 3, 2), atol=1e-12)

    # the irreducible k-points keep the data they were read with
    assert np.array_equal(spn.spn[sym.iks2ik], raw.spn)

    # and the trace, which does not depend on the band basis, follows the rotation.
    # Only the lowest 8 bands are used: they are separated from the rest by 2 eV at
    # every k-point, so they form a subspace the symmetry maps onto itself. Taking
    # all 24 would cut a degenerate multiplet at the top, which is the defect
    # Sym.check_repmat reports for this system, and the trace is then not symmetric
    # in the data to begin with.
    n_gapped = 8
    tr = np.einsum("kamm->ka", spn.spn[:, :, :n_gapped, :n_gapped]).real
    for isym in range(sym.nsym):
        rot = spn.spin_rotation(isym)
        sign = -1 if sym.t_rev[isym] == 1 else 1
        for ik in range(sym.nkf):
            jk = sym.search_ik_full(np.dot(sym.s[isym], nnkp.kpoints[ik]) * sign)
            assert np.allclose(tr[jk], sign * (rot @ tr[ik]), atol=1e-9)


def test_spn_round_trip(test_data_dir, tmp_path):
    """What write_spn produces reads back as the same matrices."""
    spn, _, nnkp = read_spn(test_data_dir, "Sn", expand=True)
    path = tmp_path / "Sn.spn"
    spn.write_spn(str(path))

    again = Spn(file_spn=str(path), nnkp=nnkp)
    assert again.num_bands == spn.num_bands
    assert again.nk == spn.nk
    assert np.allclose(again.spn, spn.spn, atol=1e-9)


def test_spn_rejects_unformatted(test_data_dir, tmp_path):
    """A file written without spn_formatted cannot be read, and says so."""
    path = tmp_path / "bad.ispn"
    path.write_text("header\n4 2\n1.0 0.0\n")
    with pytest.raises(ValueError, match="spn_formatted"):
        Spn(file_spn=str(path), nnkp=None)
