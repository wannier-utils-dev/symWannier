import numpy as np
import pytest

from symwannier.nnkp import Nnkp
from symwannier.sym import Sym
from symwannier.amn import Amn
from symwannier.dmn import Dmn


def build(test_data_dir, material):
    nnkp = Nnkp(file_nnkp=str(test_data_dir / f"{material}.nnkp"))
    sym = Sym(file_sym=str(test_data_dir / f"{material}.isym"), nnkp=nnkp)
    amn = Amn(file_amn=str(test_data_dir / f"{material}.iamn"), nnkp=nnkp, sym=sym)
    return Dmn(nnkp=nnkp, sym=sym, amn=amn), sym, amn


@pytest.mark.parametrize("material, nsym, nkirr", [("diamond", 48, 8), ("H", 48, 10), ("Sn", 48, 8)])
def test_dmn_reproduces_amn(test_data_dir, material, nsym, nkirr):
    """The d matrices connect Amn of the full BZ, which is what wannier90 uses them for.

    Sn is a spinor case, where the sign of the double group has to cancel between
    the two matrices instead of being applied to one of them.

    The relation can only hold as well as the Amn written by pw2wannier90 is itself
    symmetric under the little group, which for these inputs is 6e-09 (diamond),
    6e-08 (H) and 1e-07 (Sn); the deviation found here follows those. The threshold
    is the one Dmn.check() warns at, and it still separates these cases from a
    broken one by four orders of magnitude.
    """
    dmn, sym, amn = build(test_data_dir, material)

    assert dmn.nsym == nsym                     # only the spatial operations
    assert dmn.nsym * 2 == sym.nsym             # prefix.isym also holds their time-reversed copies
    assert dmn.nkirr == nkirr
    assert dmn.nk == sym.nkf
    assert dmn.check(amn.amn) < 1e-6


def test_dmn_structure(test_data_dir):
    """kptsym, ir2ik and ik2ir have to agree with the k-point mapping of prefix.isym."""
    dmn, sym, _ = build(test_data_dir, "diamond")

    assert dmn.isym_list[0] == 0 or np.array_equal(sym.s[dmn.isym_list[0]], np.eye(3, dtype=int))
    assert np.all(sym.t_rev[dmn.isym_list] == 0)

    for iks in range(dmn.nkirr):
        assert dmn.ik2ir[dmn.ir2ik[iks]] == iks
        k = sym.irr_kpoints[iks]
        for i, isym in enumerate(dmn.isym_list):
            kdiff = np.dot(sym.s[isym], k) - sym.full_kpoints[dmn.kptsym[i, iks]]
            assert np.allclose(kdiff, np.round(kdiff))
        # the identity leaves the irreducible k-point where it is
        assert dmn.kptsym[0, iks] == dmn.ir2ik[iks]

    # the d matrices are unitary wherever the representation matrices are
    for i in range(dmn.nsym):
        for iks in range(dmn.nkirr):
            d = dmn.d_matrix_wann[i, iks]
            assert np.allclose(d.conj().T @ d, np.eye(dmn.num_wann), atol=1e-8)


def test_dmn_file_round_trip(test_data_dir, tmp_path):
    """The file is written in the order wannier90's list-directed read expects."""
    dmn, _, _ = build(test_data_dir, "diamond")
    path = tmp_path / "diamond.dmn"
    dmn.write_dmn(str(path))

    values = path.read_text().split("\n", 1)[1]
    values = values.replace("(", " ").replace(")", " ").replace(",", " ").split()
    nb, nsym, nkirr, nk = (int(x) for x in values[:4])
    assert (nb, nsym, nkirr, nk) == (dmn.num_bands, dmn.nsym, dmn.nkirr, dmn.nk)

    p = 4
    assert np.array_equal(np.array(values[p:p+nk], dtype=int), dmn.ik2ir + 1)
    p += nk
    assert np.array_equal(np.array(values[p:p+nkirr], dtype=int), dmn.ir2ik + 1)
    p += nkirr
    assert np.array_equal(np.array(values[p:p+nsym*nkirr], dtype=int).reshape(nkirr, nsym).T,
                          dmn.kptsym + 1)
    p += nsym * nkirr

    def read_mat(p, n):
        a = np.array(values[p:p+2*n*n*nsym*nkirr], dtype=float)
        c = (a[0::2] + 1j*a[1::2]).reshape(nkirr, nsym, n, n).transpose(1, 0, 3, 2)
        return c, p + 2*n*n*nsym*nkirr

    dw, p = read_mat(p, dmn.num_wann)
    db, p = read_mat(p, dmn.num_bands)
    assert p == len(values)
    assert np.allclose(dw, dmn.d_matrix_wann, atol=1e-9)
    assert np.allclose(db, dmn.d_matrix_band, atol=1e-9)


def test_dmn_rejects_time_reversal(test_data_dir):
    """Fe_atom_proj needs time reversal to cover the mesh, which the dmn cannot express."""
    with pytest.raises(ValueError, match="time reversal"):
        build(test_data_dir, "Fe_atom_proj")


@pytest.mark.parametrize("material", ["graphene", "Ni_atom_proj"])
def test_dmn_reports_incomplete_multiplets(test_data_dir, material, caplog):
    """A multiplet cut by num_bands makes repmat non-unitary, and the dmn inherits it."""
    dmn, _, amn = build(test_data_dir, material)
    with caplog.at_level("WARNING"):
        assert dmn.check(amn.amn) > 1e-3
    assert "not consistent with Amn" in caplog.text
