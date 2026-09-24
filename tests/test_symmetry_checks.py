import numpy as np
import pytest

from symwannier.nnkp import Nnkp
from symwannier.sym import Sym


def read_sym(test_data_dir, material):
    nnkp = Nnkp(file_nnkp=str(test_data_dir / f"{material}.nnkp"))
    return Sym(file_sym=str(test_data_dir / f"{material}.isym"), nnkp=nnkp)


def test_little_group_sizes(test_data_dir):
    """The little group of a high-symmetry k-point contains the expected operations."""
    sym = read_sym(test_data_dir, "graphene")

    gamma = np.flatnonzero(np.all(np.isclose(sym.irr_kpoints, 0.0), axis=1))[0]
    assert len(sym.little_group(gamma)) == sym.nsym          # G_Gamma is the full group

    kk = np.flatnonzero(np.all(np.isclose(sym.irr_kpoints, [1/3, 1/3, 0.0]), axis=1))[0]
    assert len(sym.little_group(kk)) == 24

    # every operation of G_k really leaves k invariant (up to a reciprocal lattice vector)
    for iks, k in enumerate(sym.irr_kpoints):
        for isym in sym.little_group(iks):
            sk = np.dot(sym.s[isym], k)
            if sym.t_rev[isym] == 1:
                sk = -sk
            assert np.allclose(k - sk, np.round(k - sk))


@pytest.mark.parametrize("material", ["diamond", "H", "Fe_atom_proj"])
def test_repmat_unitary(test_data_dir, material, caplog):
    """Complete multiplets: the representation matrices stay unitary."""
    with caplog.at_level("WARNING"):
        read_sym(test_data_dir, material)
    assert "not unitary" not in caplog.text


@pytest.mark.parametrize(
    "material, kpoint, bands",
    [
        ("graphene", (1/3, 1/3, 0.0), [16]),
        ("Sn", (0.0, 0.0, 0.0), [23, 24]),
        ("Ni_atom_proj", (0.0, 0.0, 0.25), [25]),
    ],
)
def test_repmat_not_unitary_is_reported(test_data_dir, material, kpoint, bands, caplog):
    """A multiplet cut by num_bands makes repmat non-unitary and must be reported.

    The G_k average of Amn/Umat is then not a projector any more, so the affected
    bands are shrunk instead of symmetrized.
    """
    with caplog.at_level("WARNING"):
        sym = read_sym(test_data_dir, material)

    assert "not unitary" in caplog.text
    for n in bands:
        assert f"band {n} keeps only" in caplog.text

    # the reported bands are exactly the ones losing norm under the little group
    iks = np.flatnonzero(np.all(np.isclose(sym.irr_kpoints, kpoint), axis=1))[0]
    d = sym.repmat[iks, sym.little_group(iks), :, :]
    norm = np.real(np.einsum("hmn,hmn->hn", np.conj(d), d))
    assert np.flatnonzero(np.min(norm, axis=0) < 1 - 1e-6).tolist() == [n - 1 for n in bands]


def test_amn_is_checked_not_symmetrized_on_read(test_data_dir, caplog):
    """Reading Amn expands the IBZ data but leaves the values as computed.

    At the K point of graphene the multiplet of band 16 is cut by num_bands, so the
    G_k average is not a projector there: it used to shrink that band by 39%. The
    band must now keep the norm it has in the file.
    """
    from symwannier.amn import Amn

    nnkp = Nnkp(file_nnkp=str(test_data_dir / "graphene.nnkp"))
    sym = Sym(file_sym=str(test_data_dir / "graphene.isym"), nnkp=nnkp)
    with caplog.at_level("WARNING"):
        amn = Amn(file_amn=str(test_data_dir / "graphene.iamn"), nnkp=nnkp, sym=sym)

    iks = np.flatnonzero(np.all(np.isclose(sym.irr_kpoints, [1/3, 1/3, 0.0]), axis=1))[0]
    ik = sym.iks2ik[iks]

    # values from the file; the expansion only mixes the projections, so the norm
    # over the Wannier index is preserved
    assert np.isclose(np.linalg.norm(amn.amn[ik, 15, :]), 0.045659, atol=1e-6)

    # a band with a complete multiplet is unaffected either way
    assert np.isclose(np.linalg.norm(amn.amn[ik, 2, :]), 0.526530, atol=1e-6)

    # the deviation from G_k symmetry is reported
    assert "Amn is not symmetric under G_k" in caplog.text
