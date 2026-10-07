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
