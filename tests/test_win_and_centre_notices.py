import os

import numpy as np
import pytest

from symwannier.nnkp import Nnkp
from symwannier.sym import Sym
from symwannier.amn import Amn
from symwannier.win import Win
from symwannier.wannierize import Wannierize


@pytest.mark.parametrize(
    "line, expected",
    [
        ("", None),                                  # absent, and the default is not False
        ("guiding_centres = .true.", True),
        ("guiding_centres = true", True),
        ("guiding_centres = .TRUE.", True),
        ("guiding_centres = .false.", False),
    ],
)
def test_win_reads_guiding_centres(tmp_path, line, expected):
    (tmp_path / "a.win").write_text(
        "num_wann = 4\nnum_bands = 4\nmp_grid : 2 2 2\n" + line + "\n")
    assert Win(str(tmp_path / "a")).guiding_centres is expected


@pytest.mark.parametrize(
    "line, guided",
    [
        ("", True),                                  # on unless the win says otherwise
        ("guiding_centres = .true.", True),
        ("guiding_centres = .false.", False),
    ],
)
def test_guide_centres_follow_the_win(test_data_dir, tmp_path, copy_inputs, line, guided):
    """Only guiding_centres = .false. puts the principal branch back.

    The default differs from wannier90's, whose own default is .false.
    """
    copy_inputs("diamond", tmp_path)
    win = tmp_path / "diamond.win"
    if line:
        win.write_text(win.read_text() + "\n" + line + "\n")

    cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        wann = Wannierize(prefix="diamond", lsym=True, log_level="ERROR")
    finally:
        os.chdir(cwd)

    assert (wann.guide_centres is not None) is guided
    if guided:
        # the guide is the projection centers, in cartesian coordinates
        assert np.allclose(wann.guide_centres,
                           np.array(wann.amn.projection_centers) @ wann.nnkp.a)


def test_imln_takes_the_branch_near_the_guide(test_data_dir, tmp_path, copy_inputs):
    """Im ln is multivalued, and the principal branch is only safe within pi of zero.

    Built here from a phase that the principal branch would wrap but which is within
    pi of the guide, which is the situation a projection written at a periodic image
    far from the origin produces.
    """
    copy_inputs("diamond", tmp_path)
    cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        wann = Wannierize(prefix="diamond", lsym=True, log_level="ERROR")
    finally:
        os.chdir(cwd)

    # diamond's own centers are well inside the safe window, so the guide is pushed
    # out to where a periodic image far from the origin would put it
    wann.guide_centres = wann.guide_centres + 4.0 * wann.nnkp.a[0]
    phase = -np.einsum("ba,na->bn", wann.bvec, wann.guide_centres)
    assert np.max(np.abs(phase)) > np.pi          # the principal branch would wrap
    mnn = np.exp(1j * phase)[np.newaxis, :, :] * np.ones([wann.nk, 1, 1])

    assert np.allclose(wann._imln_mnn(mnn), phase[np.newaxis, :, :], atol=1e-10)
    wann.guide_centres = None
    assert not np.allclose(wann._imln_mnn(mnn), phase[np.newaxis, :, :], atol=1e-10)
    assert np.allclose(wann._imln_mnn(mnn), np.log(mnn).imag)


@pytest.mark.parametrize(
    "material, n_moved, n_centres",
    [
        ("Sn", 0, 1),            # the single center is the origin
        ("GaAs", 36, 1),         # one center, symmorphic group, but it sits at (1/4,1/4,1/4)
        ("diamond", 48, 4),      # several centers and a nonsymmorphic group
    ],
)
def test_lattice_shift_report(test_data_dir, material, n_moved, n_centres, caplog):
    """How often an operation moves a projection by a lattice vector, which costs a phase.

    It is not about how many centers there are: GaAs has one, shared by every
    projection, and a symmorphic group, and still 36 of its 48 operations move it,
    because the center is not a point they leave alone.
    """
    nnkp = Nnkp(file_nnkp=str(test_data_dir / f"{material}.nnkp"))
    sym = Sym(file_sym=str(test_data_dir / f"{material}.isym"), nnkp=nnkp)
    with caplog.at_level("INFO"):
        amn = Amn(file_amn=str(test_data_dir / f"{material}.iamn"), nnkp=nnkp, sym=sym)

    _, Rshift, pos = amn.projection_sym_mat()
    moved = np.sum(np.any(np.any(np.abs(Rshift) > 1e-6, axis=2), axis=1))
    assert moved == n_moved
    assert len(np.unique(np.round(pos, 6), axis=0)) == n_centres

    if n_moved == 0:
        assert "rotation alone" in caplog.text
    else:
        assert f"{n_moved} of {sym.nsym} symmetry operations move" in caplog.text
