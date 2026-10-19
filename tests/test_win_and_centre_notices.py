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
        ("", False),
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


def test_wannierize_warns_about_guiding_centres(test_data_dir, tmp_path, copy_inputs, caplog):
    """wannier90 uses it to re-choose the branch of Im ln M_nn; this code does not.

    Silently ignoring it would make the two minimize different things while looking
    comparable.
    """
    copy_inputs("diamond", tmp_path)
    win = tmp_path / "diamond.win"
    win.write_text(win.read_text() + "\nguiding_centres = .true.\n")

    cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        with caplog.at_level("WARNING"):
            Wannierize(prefix="diamond", lsym=True, log_level="WARNING")
    finally:
        os.chdir(cwd)

    assert "guiding_centres" in caplog.text
    assert "no effect here" in caplog.text


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
