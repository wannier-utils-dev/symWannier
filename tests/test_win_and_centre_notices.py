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
        ("guiding_centres = T", True),
        ("guiding_centres = .false.", False),
        ("guiding_centres = .f.", False),
        # wannier90 takes ! and # as comments anywhere in the line, and lets the
        # keyword be indented. Reading the raw line made the first spelling below
        # come out False, the opposite of what it says
        ("guiding_centres = .true. ! on", True),
        ("guiding_centres = .true.  # on", True),
        ("guiding_centres = .false. ! off", False),
        ("    guiding_centres = .false.", False),
        ("\tguiding_centres = .true.", True),
    ],
)
def test_win_reads_guiding_centres(tmp_path, line, expected):
    (tmp_path / "a.win").write_text(
        "num_wann = 4\nnum_bands = 4\nmp_grid : 2 2 2\n" + line + "\n")
    assert Win(str(tmp_path / "a")).guiding_centres is expected


def test_win_rejects_a_value_that_is_not_a_logical(tmp_path):
    """A typo must not quietly become False, which is the branch the default avoids."""
    (tmp_path / "a.win").write_text(
        "num_wann = 4\nnum_bands = 4\nmp_grid : 2 2 2\nguiding_centres = treu\n")
    with pytest.raises(ValueError, match="guiding_centres = treu is not a logical"):
        Win(str(tmp_path / "a"))


@pytest.mark.parametrize(
    "line, attr, expected",
    [
        ("mp_grid : 2 2 2 ! the mesh", "mp_grid", [2, 2, 2]),
        ("num_iter = 50 ! iterations", "num_iter", 50),
        ("  dis_froz_max = 7.5  # the window", "dis_froz_max", 7.5),
    ],
)
def test_win_keywords_survive_a_comment(tmp_path, line, attr, expected):
    """The same defect reached every keyword, loudly for the numeric ones.

    int(" 50 ! iterations") raises, so a win file with a commented num_iter could
    not be read at all.
    """
    body = "num_wann = 4\nnum_bands = 8\nmp_grid : 2 2 2\n"
    if attr == "mp_grid":
        body = "num_wann = 4\nnum_bands = 8\n"
    (tmp_path / "a.win").write_text(body + line + "\n")
    assert np.all(getattr(Win(str(tmp_path / "a")), attr) == expected)


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


def _run_sb(copy_inputs, tmp_path, name, line=None):
    """Wannierize the bundled Sb inputs, optionally with a line added to the win."""
    work = tmp_path / name
    work.mkdir()
    copy_inputs("Sb", work)
    if line:
        win = work / "Sb.win"
        win.write_text(win.read_text() + "\n" + line + "\n")

    cwd = os.getcwd()
    try:
        os.chdir(work)
        # L-BFGS: the test compares with the converged value, which conjugate
        # gradients, the default, approach to 2e-4 in the 200 iterations of the win
        # file (13.155809, as wannier90 does) and reach in 374
        wann = Wannierize(prefix="Sb", lsym=True, log_level="ERROR", minimizer="lbfgs")
        wann.run()
    finally:
        os.chdir(cwd)
    return wann


def test_sb_projections_far_from_the_origin(tmp_path, copy_inputs):
    """Sb at 4x4x4, the case the guided branch exists for.

    Its win file, as cif2qewan writes it, puts the first atom at (0.767,0.767,0.767).
    That is a correct crystallographic coordinate, but it is most of a lattice vector
    from the origin, and on this mesh the phase b.r reaches 0.87 pi. The principal
    branch wraps there, and the four Wannier functions on that site land in the wrong
    cell: 40.6 A^2 each, 169.27 A^2 in total, which is what wannier90 gives for the
    same input without guiding_centres. With the guide the same files give 13.1556
    A^2, as wannier90 does with guiding_centres = .true. (13.15564 after 500
    iterations).

    Only the original structure is bundled. The periodic image is the other cure and
    agrees with the guide to 3e-8 in the total, but it is a different input file
    rather than a different code path, so it would not test anything this does not.
    """
    guided = _run_sb(copy_inputs, tmp_path, "guided")
    principal = _run_sb(copy_inputs, tmp_path, "principal", "guiding_centres = .false.")

    assert np.isclose(np.sum(guided.spreads), 13.15564759, rtol=1e-6)
    assert np.isclose(np.sum(principal.spreads), 169.27045529, rtol=1e-6)

    # the eight functions are four on each of the two sites, so with the branch taken
    # consistently the spreads come in pairs
    order = np.sort(guided.spreads)
    assert np.allclose(order[0::2], order[1::2], atol=1e-6)
    assert np.max(guided.spreads) < 2.0

    # and every center sits on an Sb site. Without the guide the four on the far site
    # are a lattice vector away from either of them
    a = guided.nnkp.a
    sites = np.array([[0.767]*3, [0.233]*3]) @ a
    to_site = lambda r: np.min(np.linalg.norm(r[:,None,:] - sites[None,:,:], axis=2), axis=1)
    assert np.max(to_site(guided.r)) < 1.0
    assert np.sum(to_site(principal.r) > 4.0) == 4
