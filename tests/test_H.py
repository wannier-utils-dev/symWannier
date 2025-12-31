import numpy as np
import pytest

from symwannier.nnkp import Nnkp
from symwannier.mmn import Mmn
from symwannier.amn import Amn
from symwannier.eig import Eig
from symwannier.sym import Sym
from symwannier.win import Win
from symwannier.wannierize import Wannierize


def test_nnkp_parsing(test_data_dir):
    """Test parsing of H.nnkp file."""
    nnkp_file = test_data_dir / "H.nnkp"
    nnkp = Nnkp(str(nnkp_file))

    assert nnkp.nk == 64  # 4x4x4 k-point grid
    assert nnkp.nb == 6
    assert nnkp.num_wann == 1
    assert nnkp.bvec.shape == (6, 3)
    assert nnkp.kpoints.shape == (64, 3)
    # Check first k-point is Gamma
    assert np.allclose(nnkp.kpoints[0], [0, 0, 0])


def test_eig_parsing(test_data_dir):
    """Test parsing of H.ieig file."""
    eig_file = test_data_dir / "H.ieig"
    eig = Eig(str(eig_file))

    assert eig.nk == 10
    assert eig.num_bands == 1
    assert eig.eig.shape == (10, 1)


def test_mmn_parsing_with_symmetry(test_data_dir):
    """Test parsing of H.immn file (IBZ Mmn)."""
    mmn_file = test_data_dir / "H.immn"
    isym_file = test_data_dir / "H.isym"
    nnkp_file = test_data_dir / "H.nnkp"
    
    nnkp = Nnkp(str(nnkp_file))
    sym = Sym(file_sym=str(isym_file), nnkp=nnkp)
    mmn = Mmn(str(mmn_file), nnkp=nnkp, sym=sym)

    # Check that Mmn was expanded to full BZ
    assert mmn.nk == 64  # 4x4x4 k-point grid
    assert mmn.num_bands == 1
    assert mmn.nb == 6
    assert mmn.mmn.shape == (64, 6, 1, 1)
    # Check kb2k mapping
    assert mmn.kb2k.shape == (64, 6)


def test_amn_parsing_with_symmetry(test_data_dir):
    """Test parsing of H.iamn file and Umat generation."""
    amn_file = test_data_dir / "H.iamn"
    isym_file = test_data_dir / "H.isym"
    nnkp_file = test_data_dir / "H.nnkp"
    
    nnkp = Nnkp(str(nnkp_file))
    sym = Sym(file_sym=str(isym_file), nnkp=nnkp)
    amn = Amn(str(amn_file), nnkp=nnkp, sym=sym)

    # Check that Amn was expanded to full BZ
    assert amn.nk == 64  # 4x4x4 k-point grid
    assert amn.num_bands == 1
    assert amn.num_wann == 1
    assert amn.amn.shape == (64, 1, 1)
    
    # Generate Umat and check unitarity
    umat = amn.Umat()
    assert umat.shape == (64, 1, 1)
    for k in range(64):
        assert np.allclose(np.conj(umat[k].T) @ umat[k], np.eye(1), atol=1e-5)


def test_sym_parsing(test_data_dir):
    """Test parsing of H.isym symmetry file."""
    isym_file = test_data_dir / "H.isym"
    nnkp_file = test_data_dir / "H.nnkp"
    nnkp = Nnkp(str(nnkp_file))
    sym = Sym(file_sym=str(isym_file), nnkp=nnkp)

    # Check number of symmetry operations (96 for cubic)
    assert sym.nsym == 96
    assert sym.nks == 10  # irreducible k-points
    assert sym.nkf == 64  # full k-points (4x4x4 grid)
    assert sym.nbnd == 1
    
    # Check symmetry matrices
    assert sym.s.shape == (96, 3, 3)


def test_win_parsing(test_data_dir):
    """Test parsing of H.win file."""
    prefix = str(test_data_dir / "H")
    win = Win(prefix)
    
    assert win.num_wann == 1
    assert win.num_iter == 20
    assert hasattr(win, 'mp_grid')


def test_wannierize_initialization_with_symmetry(copy_inputs, tmp_path):
    """Test Wannierize class initialization with symmetry."""
    copy_inputs("H", tmp_path)

    cwd = __import__("os").getcwd()
    try:
        __import__("os").chdir(tmp_path)
        wann = Wannierize(prefix="H", lsym=True, lsite_sym=False)

        # Check initialization
        assert wann.num_wann == 1
        assert wann.num_bands == 1
        assert wann.nk == 64
        assert wann.nb == 6
        assert wann.lsym is True
        assert wann.lsite_sym is False
        assert wann.sym is not None
        assert wann.mmn0.shape == (64, 6, 1, 1)
        assert wann.amn.amn.shape == (64, 1, 1)

    finally:
        __import__("os").chdir(cwd)


def test_wannierize_run_basic(run_wannier):
    """Test basic Wannierize.run() execution for H inputs."""
    wann, workdir = run_wannier("H", lsym=True, num_iter=2)

    # Check output files were created
    assert (workdir / "H_py_hr.dat").exists()
    assert (workdir / "H_py_tb.dat").exists()

    # Check spreads and centers were calculated
    assert hasattr(wann, 'spreads')
    assert hasattr(wann, 'r')
    assert wann.spreads.shape == (1,)
    assert wann.r.shape == (1, 3)

    # Check specific values from the wannierization
    omega_tot = np.sum(wann.spreads)
    assert np.isclose(omega_tot, 1.01101519, rtol=1e-6)
    assert np.allclose(wann.r[0], [0, 0, 0], atol=1e-5)
    assert np.isclose(wann.spreads[0], 1.01101519, rtol=1e-6)
