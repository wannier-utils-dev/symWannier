"""水素系入力の読み込みと Wannier 化を検証する統合テスト群。

個別ファイルパーサーの挙動、対称性を使った全 BZ への展開、
`Wannierize` の初期化と実行結果を一通り確認する。
"""

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
    """`H.nnkp` の基本情報が正しく読み込まれることを確認する。

    k 点数、近接ベクトル数、Wannier 関数数、配列形状に加え、
    先頭 k 点が Gamma 点であることを検証する。
    """
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
    """`H.ieig` の固有値データが期待する形で読み込まれることを確認する。

    IBZ 上の k 点数、バンド数、固有値配列の形状を検証する。
    """
    eig_file = test_data_dir / "H.ieig"
    eig = Eig(str(eig_file))

    assert eig.nk == 10
    assert eig.num_bands == 1
    assert eig.eig.shape == (10, 1)


def test_mmn_parsing_with_symmetry(test_data_dir):
    """`H.immn` が対称操作を使って全 BZ に展開されることを確認する。

    `H.isym` と `H.nnkp` を併用して Mmn を構築し、展開後の k 点数、
    近接数、行列配列形状、`kb2k` 対応表の形状を検証する。
    """
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
    """`H.iamn` の全 BZ 展開と `Umat()` の性質を確認する。

    対称性込みで展開した Amn の配列形状を確認し、生成した `Umat`
    が各 k 点でユニタリになることを検証する。
    """
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
    """`H.isym` に含まれる対称操作情報の読み込み結果を確認する。

    対称操作数、既約 k 点数、全 k 点数、バンド数、
    および対称行列配列の形状を検証する。
    """
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
    """`H.win` の主要設定値が正しく読み込まれることを確認する。

    Wannier 関数数、反復回数、`mp_grid` 属性の有無を検証する。
    """
    prefix = str(test_data_dir / "H")
    win = Win(prefix)
    
    assert win.num_wann == 1
    assert win.num_iter == 20
    assert hasattr(win, 'mp_grid')


def test_wannierize_initialization_with_symmetry(copy_inputs, tmp_path):
    """対称性を有効にした `Wannierize` 初期化時の内部状態を確認する。

    入力ファイルを一時ディレクトリへコピーして初期化を行い、
    格子点数、近接数、対称性フラグ、内部配列形状を検証する。
    """
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
    """H 入力に対する `Wannierize.run()` の実行結果を確認する。

    出力ファイル生成、spread と中心座標の配列生成、
    単一 Wannier 関数に対する総 spread と中心位置の一致を検証する。
    """
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
