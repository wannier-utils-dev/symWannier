"""Fe SW+PD ケース向けの回帰テスト。

このファイルでは、`symwan_proj` から取り込んだ Fe の `SW+PD`
(`irr_bz + atom_proj + projectability disentanglement`) ケースを、
`symWannier` 側で継続的に検証する。

確認したい観点は大きく 4 つある。

1. 静的入力のサイズやヘッダーが壊れていないこと
2. projectability から作る disentanglement window が期待どおりであること
3. `iamn` の little-group 共変性と、そこから作る U 行列のユニタリ性が保たれること
4. `-P -S` の end-to-end 実行が数値的に破綻しないこと

既存の小規模テストでは拾いにくい、`atom_proj + irr_bz` 経路の回帰を
補うのが目的である。
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pytest
import scipy.linalg

from symwannier.amn import Amn
from symwannier.eig import Eig
from symwannier.mmn import Mmn
from symwannier.nnkp import Nnkp
from symwannier.sym import Sym


def _prefix(test_data_dir: Path) -> Path:
    """Fe SW+PD 入力群の共通 prefix を返す。"""
    return test_data_dir / "fe_sw_pd"


def _iamn_header(path: Path) -> tuple[int, int, int]:
    """`.iamn` ヘッダー 2 行目から `(num_bands, nks, nproj)` を読む。

    平文と `.gz` の両方に対応しておき、テストデータの持ち方を変えても
    同じ helper を使えるようにする。
    """
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as fp:
        fp.readline()
        return tuple(int(x) for x in fp.readline().split())


def _svd_umat(amn_data: np.ndarray) -> np.ndarray:
    """AMN から各 k 点の初期 U 行列を SVD で構成する。

    `Amn.Umat()` の内部で使っている初期化と同じ発想を、
    テスト側から段階的に検証できるよう切り出している。
    """
    nk, num_bands, num_wann = amn_data.shape
    umat = np.zeros((nk, num_bands, num_wann), dtype=complex)
    for k in range(nk):
        u, _s, vh = scipy.linalg.svd(amn_data[k, :, :])
        umat[k, :, :] = np.matmul(u[:, :num_wann], vh)
    return umat


def _unitarity_error(umat: np.ndarray) -> tuple[float, float]:
    """U 行列集合のユニタリ誤差を `(max, mean)` で返す。

    各 k 点で `U^dagger U - I` のノルムを計算し、
    対称化の各段階でユニタリ性が崩れていないかを見る。
    """
    errs = []
    for k in range(umat.shape[0]):
        gram = umat[k].conj().T @ umat[k]
        errs.append(np.linalg.norm(gram - np.eye(gram.shape[0])))
    return float(np.max(errs)), float(np.mean(errs))


def _little_group_covariance_stats(prefix: Path) -> tuple[float, float, float]:
    """raw `iamn` の little-group 共変残差を統計量で返す。

    既約 k 点 `k` を little-group で写したときに、
    `repmat`、Wannier 側回転行列、位相因子を通した AMN が元の AMN に
    戻るかを確認する。返り値は `(mean, p95, max)`。
    """
    nnkp = Nnkp(str(prefix) + ".nnkp")
    sym = Sym(str(prefix) + ".isym", nnkp=nnkp)
    amn_raw = Amn(str(prefix) + ".iamn", nnkp=nnkp, sym=None).amn
    amn_ctx = Amn(str(prefix) + ".iamn", nnkp=nnkp, sym=sym)
    rmat, rshift, _ = amn_ctx.projection_sym_mat()

    residuals = []
    for ik, kpt in enumerate(sym.irr_kpoints):
        base = np.linalg.norm(amn_raw[ik])
        for isym, smat in enumerate(sym.s):
            sk = np.dot(smat, kpt)
            if sym.t_rev[isym] == 1:
                sk = -sk
            if not np.allclose(kpt - sk, np.round(kpt - sk)):
                continue

            phase = np.exp(-1j * 2 * np.pi * np.einsum("a,na->n", kpt, rshift[isym]))
            transformed = np.einsum(
                "ml,ln,n->mn", amn_raw[ik], rmat[isym], phase, optimize=True
            )
            if sym.t_rev[isym] == 1:
                transformed = np.conj(transformed)
            transformed = np.einsum(
                "ml,ln->mn", sym.repmat[ik, isym], transformed, optimize=True
            )
            residuals.append(np.linalg.norm(transformed - amn_raw[ik]) / (base + 1e-15))

    residuals = np.asarray(residuals)
    return (
        float(np.mean(residuals)),
        float(np.percentile(residuals, 95)),
        float(np.max(residuals)),
    )


def test_fe_sw_pd_static_inputs(test_data_dir):
    """静的入力の基本メタデータを回帰固定する。

    ここではファイルの物理量そのものではなく、
    取り込みの前提になる配列サイズとヘッダー値を確認する。
    `nk/nks/nsym/num_bands/num_wann` がずれると、
    以降の対称化や disentanglement の失敗原因が切り分けにくくなるため、
    まず最初にこのテストで壊れを止める。
    """
    prefix = _prefix(test_data_dir)
    nnkp = Nnkp(str(prefix) + ".nnkp")
    sym = Sym(str(prefix) + ".isym", nnkp=nnkp)
    amn = Amn(str(prefix) + ".iamn", nnkp=nnkp, sym=sym)
    mmn = Mmn(str(prefix) + ".immn", nnkp=nnkp, sym=sym)
    eig = Eig(str(prefix) + ".ieig", sym=sym)

    assert nnkp.nk == 512
    assert nnkp.nb == 12
    assert nnkp.num_wann == 0
    assert sym.nks == 59
    assert sym.nsym == 16
    assert amn.num_bands == 50
    assert amn.num_wann == 18
    assert amn.amn.shape == (512, 50, 18)
    assert mmn.mmn.shape == (512, 12, 50, 50)
    assert eig.eig.shape == (512, 50)
    assert _iamn_header(prefix.with_suffix(".iamn")) == (50, 59, 18)


def test_fe_sw_pd_projectability_windows(copy_inputs, tmp_path):
    """projectability ベース window の自動決定結果を固定する。

    取り込み対象の Python 側変更の中心は、
    `dis_window_projectability()` が各 k 点の projectability を見て
    `inner/outer` window を決める部分にある。
    そのため、このテストでは

    - projectability に NaN/inf/負値がないこと
    - 最大値が期待したレンジにあること
    - `ndimfroz` と `ndimwin` の min/max/mean が変わっていないこと

    を直接チェックする。
    """
    from symwannier.wannierize import Wannierize

    copy_inputs("fe_sw_pd", tmp_path)
    prefix = tmp_path / "fe_sw_pd"
    wann = Wannierize(
        prefix=str(prefix),
        lsym=True,
        lsite_sym=True,
        projectability_disentangle=True,
    )
    wann.dis_window_projectability()

    assert np.isfinite(wann.projectability).all()
    assert (wann.projectability >= 0).all()
    assert wann.projectability.min() >= 0.0
    assert np.isclose(wann.projectability.max(), 0.9999648562915521, atol=1e-12)
    assert wann.ndimfroz.min() == 10
    assert wann.ndimfroz.max() == 16
    assert np.isclose(np.mean(wann.ndimfroz), 12.671875)
    assert wann.ndimwin.min() == 18
    assert wann.ndimwin.max() == 22
    assert np.isclose(np.mean(wann.ndimwin), 20.12109375)


def test_fe_sw_pd_iamn_little_group_covariance(test_data_dir):
    """raw `iamn` が little-group 共変性を高精度に満たすことを確認する。

    このケースでは `iamn` 自体が対称性情報の質をかなり直接反映する。
    そのため、mean/p95/max の残差を十分小さく抑えることで、
    `atom_proj + irr_bz` 経路で生成された入力が Python 側でも
    整合して読めていることを確認する。
    """
    mean_r, p95_r, max_r = _little_group_covariance_stats(_prefix(test_data_dir))
    assert mean_r < 1e-6
    assert p95_r < 1e-6
    assert max_r < 1e-6


def test_fe_sw_pd_symmetrized_umat_is_unitary(test_data_dir):
    """SVD と対称化の各段階で U 行列のユニタリ性が崩れないことを確認する。

    段階を分けて

    - raw SVD 直後
    - irreducible k 点だけを抜き出した後
    - `symmetrize_Gk()` 後
    - `symmetrize_expand()` 後

    の誤差を見ておくと、将来崩れた場合に
    「どの段階で壊れたか」をこのテストだけで判断しやすい。
    """
    prefix = _prefix(test_data_dir)
    nnkp = Nnkp(str(prefix) + ".nnkp")
    sym = Sym(str(prefix) + ".isym", nnkp=nnkp)
    amn = Amn(str(prefix) + ".iamn", nnkp=nnkp, sym=sym)

    umat_raw = _svd_umat(amn.amn)
    emax_raw, _ = _unitarity_error(umat_raw)
    assert emax_raw < 1e-10

    umat_irk = umat_raw[amn.sym.iks2ik[:]]
    emax_irk, _ = _unitarity_error(umat_irk)
    assert emax_irk < 1e-10

    umat_irk_sym = amn.symmetrize_Gk(umat_irk)
    emax_irk_sym, _ = _unitarity_error(umat_irk_sym)
    assert emax_irk_sym < 1e-10

    umat_full = amn.symmetrize_expand(umat_irk_sym)
    emax_full, _ = _unitarity_error(umat_full)
    assert emax_full < 1e-10


@pytest.mark.slow
def test_fe_sw_pd_full_run_regression(run_wannier):
    """`-P -S` の end-to-end 実行が数値的に健全であることを確認する。

    archived log の厳密一致ではなく、現行コードで再計算した結果が
    少なくとも次を満たすことを slow テストで確認する。

    - `hr/tb` 出力が生成される
    - `spreads` と中心座標が有限値である
    - 最大 spread が明らかな外れ値になっていない
    - 中心座標が不自然に大きくずれていない

    これにより、Fe SW+PD の重い経路を CI でも最低限監視できる。
    """
    wann, workdir = run_wannier(
        "fe_sw_pd",
        lsym=True,
        lsite_sym=True,
        projectability_disentangle=True,
    )

    assert (workdir / "fe_sw_pd_py_hr.dat").exists()
    assert (workdir / "fe_sw_pd_py_tb.dat").exists()
    assert np.isfinite(wann.spreads).all()
    assert np.isfinite(wann.r).all()
    assert np.max(wann.spreads) < 21.0
    assert np.max(np.abs(wann.r)) < 0.1
