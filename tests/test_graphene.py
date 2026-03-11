"""グラフェン入力に対する Wannier 化の回帰テスト。

出力ファイル生成と 5 個の Wannier 関数の広がり・中心座標を
参照値と比較し、総 spread も含めて再現性を確認する。
"""

import numpy as np
import pytest
from symwannier.wannierize import Wannierize


def test_wannierize_graphene_runs(run_wannier):
    """グラフェン入力での Wannier 化結果が期待値どおりかを検証する。

    出力ファイルの生成、spread 配列と中心座標の形状、
    参照計算に基づく各値および総 spread の一致を確認する。
    """
    wann, workdir = run_wannier("graphene", lsym=True)

    # Output files should be written
    assert (workdir / "graphene_py_hr.dat").exists()
    assert (workdir / "graphene_py_tb.dat").exists()

    # Expected reference values
    expected_spreads = np.array([0.96090325, 0.96090325, 0.57808270, 0.57808270, 0.57808270])
    expected_centers = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 1.406006, 0.0],
        [0.0, 0.703003, 0.0],
        [-0.608818, -0.351501, 0.0],
        [0.608818, -0.351501, 0.0],
    ])

    assert wann.spreads.shape == (5,)
    assert wann.r.shape == (5, 3)
    assert np.allclose(wann.spreads, expected_spreads, rtol=1e-6)
    assert np.allclose(wann.r, expected_centers, atol=1e-5)

    omega_tot = np.sum(wann.spreads)
    assert np.isclose(omega_tot, 3.65605461, rtol=1e-6)
