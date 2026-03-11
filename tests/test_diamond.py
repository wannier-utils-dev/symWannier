"""ダイヤモンド構造入力に対する Wannier 化の回帰テスト。

Wannier 化の実行後に生成される出力ファイルと、
4 個の Wannier 関数に対する広がり・中心座標の参照一致を確認する。
"""

import numpy as np
import pytest
from symwannier.wannierize import Wannierize


def test_wannierize_diamond_basic(run_wannier):
    """ダイヤモンド入力での Wannier 化が期待どおり収束することを確認する。

    `diamond_py_hr.dat` と `diamond_py_tb.dat` の生成を確認し、
    spread、中心座標、総 spread が参照値に一致するかを検証する。
    """
    wann, workdir = run_wannier("diamond", lsym=True)

    # Files written
    assert (workdir / "diamond_py_hr.dat").exists()
    assert (workdir / "diamond_py_tb.dat").exists()

    # Quantitative checks (values provided by reference run)
    expected_spread = np.array([0.58294612] * 4)
    expected_centers = np.array(
        [
            [0.0, 0.0, 0.0],
            [-0.806995, 0.806995, 0.0],
            [0.0, 0.806995, 0.806995],
            [-0.806995, 0.0, 0.806995],
        ]
    )

    assert wann.spreads.shape == (4,)
    assert wann.r.shape == (4, 3)
    assert np.allclose(wann.spreads, expected_spread, rtol=1e-6)
    assert np.allclose(wann.r, expected_centers, atol=1e-6)

    omega_tot = np.sum(wann.spreads)
    assert np.isclose(omega_tot, 2.33178448, rtol=1e-6)
