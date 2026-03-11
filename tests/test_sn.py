"""Sn 系入力に対する Wannier 化の回帰テスト。

出力ファイルの生成に加え、Wannier 関数の広がりと中心座標、
および総 spread が参照値と一致することを確認する。
"""

import numpy as np
import pytest


def test_wannierize_sn(run_wannier):
    """Sn 入力での Wannier 化結果が参照データと一致することを確認する。

    `Sn_py_hr.dat` と `Sn_py_tb.dat` の生成を確認し、8 個の
    Wannier 関数に対する spread、中心座標、総 spread を検証する。
    """
    wann, workdir = run_wannier("Sn", lsym=True)

    # Output files
    assert (workdir / "Sn_py_hr.dat").exists()
    assert (workdir / "Sn_py_tb.dat").exists()

    # Expected reference values
    expected_spreads = np.array([
        1.68328577,
        1.68328577,
        2.44395825,
        2.44395825,
        2.44402472,
        2.44402472,
        2.44402472,
        2.44402472,
    ])
    expected_centers = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    ])

    assert wann.spreads.shape == (8,)
    assert wann.r.shape == (8, 3)
    assert np.allclose(wann.spreads, expected_spreads, rtol=1e-6)
    assert np.allclose(wann.r, expected_centers, atol=1e-6)

    omega_tot = np.sum(wann.spreads)
    assert np.isclose(omega_tot, 18.03058690, rtol=1e-6)

    # Optional: check Omega components if exposed (not directly exposed by class; rely on logs otherwise)
    # If needed later, we can parse stdout via capsys, but keep it minimal for now.
