"""Regression test for wannierization on the diamond input.

The test checks the generated output files and the reference spreads and
centers for the four Wannier functions.
"""

import numpy as np
import pytest
from symwannier.wannierize import Wannierize


def test_wannierize_diamond_basic(run_wannier):
    """Check that wannierization of the diamond input converges as expected.

    This verifies creation of `diamond_py_hr.dat` and `diamond_py_tb.dat`, and
    checks the spreads, centers, and total spread against reference values.
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
