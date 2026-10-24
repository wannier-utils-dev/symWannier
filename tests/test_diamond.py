import numpy as np
import pytest
from symwannier.wannierize import Wannierize


def test_wannierize_diamond_basic(run_wannier):
    """Run wannierization for diamond inputs and check reported spreads/centers."""
    wann, workdir = run_wannier("diamond", lsym=True)

    # Files written
    assert (workdir / "diamond_py_hr.dat").exists()
    assert (workdir / "diamond_py_tb.dat").exists()

    # Quantitative checks. The minimum, reached in 4 of the 20 iterations the win
    # file allows; wannier90 gives 2.33178042 in total. The earlier steepest descent,
    # whose step was capped at 1, stopped at 0.58294612 each after the 20 and needed
    # 35 to get here.
    expected_spread = np.array([0.58294511] * 4)
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
    assert np.isclose(omega_tot, 2.33178045, rtol=1e-6)

