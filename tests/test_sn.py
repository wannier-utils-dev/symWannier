import numpy as np
import pytest


def test_wannierize_sn(run_wannier):
    """Run Wannierize on Sn inputs and validate against reference output."""
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
