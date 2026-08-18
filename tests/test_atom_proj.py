import numpy as np
import pytest


@pytest.mark.parametrize(
    (
        "material",
        "num_bands",
        "num_wann",
        "num_irr_kpoints",
        "num_symmetries",
        "spinors",
        "expected_total_spread",
    ),
    [
        ("Fe_atom_proj", 50, 18, 13, 16, True, 14.94006918),
        ("Ni_atom_proj", 25, 9, 8, 96, False, 6.73101300),
    ],
)
def test_atom_proj_ibz_wannierize(
    run_wannier,
    material,
    num_bands,
    num_wann,
    num_irr_kpoints,
    num_symmetries,
    spinors,
    expected_total_spread,
):
    """Expand QE atom_proj IBZ data and reproduce the reference spread."""
    wann, workdir = run_wannier(material, lsym=True)

    assert wann.num_bands == num_bands
    assert wann.num_wann == num_wann
    assert wann.nk == 64
    assert wann.amn.amn.shape == (64, num_bands, num_wann)

    assert wann.sym.nks == num_irr_kpoints
    assert wann.sym.nsym == num_symmetries
    assert wann.sym.spinors is spinors
    assert wann.sym.centers.shape == (num_wann, 3)
    assert np.allclose(wann.sym.centers, 0.0, atol=1e-12)

    assert (workdir / f"{material}_py_hr.dat").exists()
    assert (workdir / f"{material}_py_tb.dat").exists()
    assert np.isclose(
        np.sum(wann.spreads),
        expected_total_spread,
        rtol=1e-7,
        atol=2e-6,
    )


def test_fe_atom_proj_ibz_projectability(run_wannier):
    """Use the same energy/projectability window union as Wannier90."""
    wann, workdir = run_wannier(
        "Fe_atom_proj",
        lsym=True,
        projectability_disentangle=True,
    )

    assert wann.win.dis_proj_min == pytest.approx(0.01)
    assert wann.win.dis_proj_max == pytest.approx(0.90)
    assert wann.win.has_dis_froz_window

    assert (np.min(wann.ndimfroz), np.max(wann.ndimfroz)) == (12, 16)
    assert (np.min(wann.ndimwin), np.max(wann.ndimwin)) == (18, 28)

    assert (workdir / "Fe_atom_proj_py_hr.dat").exists()
    assert (workdir / "Fe_atom_proj_py_tb.dat").exists()
    assert np.isclose(
        np.sum(wann.spreads),
        17.33392368582627,
        rtol=1e-7,
        atol=2e-6,
    )
