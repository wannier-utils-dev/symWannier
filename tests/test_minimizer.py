import io
import os
import contextlib

import numpy as np
import pytest

from symwannier.wannierize import Wannierize


def _ready(material, tmp_path, copy_inputs, **kw):
    """A Wannierize at the start of the spread minimization."""
    copy_inputs(material, tmp_path)
    cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        with contextlib.redirect_stdout(io.StringIO()):
            wann = Wannierize(prefix=material, lsym=True, log_level="ERROR", **kw)
            if wann.disentangle:
                wann.dis_window()
                wann.dis_project()
                wann.dis_extract()
            wann.init_Umat_and_Mmn()
            wann._reset_minimizer()
    finally:
        os.chdir(cwd)
    return wann


def _random_tangent(shape, seed):
    rng = np.random.default_rng(seed)
    a = rng.standard_normal(shape) + 1j * rng.standard_normal(shape)
    return (a - np.conj(np.swapaxes(a, 1, 2))) / 2


def test_slope_along_the_descent_direction(tmp_path, copy_inputs):
    """The line search takes the slope of Omega along d as -<dw, d>.

    That holds only if calc_dw returns the gradient with its 1/N, which a finite
    difference checks.
    """
    wann = _ready("diamond", tmp_path, copy_inputs)
    omega0 = wann.calc_omega(wann.mmn)
    dw = wann.calc_dw()
    eps = 1e-4
    fd = (wann._omega_along(dw, eps)[0] - wann._omega_along(dw, -eps)[0]) / (2 * eps)
    assert fd == pytest.approx(-wann._ip(dw, dw), rel=1e-6)
    assert omega0 == pytest.approx(wann.calc_omega(wann.mmn))


@pytest.mark.parametrize("material", ["diamond", "GaAs"])
def test_tangent_symmetrizer_is_an_orthogonal_projection(tmp_path, copy_inputs, material):
    """Idempotent and self-adjoint, so a projected descent direction stays one."""
    wann = _ready(material, tmp_path, copy_inputs)
    shape = wann.Umat.shape[:1] + (wann.num_wann, wann.num_wann)
    a, b = _random_tangent(shape, 0), _random_tangent(shape, 1)
    P = wann._symmetrize_tangent
    pa = P(a)
    norm = lambda x: np.sqrt(wann._ip(x, x))

    assert norm(P(pa) - pa) < 1e-12 * norm(pa)
    assert abs(wann._ip(pa, b) - wann._ip(a, P(b))) < 1e-12 * norm(a) * norm(b)
    assert norm(pa + np.conj(np.swapaxes(pa, 1, 2))) < 1e-12 * norm(pa)
    # a random direction is mostly not symmetric
    assert norm(pa) < 0.5 * norm(a)


@pytest.mark.parametrize("material", ["diamond", "GaAs"])
def test_gradient_at_a_symmetric_point_is_symmetric(tmp_path, copy_inputs, material):
    """The gradient at symmetric rotations is itself symmetric, which the
    projection must leave alone; a wrong transformation rule would not."""
    wann = _ready(material, tmp_path, copy_inputs)
    wann.calc_omega(wann.mmn)
    dw = wann.calc_dw()
    p = wann._symmetrize_tangent(dw)
    assert np.sqrt(wann._ip(p - dw, p - dw)) < 1e-5 * np.sqrt(wann._ip(dw, dw))


@pytest.mark.parametrize("material, omega", [("diamond", 2.33178045), ("GaAs", 7.16057671)])
def test_both_minimizers_reach_the_same_minimum(tmp_path, copy_inputs, material, omega):
    """Steepest descent (5 and 18 iterations here), conjugate gradients (4 and 8)
    and L-BFGS (4 and 9) agree."""
    found = {}
    for minimizer in ("sd", "cg", "lbfgs"):
        work = tmp_path / minimizer
        work.mkdir()
        wann = _ready(material, work, copy_inputs, minimizer=minimizer)
        cwd = os.getcwd()
        try:
            os.chdir(work)
            with contextlib.redirect_stdout(io.StringIO()):
                wann.wannierization()
        finally:
            os.chdir(cwd)
        found[minimizer] = np.sum(wann.spreads)
    for minimizer, value in found.items():
        assert value == pytest.approx(omega, abs=2e-8), minimizer


def test_default_minimizer_and_a_wrong_one(tmp_path, copy_inputs):
    assert _ready("diamond", tmp_path, copy_inputs).minimizer == "cg"
    other = tmp_path / "wrong"
    other.mkdir()
    with pytest.raises(ValueError, match="minimizer must be 'sd', 'cg' or 'lbfgs'"):
        _ready("diamond", other, copy_inputs, minimizer="newton")
