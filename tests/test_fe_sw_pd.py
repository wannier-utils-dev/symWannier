"""Regression tests for the Fe SW+PD case.

This file exercises the Fe `SW+PD` workflow:
`irr_bz + atom_proj + projectability disentanglement`.

The main checks are:

1. Static input sizes and headers remain intact.
2. The disentanglement windows derived from projectability remain stable.
3. The `iamn` data preserves little-group covariance and yields unitary
   U matrices.
4. The `-P -S` end-to-end run remains numerically healthy.

The goal is to cover regressions in the `atom_proj + irr_bz` path that are
hard to catch with the smaller legacy tests.
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


def _case_prefix(test_data_dir: Path) -> Path:
    """Return the shared prefix for the Fe SW+PD input set."""
    return test_data_dir / "fe_sw_pd"


def _iamn_header(path: Path) -> tuple[int, int, int]:
    """Read `(num_bands, nks, nproj)` from the second line of an `.iamn` file.

    The helper supports both plain-text and `.gz` inputs so the tests keep
    working if the storage format changes.
    """
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as fp:
        fp.readline()
        return tuple(int(x) for x in fp.readline().split())


def _svd_umat(amn_data: np.ndarray) -> np.ndarray:
    """Build the initial U matrix at each k-point from AMN via SVD.

    This mirrors the initialization logic used inside `Amn.Umat()`, but keeps
    it explicit so the tests can inspect each stage separately.
    """
    nk, num_bands, num_wann = amn_data.shape
    umat = np.zeros((nk, num_bands, num_wann), dtype=complex)
    for k in range(nk):
        u, _s, vh = scipy.linalg.svd(amn_data[k, :, :])
        umat[k, :, :] = np.matmul(u[:, :num_wann], vh)
    return umat


def _unitarity_error(umat: np.ndarray) -> tuple[float, float]:
    """Return the unitary error of a U-matrix set as `(max, mean)`.

    For each k-point, this computes the norm of `U^dagger U - I` and is used to
    check that unitarity is preserved through the symmetrization stages.
    """
    errs = []
    for k in range(umat.shape[0]):
        gram = umat[k].conj().T @ umat[k]
        errs.append(np.linalg.norm(gram - np.eye(gram.shape[0])))
    return float(np.max(errs)), float(np.mean(errs))


def _little_group_covariance_stats(prefix: Path) -> tuple[float, float, float]:
    """Return `(mean, p95, max)` residuals for little-group covariance of `iamn`.

    For each irreducible k-point, the test applies the little-group operation
    through `repmat`, the Wannier-side rotation matrix, and the phase factor,
    then measures how closely the transformed AMN returns to the original data.
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


def test_fe_sw_pd_input_metadata(test_data_dir):
    """Lock down the basic metadata of the static Fe SW+PD inputs.

    This test focuses on array sizes and header values rather than physical
    observables. If `nk`, `nks`, `nsym`, `num_bands`, or `num_wann` drift, later
    symmetry and disentanglement failures become much harder to diagnose.
    """
    prefix = _case_prefix(test_data_dir)
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


def test_fe_sw_pd_projectability_windows(copy_inputs, work_dir):
    """Lock down the automatically chosen projectability-based windows.

    The main Python-side change imported for this case sits in
    `dis_window_projectability()`, where the inner and outer windows are chosen
    from the projectability at each k-point. This test therefore checks:

    - projectability contains no NaN, inf, or negative values
    - the maximum projectability stays in the expected range
    - the min, max, and mean of `ndimfroz` and `ndimwin` remain unchanged
    """
    from symwannier.wannierize import Wannierize

    copy_inputs("fe_sw_pd", work_dir)
    prefix = work_dir / "fe_sw_pd"
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
    """Check that raw `iamn` satisfies little-group covariance to high accuracy.

    In this case, `iamn` reflects the symmetry quality quite directly. Keeping
    the mean, p95, and max residuals very small confirms that the
    `atom_proj + irr_bz` inputs are read consistently on the Python side.
    """
    mean_r, p95_r, max_r = _little_group_covariance_stats(_case_prefix(test_data_dir))
    assert mean_r < 1e-6
    assert p95_r < 1e-6
    assert max_r < 1e-6


def test_fe_sw_pd_symmetrized_umat_is_unitary(test_data_dir):
    """Check that U-matrix unitarity survives each SVD and symmetrization step.

    The test records the error after:

    - the raw SVD
    - selecting only irreducible k-points
    - `symmetrize_Gk()`
    - `symmetrize_expand()`

    Keeping these stages separate makes it easier to pinpoint where a future
    regression first breaks unitarity.
    """
    prefix = _case_prefix(test_data_dir)
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
    """Check that the `-P -S` end-to-end run stays numerically healthy.

    Rather than requiring an exact match to an archived log, this slow test
    checks that the current code still satisfies the basic health conditions:

    - `hr/tb` outputs are generated
    - `spreads` and centers are finite
    - the maximum spread is not an obvious outlier
    - the centers do not drift to unphysical values

    This gives CI at least a basic guardrail for the heavy Fe SW+PD path.
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
