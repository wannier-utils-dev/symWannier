"""The k-point distribution under a real MPI launcher.

These run `symwannier wannierize` under mpirun and compare its output with the same
calculation done serially in this process. They need mpi4py importable by the
interpreter running the tests and an mpirun on PATH, and are skipped otherwise.
"""
import os
import re
import shutil
import subprocess
import sys

import numpy as np
import pytest

from symwannier.wannierize import Wannierize
from symwannier.expand_wannier_inputs import main as expand_main

pytest.importorskip("mpi4py")       # the package only; importing it does not start MPI
# Called by name, not by the path which() finds: Open MPI takes an absolute path to
# mpirun as --prefix, and refuses /bin/mpirun as the unsafe prefix "/"
MPIRUN = next((name for name in ("mpirun", "mpiexec") if shutil.which(name)), None)
pytestmark = pytest.mark.skipif(MPIRUN is None, reason="no mpirun on PATH")

_NUMBER = re.compile(r"-?\d+\.\d+(?:[eE][-+]?\d+)?")


def _numbers(path):
    """Every number in an hr/tb file after its first line, which holds a timestamp."""
    with open(path) as fp:
        fp.readline()
        return np.array([float(x) for x in _NUMBER.findall(fp.read())])


def _mpirun(n, args, cwd):
    cmd = [MPIRUN]
    if "Open MPI" in subprocess.run([MPIRUN, "--version"], capture_output=True, text=True).stdout:
        cmd.append("--oversubscribe")
    cmd += ["-n", str(n), sys.executable, "-c", "from symwannier.cli import main; main()"] + args
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               PYTHONPATH=os.pathsep.join(sys.path))
    r = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    return r.stdout


def _serial(material, cwd, **kw):
    here = os.getcwd()
    try:
        os.chdir(cwd)
        wann = Wannierize(prefix=material, prec=True, log_level="ERROR", **kw)
        wann.run()
    finally:
        os.chdir(here)
    return wann


@pytest.mark.parametrize(
    "material, nproc, flags, kw",
    [
        ("diamond", 3, ["-s"], dict(lsym=True)),                     # 64 k-points, uneven blocks
        ("diamond", 2, ["-S"], dict(lsite_sym=True)),                # site symmetry gathers Umat
        ("Ni_atom_proj", 2, ["-s", "-P"], dict(lsym=True, projectability_disentangle=True)),
        ("diamond", 4, [], dict()),                                  # plain full-BZ files
    ],
)
def test_mpi_matches_serial(tmp_path, copy_inputs, material, nproc, flags, kw):
    """Distributed over k-points, the result is the serial one up to the order of the sums.

    Measured over all bundled systems and 2, 3, 4 and 7 ranks, the spreads agree to
    3e-13 and hr.dat to 3e-12. The bound below leaves room for the summation order
    without hiding a missing or doubled k-point, which would move them by 1/nk.
    """
    ser, par = tmp_path / "serial", tmp_path / "mpi"
    for d in (ser, par):
        d.mkdir()
        copy_inputs(material, d)
        if not flags:
            # the path without symmetry reads full-BZ files, made here by the expander
            here = os.getcwd()
            try:
                os.chdir(d)
                expand_main([material])
            finally:
                os.chdir(here)

    wann = _serial(material, ser, **kw)
    out = _mpirun(nproc, ["wannierize", "-H"] + flags + [material], par)

    # only rank 0 reports
    assert out.count("! Final:") == 1
    omega = float(re.search(r"! Final:\s+omega_tot =\s+(\S+)", out).group(1))
    assert omega == pytest.approx(np.sum(wann.spreads), abs=1e-8)

    for suffix in ("_py_hr.dat", "_py_tb.dat"):
        a = _numbers(ser / (material + suffix))
        b = _numbers(par / (material + suffix))
        assert a.shape == b.shape
        assert np.max(np.abs(a - b)) < 1e-9, suffix


def test_an_error_on_one_rank_stops_every_rank(tmp_path, copy_inputs):
    """An error raised at one k-point is raised only on the rank that holds it.

    Here the outer window keeps seven states at k-point 43 and at least eight
    everywhere else, which is too few for eight Wannier functions, and k-point 43 is
    on rank 1 of 2. Without the abort in main(), rank 0 waits for rank 1 in its next
    exchange for as long as the job runs; measured, it was still waiting when killed
    after 90 s. With it the job ends at once with the error of rank 1.
    """
    copy_inputs("Sb", tmp_path)
    win = tmp_path / "Sb.win"
    lines = [l for l in win.read_text().splitlines() if not l.startswith("dis_froz_max")]
    win.write_text("\n".join(lines) + "\ndis_froz_max = -100\ndis_win_min = -3.75\ndis_win_max = 13.75\n")

    cmd = [MPIRUN]
    if "Open MPI" in subprocess.run([MPIRUN, "--version"], capture_output=True, text=True).stdout:
        cmd.append("--oversubscribe")
    cmd += ["-n", "2", sys.executable, "-c", "from symwannier.cli import main; main()",
            "wannierize", "-s", "Sb"]
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               PYTHONPATH=os.pathsep.join(sys.path))
    r = subprocess.run(cmd, cwd=tmp_path, env=env, capture_output=True, text=True, timeout=120)
    assert r.returncode != 0
    assert "could not be broadcast" in r.stderr
