import os
import subprocess
import sys

import numpy as np
import pytest

from symwannier.kdist import KDist


class _FakeComm:
    """Enough of a communicator for KDist to work out the distribution."""

    def __init__(self, rank, size):
        self.rank, self.size = rank, size

    def Get_rank(self):
        return self.rank

    def Get_size(self):
        return self.size


def test_serial_is_the_identity():
    """Without a communicator every operation must leave the arithmetic untouched."""
    rng = np.random.default_rng(0)
    a = rng.standard_normal((10, 3, 2)) + 1j * rng.standard_normal((10, 3, 2))
    kd = KDist(10)

    assert kd.serial and kd.is_root
    assert (kd.k0, kd.k1, kd.nk_loc) == (0, 10, 10)
    assert kd.local(a) is not None and np.shares_memory(kd.local(a), a)
    assert kd.allgather(a) is a
    x = np.float64(1.25)
    assert kd.sum(x) is x
    idx = rng.integers(0, 10, size=(4, 6))
    assert np.array_equal(kd.plan(idx).fetch(a), a[idx])


@pytest.mark.parametrize("nk, size", [(64, 2), (64, 3), (37, 7), (5, 5)])
def test_blocks_cover_the_mesh(nk, size):
    """Contiguous blocks, sizes differing by at most one, every k-point owned once."""
    blocks = [KDist(nk, _FakeComm(r, size)) for r in range(size)]
    assert blocks[0].k0 == 0 and blocks[-1].k1 == nk
    for a, b in zip(blocks, blocks[1:]):
        assert a.k1 == b.k0
    sizes = [kd.nk_loc for kd in blocks]
    assert max(sizes) - min(sizes) <= 1
    for kd in blocks:
        assert np.array_equal(np.flatnonzero(kd.owner == kd.rank), np.arange(kd.k0, kd.k1))


def test_more_ranks_than_k_points():
    with pytest.raises(ValueError, match="at most one rank per k-point"):
        KDist(4, _FakeComm(0, 5))


def _run(code, env):
    return subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True)


def _clean_env(**extra):
    env = {k: v for k, v in os.environ.items()
           if not k.startswith(("OMPI_", "PMI_", "PMIX_", "MV2_", "SLURM_"))}
    env["PYTHONPATH"] = os.pathsep.join(sys.path)
    env.update(extra)
    return env


# mpi4py made unimportable, or replaced by one whose COMM_WORLD has the given size,
# so that every branch runs without MPI
_NO_MPI4PY = "import sys; sys.modules['mpi4py'] = None\n"
_FAKE_MPI4PY = (
    "import sys, types\n"
    "MPI = types.ModuleType('mpi4py.MPI')\n"
    "MPI.COMM_WORLD = type('Comm', (), {{'Get_size': lambda self: {size}}})()\n"
    "pkg = types.ModuleType('mpi4py'); pkg.MPI = MPI\n"
    "sys.modules['mpi4py'] = pkg; sys.modules['mpi4py.MPI'] = MPI\n")
_CALL = ("from symwannier.kdist import world_comm\n"
         "try:\n    c = world_comm()\n    print('comm', None if c is None else c.Get_size())\n"
         "except RuntimeError as e:\n    print('raised', e)\n"
         "print('imported' if sys.modules.get('mpi4py') not in (None,) else 'not imported')\n")


def _world_comm(prelude, **env):
    r = _run("import sys\n" + prelude + _CALL, _clean_env(**env))
    assert r.returncode == 0, r.stderr
    return r.stdout


def test_world_comm_leaves_mpi_alone_outside_a_launcher():
    """A serial run must not import mpi4py, which would initialize MPI."""
    out = _world_comm("")
    assert "comm None" in out and "not imported" in out


def test_world_comm_ignores_the_scheduler():
    """Only MPI's own variables count; a scheduler's describe the job, not the processes.

    These are what a serial command in a Slurm batch script with --ntasks=4 sees, and
    what srun sets without an MPI plugin. Neither starts an MPI process.
    """
    out = _world_comm("", SLURM_PROCID="0", SLURM_NTASKS="4", SLURM_STEP_NUM_TASKS="4")
    assert "comm None" in out and "not imported" in out


@pytest.mark.parametrize(
    "env, expected",
    [
        (dict(OMPI_COMM_WORLD_RANK="0", OMPI_COMM_WORLD_SIZE="4"), "as 4 processes"),
        (dict(PMI_RANK="0", PMI_SIZE="2"), "as 2 processes"),
        (dict(PMIX_RANK="2"), "as process 2"),       # PMIx reports the rank only
    ],
)
def test_world_comm_refuses_several_processes_without_mpi4py(env, expected):
    """Every process would run the whole calculation and write the same files."""
    out = _world_comm(_NO_MPI4PY, **env)
    assert "raised started by an MPI launcher " + expected + ", but mpi4py is not installed" in out


def test_world_comm_process_0_of_an_unknown_count_runs_serially_without_mpi4py():
    """PMIx rank 0 cannot tell whether it is alone; the others refuse, so it writes alone."""
    out = _world_comm(_NO_MPI4PY, PMIX_RANK="0")
    assert "comm None" in out


def test_world_comm_refuses_a_communicator_of_one_for_several_processes():
    """MPI that does not belong to the launcher starts each process as a run of its own."""
    out = _world_comm(_FAKE_MPI4PY.format(size=1), OMPI_COMM_WORLD_RANK="0", OMPI_COMM_WORLD_SIZE="4")
    assert "raised started by an MPI launcher as 4 processes, but MPI sees this process as a run of its own" in out


def test_world_comm_returns_the_communicator():
    out = _world_comm(_FAKE_MPI4PY.format(size=4), OMPI_COMM_WORLD_RANK="0", OMPI_COMM_WORLD_SIZE="4")
    assert "comm 4" in out
