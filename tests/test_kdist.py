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


def _clean_env():
    env = {k: v for k, v in os.environ.items()
           if not k.startswith(("OMPI_", "PMI_", "PMIX_", "MV2_", "SLURM_"))}
    env["PYTHONPATH"] = os.pathsep.join(sys.path)
    return env


def test_world_comm_leaves_mpi_alone_outside_a_launcher():
    """A serial run must not import mpi4py, which would initialize MPI."""
    r = _run("import sys; from symwannier.kdist import world_comm; "
             "assert world_comm() is None; assert 'mpi4py' not in sys.modules", _clean_env())
    assert r.returncode == 0, r.stderr


def test_world_comm_refuses_several_ranks_without_mpi4py():
    """Every rank would run the whole calculation and write the same files."""
    env = _clean_env()
    env.update(OMPI_COMM_WORLD_RANK="0", OMPI_COMM_WORLD_SIZE="4")
    # make mpi4py unimportable, whether or not it is installed
    r = _run("import sys; sys.modules['mpi4py'] = None\n"
             "from symwannier.kdist import world_comm\n"
             "try:\n    world_comm()\nexcept RuntimeError as e:\n    print(e)\n", env)
    assert r.returncode == 0, r.stderr
    assert "started on 4 MPI ranks, but mpi4py is not installed" in r.stdout
