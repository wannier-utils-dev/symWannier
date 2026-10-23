#!/usr/bin/env python

import os
import numpy as np

# Set for every process an MPI launcher starts - mpirun, or a scheduler's launcher
# with an MPI plugin, such as srun --mpi=pmix. They belong to MPI's process
# management, not to the scheduler: the scheduler's own variables differ from one
# system to the next and describe the job rather than the processes MPI started.
_MPI_RANK_VARS = ("OMPI_COMM_WORLD_RANK", "PMI_RANK", "PMIX_RANK", "MV2_COMM_WORLD_RANK")
_MPI_SIZE_VARS = ("OMPI_COMM_WORLD_SIZE", "PMI_SIZE", "MV2_COMM_WORLD_SIZE")


def _int_env(names):
    """The first of the variables that is set to an integer, or None."""
    for v in names:
        if os.environ.get(v, "").isdigit():
            return int(os.environ[v])
    return None


def world_comm():
    """MPI.COMM_WORLD when the process was started by an MPI launcher as one of several.

    Returns None otherwise, without importing mpi4py outside an MPI launcher:
    importing it initializes MPI, which a serial run has no reason to do. Started
    as one of several processes, every process would run the whole calculation and
    write the same files unless they share a communicator, so it raises when they
    cannot - when mpi4py is missing, and when MPI sees the process as a run of its
    own. Several processes are recognized from the size the launcher reports, or,
    where it reports only the rank (PMIx), from a rank above 0.

    A launcher that starts the processes without telling MPI - a scheduler's launcher
    with no MPI plugin - leaves nothing to recognize them by, and each runs serially.
    """
    if not any(v in os.environ for v in _MPI_RANK_VARS):
        return None
    size = _int_env(_MPI_SIZE_VARS)
    rank = _int_env(_MPI_RANK_VARS)
    several = (size is not None and size > 1) or (rank is not None and rank > 0)
    what = "{} processes".format(size) if size is not None else "process {}".format(rank)
    try:
        from mpi4py import MPI
    except ImportError:
        if several:
            raise RuntimeError(
                "started by an MPI launcher as {}, but mpi4py is not installed; install "
                "it with pip install mpi4py, built against the MPI of the launcher, or "
                "run as a single process".format(what))
        return None
    comm = MPI.COMM_WORLD
    if several and comm.Get_size() == 1:
        raise RuntimeError(
            "started by an MPI launcher as {}, but MPI sees this process as a run of "
            "its own, so every process would compute everything and write the same "
            "files; mpi4py and the launcher have to belong to the same MPI".format(what))
    return comm if comm.Get_size() > 1 else None


class KDist:
    """Distribution of the k-points of the full mesh over MPI ranks.

    Every rank owns one contiguous block [k0, k1) of the full k list. In the order
    wannier90 lists the mesh, with the last index fastest, a contiguous block is a
    slab of the mesh, so most of the neighbours k+b of a k-point are in the same
    block and the rest are in the two adjacent ones.

    Without a communicator - mpi4py missing, or a single rank - the rank owns every
    k-point and every operation below is the identity, so the serial path does the
    same arithmetic as before the distribution was introduced, bit for bit.

    Sums over ranks are reduced on rank 0 and broadcast rather than all-reduced. That
    makes the result the same on every rank to the last bit, which the iterations
    rely on: they decide when to stop, and in the line search where to evaluate
    next, from these sums, and ranks that disagreed by one ulp could take different
    branches and deadlock.
    """

    def __init__(self, nk, comm=None):
        self.nk = nk
        self.comm = comm
        self.rank = 0 if comm is None else comm.Get_rank()
        self.size = 1 if comm is None else comm.Get_size()
        if self.size > nk:
            raise ValueError("{} MPI ranks for {} k-points; use at most one rank per k-point"
                             .format(self.size, nk))
        self.counts = np.array([ nk // self.size + (1 if r < nk % self.size else 0)
                                 for r in range(self.size) ], dtype=int)
        self.offsets = np.concatenate([[0], np.cumsum(self.counts)[:-1]]).astype(int)
        self.k0 = int(self.offsets[self.rank])
        self.k1 = self.k0 + int(self.counts[self.rank])
        self.nk_loc = self.k1 - self.k0
        self.owner = np.repeat(np.arange(self.size), self.counts)   # rank owning each k
        self.is_root = self.rank == 0

    @property
    def serial(self):
        return self.size == 1

    def local(self, a):
        """The block of a full-mesh array [nk, ...] that this rank owns."""
        return a[self.k0:self.k1]

    def is_local(self, ik):
        return self.k0 <= ik < self.k1

    def sum(self, x):
        """Sum of x over the ranks, identical on every rank. x is a scalar or an array."""
        if self.serial:
            return x
        from mpi4py import MPI
        a = np.asarray(x)
        send = np.ascontiguousarray(a).reshape(-1)
        recv = np.empty_like(send)
        self.comm.Reduce(send, recv, op=MPI.SUM, root=0)
        self.comm.Bcast(recv, root=0)
        if a.ndim == 0:
            return recv[0]
        return recv.reshape(a.shape)

    def allgather(self, a_loc):
        """The full-mesh array [nk, ...] assembled from every rank's block."""
        if self.serial:
            return a_loc
        a_loc = np.ascontiguousarray(a_loc)
        out = np.empty((self.nk,) + a_loc.shape[1:], dtype=a_loc.dtype)
        rowtype = self._rowtype(a_loc)
        try:
            self.comm.Allgatherv([a_loc, a_loc.shape[0], rowtype],
                                 [out, (self.counts, self.offsets), rowtype])
        finally:
            rowtype.Free()
        return out

    def plan(self, idx):
        """Prepare to fetch, from the distributed array, the rows at the global k indices idx.

        idx can have any shape; the fetched array has shape idx.shape + row shape. The
        communication pattern is worked out once here, so a plan made for kb2k can be
        used every time the rotations change.
        """
        return _Fetch(self, idx)

    @staticmethod
    def _rowtype(a):
        """One row a[i] as a single MPI datatype, so that counts stay counts of rows."""
        from mpi4py.util.dtlib import from_numpy_dtype
        rowsize = int(np.prod(a.shape[1:], dtype=int))
        return from_numpy_dtype(a.dtype).Create_contiguous(rowsize).Commit()


class _Fetch:
    """Rows of a k-distributed array at a fixed set of global k indices.

    The indices asked for repeat - kb2k names every neighbour once per k-point
    it neighbours - and many are local. Only the distinct remote ones travel, once
    each; the requested layout is assembled from those and the local block here.
    """

    def __init__(self, kd, idx):
        self.kd = kd
        idx = np.asarray(idx, dtype=int)
        self.shape = idx.shape
        self.idx = idx.reshape(-1)
        if kd.serial:
            return

        # distinct k-points wanted, and where each requested entry finds its row
        uniq, self.inverse = np.unique(self.idx, return_inverse=True)
        self.nuniq = len(uniq)
        local = (uniq >= kd.k0) & (uniq < kd.k1)
        self.uniq_local = np.flatnonzero(local)
        self.local_rows = uniq[local] - kd.k0
        self.uniq_remote = np.flatnonzero(~local)
        remote = uniq[~local]                     # sorted, hence grouped by owner

        self.recv_counts = np.bincount(kd.owner[remote], minlength=kd.size).astype(int)
        self.recv_displs = np.concatenate([[0], np.cumsum(self.recv_counts)[:-1]]).astype(int)

        # tell every owner which of its rows are wanted
        self.send_counts = np.array(kd.comm.alltoall(list(self.recv_counts)), dtype=int)
        self.send_displs = np.concatenate([[0], np.cumsum(self.send_counts)[:-1]]).astype(int)
        wanted = np.ascontiguousarray(remote, dtype=np.int64)
        asked = np.empty(int(np.sum(self.send_counts)), dtype=np.int64)
        kd.comm.Alltoallv([wanted, (self.recv_counts, self.recv_displs)],
                          [asked, (self.send_counts, self.send_displs)])
        self.send_rows = asked - kd.k0
        assert np.all((self.send_rows >= 0) & (self.send_rows < kd.nk_loc))

    def fetch(self, a_loc):
        """a_loc is this rank's block [nk_loc, ...]; returns the rows at idx."""
        kd = self.kd
        tail = a_loc.shape[1:]
        if kd.serial:
            return a_loc[self.idx].reshape(self.shape + tail)
        send = np.ascontiguousarray(a_loc[self.send_rows])
        recv = np.empty((len(self.uniq_remote),) + tail, dtype=a_loc.dtype)
        rowtype = kd._rowtype(a_loc)
        try:
            kd.comm.Alltoallv([send, (self.send_counts, self.send_displs), rowtype],
                              [recv, (self.recv_counts, self.recv_displs), rowtype])
        finally:
            rowtype.Free()
        rows = np.empty((self.nuniq,) + tail, dtype=a_loc.dtype)
        rows[self.uniq_local] = a_loc[self.local_rows]
        rows[self.uniq_remote] = recv
        return rows[self.inverse].reshape(self.shape + tail)
