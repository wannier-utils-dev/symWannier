#!/usr/bin/env python

import numpy as np
import logging

from symwannier.io_utils import open_text_or_gz

class Spn():
    """Reader and processor for spn files
    Spn(k) = <psi_mk| sigma_a |psi_nk>,  a = x, y, z

    The file holds, for every k-point, the upper triangle of the three Pauli
    matrices in the band basis, packed as pw2wannier90 writes it: the pair index
    runs over n <= m with m outermost, and the three Cartesian components run
    fastest. Only the formatted flavour is read, i.e. pw2wannier90 has to be run
    with spn_formatted = .true.

    With a Sym instance the irreducible file prefix.ispn is read instead and
    expanded to the full BZ. The expansion needs no representation matrix: the
    states of the full BZ are defined as the images of the irreducible ones with
    the band index untouched, so only the spin operator is rotated,

        S_a(g k) = sum_b A(g)_ab S_b(k)

    with A(g) the rotation the spinor matrix of g induces on the Pauli matrices.
    An operation with time reversal reverses the spin and conjugates,

        S_a(g k) = - sum_b A(g)_ab conj(S_b(k)).
    """

    def __init__(self, file_spn, nnkp, sym=None, log=None):
        """Load Spn data from file and expand it if symmetry data is given.

        Parameters
        ----------
        file_spn : str
            Path to spn file (optionally gzipped).
        nnkp : Nnkp
            Parsed nnkp object providing the k-points.
        sym : Sym, optional
            Symmetry data. If provided, the irreducible file is read and expanded.
        log : logging.Logger, optional
            Logger to use; if not provided a module logger is created.
        """
        self.log = log or logging.getLogger(__name__)
        if not self.log.handlers:
            logging.basicConfig(level=logging.INFO, format="%(message)s")

        self.nnkp = nnkp
        self.sym = sym

        fp, used_path = open_text_or_gz(file_spn, desc="spn file")
        self.log.debug(f"Reading spn from {used_path}")
        with fp:
            self._read_spn(fp)

    def _read_spn(self, fp):
        """Read spn contents from an open file-like object."""
        self.log.info("Reading spn file")
        lines = fp.readlines()
        num_bands, nk = [ int(x) for x in lines[1].split() ]
        npair = (num_bands * (num_bands + 1)) // 2
        dat = np.fromstring("".join(lines[2:]), sep=" ")
        if dat.size != 2 * 3 * npair * nk:
            raise ValueError(
                "spn file has {} numbers, expected {} for {} bands and {} k-points; "
                "pw2wannier90 has to be run with spn_formatted = .true."
                .format(dat.size, 2*3*npair*nk, num_bands, nk))
        dat = dat.reshape(nk, npair, 3, 2)
        packed = dat[:,:,:,0] + 1j * dat[:,:,:,1]

        # unpack the upper triangle into full hermitian matrices. The pairs are
        # packed column by column, (n,m) with n <= m and m outermost, which is the
        # order tril_indices gives when its two outputs are read as (m,n)
        im, iin = np.tril_indices(num_bands)
        spn = np.zeros([nk, 3, num_bands, num_bands], dtype=complex)
        spn[:, :, iin, im] = np.transpose(packed, axes=(0,2,1))
        il = np.tril_indices(num_bands, -1)
        spn[:, :, il[0], il[1]] = np.conj(spn[:, :, il[1], il[0]])

        self.num_bands = num_bands

        ####### simple case (without symmetry) #######
        if self.sym is None:
            self.nk = nk
            self.spn = spn

        ####### symmetrized case #######
        else:
            self.nk = self.sym.nkf
            self.spn = self.symmetrize_expand(spn)

    def spin_rotation(self, isym):
        """Rotation the symmetry operation induces on the Pauli matrices.

        Returns the real 3x3 matrix A with u^dagger sigma_a u = sum_b A_ab sigma_b,
        u being the spinor matrix of the operation. For an operation with time
        reversal the spinor matrix is the one search_symop composes, u_T conj(u).
        """
        u = self.sym.u_spin[isym]
        if self.sym.t_rev[isym] == 1:
            u = np.matmul(self.sym.uspin_T, np.conj(u))
        sigma = np.array([ [[0, 1], [1, 0]],
                           [[0, -1j], [1j, 0]],
                           [[1, 0], [0, -1]] ], dtype=complex)
        rot = np.einsum("alm,mn,bnp,pl->ab", sigma, np.conj(u).T, sigma, u, optimize=True) / 2
        return np.real(rot)

    def symmetrize_expand(self, spn_irk):
        """Generate Spn(full_k) from Spn(irr_k).

        The states of the full BZ are the images of the irreducible ones, with the
        band index untouched, so the band basis does not rotate and only the spin
        operator does.
        """
        spn = np.zeros([self.sym.nkf, 3, self.num_bands, self.num_bands], dtype=complex)
        for ik in range(self.sym.nkf):
            iks = self.sym.equiv[ik]
            isym = self.sym.equiv_sym[ik]
            rot = self.spin_rotation(isym)
            if self.sym.t_rev[isym] == 1:
                # time reversal reverses the spin and conjugates the matrix elements
                spn[ik] = - np.einsum("ab,bmn->amn", rot, np.conj(spn_irk[iks]), optimize=True)
            else:
                spn[ik] = np.einsum("ab,bmn->amn", rot, spn_irk[iks], optimize=True)
        return spn

    def write_spn(self, file_spn):
        """Write Spn data back to disk in wannier90 format (formatted)."""
        im, iin = np.tril_indices(self.num_bands)
        with open(file_spn, "w") as fp:
            fp.write("spn created by spn.py\n")
            fp.write("{} {}\n".format(self.num_bands, self.nk))
            for ik in range(self.nk):
                packed = self.spn[ik][:, iin, im]          # [3, npair]
                for ipair in range(packed.shape[1]):
                    for a in range(3):
                        v = packed[a, ipair]
                        fp.write("{:20.10E}{:20.10E}\n".format(v.real, v.imag))
