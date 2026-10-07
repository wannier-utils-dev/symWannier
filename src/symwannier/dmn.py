#!/usr/bin/env python

import numpy as np
import logging

class Dmn():
    """Generator for the wannier90 dmn file (site symmetry).

    wannier90 uses the dmn file as

        U(Rk) = d(R,k) U(k) D^dagger(R,k)

    where k runs over the irreducible k-points, R over the symmetry operations,
    Rk = kptsym(R,k), d = d_matrix_band (num_bands) and D = d_matrix_wann
    (num_wann).  Both matrices follow from the IBZ data: d from the
    representation matrices repmat of prefix.isym, D from the rotation matrices
    rotmat and the lattice shifts of the projection centers, i.e. from the same
    objects that Amn.symmetrize_expand uses to generate Amn in the full BZ.  The
    file is therefore consistent with the Amn/Mmn/Eig written by
    expand_wannier_inputs.

    The dmn file has no place for an antiunitary operation, so only the spatial
    operations of prefix.isym are used, and every k-point of the full mesh has
    to be reachable from the irreducible mesh without time reversal.
    """

    def __init__(self, nnkp, sym, amn, log=None):
        """Build the dmn data from the symmetry information.

        Parameters
        ----------
        nnkp : Nnkp
            Parsed nnkp object.
        sym : Sym
            Symmetry data read from prefix.isym.
        amn : Amn
            Amn object, used for the projection rotation matrices and for
            num_bands / num_wann.
        log : logging.Logger, optional
            Logger to use; if not provided a module logger is created.
        """
        self.log = log or logging.getLogger(__name__)
        if not self.log.handlers:
            logging.basicConfig(level=logging.INFO, format="%(message)s")

        self.nnkp = nnkp
        self.sym = sym
        self.amn = amn
        self.num_bands = amn.num_bands
        self.num_wann = amn.num_wann
        self.nk = sym.nkf
        self.nkirr = sym.nks

        self._build()

    def _spatial_symops(self):
        """Symmetry operations without time reversal, the identity first."""
        isym_list = [ isym for isym in range(self.sym.nsym) if self.sym.t_rev[isym] == 0 ]
        identity = [ isym for isym in isym_list
                     if np.array_equal(self.sym.s[isym], np.eye(3, dtype=int))
                     and np.allclose(self.sym.ft[isym], 0) ]
        if len(identity) != 1:
            raise ValueError("prefix.isym does not contain a unique identity operation")
        # wannier90 assumes that the first operation is the identity
        isym_list.remove(identity[0])
        return [identity[0]] + isym_list

    def _representatives(self):
        """For each k of the full mesh, a spatial operation mapping its irreducible k to it."""
        sym = self.sym
        rep = - np.ones([self.nk], dtype=int)
        for ik, k in enumerate(sym.full_kpoints):
            ks = sym.irr_kpoints[ sym.equiv[ik] ]
            for isym in self.isym_list:
                kdiff = np.dot(sym.s[isym], ks) - k
                if np.allclose(kdiff, np.round(kdiff)):
                    rep[ik] = isym
                    break
        missing = np.flatnonzero(rep < 0)
        if len(missing) > 0:
            raise ValueError(
                "k-point {} of the full mesh is only reachable from the irreducible mesh "
                "with time reversal; the dmn file cannot represent an antiunitary "
                "operation".format(missing[0]+1)
            )
        return rep

    def _build(self):
        """Build kptsym and the d matrices for every irreducible k and operation."""
        sym = self.sym
        Rmat, Rshift, _ = self.amn.projection_sym_mat()

        self.isym_list = self._spatial_symops()
        self.nsym = len(self.isym_list)
        self.log.info("dmn: {} spatial symmetry operations out of {}".format(self.nsym, sym.nsym))

        rep = self._representatives()

        # wannier90 numbers the k-points of the full mesh; its irreducible set is
        # the one of prefix.isym
        self.ik2ir = sym.equiv
        self.ir2ik = sym.iks2ik

        self.kptsym = np.zeros([self.nsym, self.nkirr], dtype=int)
        self.d_matrix_wann = np.zeros([self.nsym, self.nkirr, self.num_wann, self.num_wann], dtype=complex)
        self.d_matrix_band = np.zeros([self.nsym, self.nkirr, self.num_bands, self.num_bands], dtype=complex)

        for iks in range(self.nkirr):
            k = sym.irr_kpoints[iks]
            ik1 = self.ir2ik[iks]
            for i, isym in enumerate(self.isym_list):
                ik2 = sym.search_ik_full(np.dot(sym.s[isym], k))
                self.kptsym[i, iks] = ik2

                # R = rep(Rk) h, with h an operation of the little group of k.
                # The sign search_symop returns for the spinor double group is not
                # used: both matrices below are expressed through the same h, so it
                # appears on either side of the relation and cancels.
                isym_h, _, _ = sym.search_symop(
                    [[rep[ik2], -1], [isym, 1], [rep[ik1], 1]] )

                # band side: d(R,k) = repmat[k, h]
                self.d_matrix_band[i, iks, :, :] = sym.repmat[iks, isym_h, :, :]

                # Wannier side: the rotation and the phase that symmetrize_expand
                # applies, built from the same decomposition. Taking the rotation of R
                # directly would be wrong whenever R k leaves the first Brillouin zone:
                # the rotation matrices pick up a phase from the fractional translation
                # when the umklapp vector is not zero, and that phase only comes out
                # right if the two factors are kept apart.
                mat = np.einsum("lm,mn->ln", self._rot(Rmat, Rshift, isym_h, k),
                                self._rot(Rmat, Rshift, rep[ik2], k), optimize=True)
                self.d_matrix_wann[i, iks, :, :] = np.conj(mat).T

    @staticmethod
    def _rot(Rmat, Rshift, isym, k):
        """Rotation of the projections for one operation, with the phase at k.

        The same matrix that Amn.symmetrize_expand contracts Amn with.
        """
        phase2 = np.einsum("a,na->n", k, Rshift[isym, :, :], optimize=True)
        phase = np.exp(-1j * 2*np.pi * phase2)
        return np.einsum("ln,n->ln", Rmat[isym,:,:], phase, optimize=True)

    def check(self, amn_full, thr=1e-6):
        """Check the dmn against Amn in the full BZ.

        Verifies the relation wannier90 relies on,

            A(Rk) = d(R,k) A(k) D^dagger(R,k),

        for every irreducible k and every symmetry operation. Returns the largest
        deviation found.
        """
        diff = 0.0
        for iks in range(self.nkirr):
            a1 = amn_full[ self.ir2ik[iks], :, :]
            for i in range(self.nsym):
                a2 = amn_full[ self.kptsym[i, iks], :, :]
                rhs = np.einsum("ml,ln,pn->mp", self.d_matrix_band[i,iks], a1,
                                np.conj(self.d_matrix_wann[i,iks]), optimize=True)
                diff = max(diff, np.max(np.abs(a2 - rhs)))
        if diff > thr:
            self.log.warning("  Warning: dmn is not consistent with Amn (max deviation %.5e)", diff)
        else:
            self.log.info("dmn vs Amn: max deviation = %.5e", diff)
        return diff

    def write_dmn(self, file_dmn):
        """Write the dmn file in wannier90 format."""
        with open(file_dmn, "w") as fp:
            fp.write("dmn created by dmn.py\n")
            fp.write("{:9d}{:9d}{:9d}{:9d}\n".format(
                self.num_bands, self.nsym, self.nkirr, self.nk))

            fp.write("\n")
            self._write_int(fp, self.ik2ir + 1)
            fp.write("\n")
            self._write_int(fp, self.ir2ik + 1)
            for iks in range(self.nkirr):
                fp.write("\n")
                self._write_int(fp, self.kptsym[:, iks] + 1)

            # wannier90 reads d_matrix_wann(num_wann, num_wann, nsym, nkptirr) and
            # d_matrix_band(num_bands, num_bands, nsym, nkptirr) with a list-directed
            # read, i.e. in Fortran order with the first index running fastest
            for mat in (self.d_matrix_wann, self.d_matrix_band):
                for iks in range(self.nkirr):
                    for i in range(self.nsym):
                        fp.write("\n")
                        for n in range(mat.shape[3]):
                            for m in range(mat.shape[2]):
                                v = mat[i, iks, m, n]
                                fp.write(" ({:18.10E},{:18.10E})\n".format(v.real, v.imag))

    @staticmethod
    def _write_int(fp, values):
        """Write an integer list, 10 per line, as wannier90's own files do."""
        for i in range(0, len(values), 10):
            fp.write("".join("{:9d}".format(x) for x in values[i:i+10]) + "\n")
