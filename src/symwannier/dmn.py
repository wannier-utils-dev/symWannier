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
    operations of prefix.isym are listed in it. Time reversal is not excluded by
    that: what the file has to express is the connection between an irreducible
    k-point and each member of its orbit, and for a spatial R that connection is
    g2 g1^-1 with g1, g2 the operations symmetrize_expand used to define the two.
    Its kind is t_rev(g1) + t_rev(g2), so it is linear when both are spatial and
    also when both are time-reversed - the two conjugations then cancel and the
    matrices are simply conjugated, which is how 16 of the 64 k-points of GaAs are
    built. Only a mismatch cannot be expressed, and _build raises there.

    Note on the phase convention. The file this writes is not identical to the
    one pw2wannier90 writes for the same system, although it is equivalent. Both
    have the form diag(phase) rotmat^dagger with the same rotation matrices
    (wws = rotmat^T to machine precision) and the same dependence of the phase on
    the Wannier index, which is the lattice shift Rshift of the projection
    centers; they differ by one scalar of modulus one per (operation, k-point).
    That scalar cancels in every use wannier90 makes of the file, because the
    relation above and everything built on it are invariant under
    (d, D) -> (lambda d, lambda D) - d and D appear in a lambda/lambda* pair in
    every one of them - and the final spreads agree to 6e-09. Of the
    48 operations of diamond, 36 of the scalars are a constant lattice offset,
    42 are that plus the phase an umklapp picks up from the fractional
    translation, and the remaining 6 - all with ft = (0,0,-1/2), the fractional
    translation along the offset of one of the projection centers - are not
    explained by either. pw2wannier90 routes the lattice vector through the
    inverse operation (vps2t(:, ips2p(ip, invs(isym)), isym)), which is the most
    likely home of the residual sign; that bookkeeping is not reproduced here.
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

    def _orbits(self):
        """Split the full mesh into orbits of the spatial operations.

        The irreducible k-points of prefix.isym are used as seeds, so when the
        spatial operations already cover the mesh the orbits are exactly their
        stars and the irreducible set of the dmn is the one of prefix.isym. When
        time reversal is needed as well - prefix.isym reduces the mesh further
        than the spatial operations can - the stars break into several orbits and
        the dmn gets more irreducible k-points than prefix.isym has.
        """
        sym = self.sym
        ik2ir = - np.ones([self.nk], dtype=int)
        ir2ik = []
        for ik in list(sym.iks2ik) + list(range(self.nk)):
            if ik2ir[ik] >= 0:
                continue
            ir = len(ir2ik)
            ir2ik.append(ik)
            for isym in self.isym_list:
                ik2 = sym.search_ik_full(np.dot(sym.s[isym], sym.full_kpoints[ik]))
                if ik2ir[ik2] < 0:
                    ik2ir[ik2] = ir
        return ik2ir, np.array(ir2ik, dtype=int)

    def _build(self):
        """Build kptsym and the d matrices for every irreducible k and operation."""
        sym = self.sym
        Rmat, Rshift, _ = self.amn.projection_sym_mat()

        self.isym_list = self._spatial_symops()
        self.nsym = len(self.isym_list)
        self.ik2ir, self.ir2ik = self._orbits()
        self.nkirr = len(self.ir2ik)
        self.log.info("dmn: {} spatial symmetry operations out of {}, {} irreducible k-points"
                      .format(self.nsym, sym.nsym, self.nkirr))

        self.kptsym = np.zeros([self.nsym, self.nkirr], dtype=int)
        self.d_matrix_wann = np.zeros([self.nsym, self.nkirr, self.num_wann, self.num_wann], dtype=complex)
        self.d_matrix_band = np.zeros([self.nsym, self.nkirr, self.num_bands, self.num_bands], dtype=complex)

        for ir in range(self.nkirr):
            # the k-point of the dmn, the irreducible k-point of prefix.isym it comes
            # from, and the operation symmetrize_expand used to get there
            ik1 = self.ir2ik[ir]
            iks = sym.equiv[ik1]
            ks = sym.irr_kpoints[iks]
            rep1 = sym.equiv_sym[ik1]
            m1 = self._rot(Rmat, Rshift, rep1, ks)

            for i, isym in enumerate(self.isym_list):
                ik2 = sym.search_ik_full(np.dot(sym.s[isym], sym.full_kpoints[ik1]))
                self.kptsym[i, ir] = ik2
                rep2 = sym.equiv_sym[ik2]

                # R rep1 = rep2 h, with h an operation of the little group of the
                # irreducible k-point. The sign search_symop returns for the spinor
                # double group is not used: both matrices below are expressed through
                # the same h, so it appears on either side of the relation and cancels.
                isym_h, _, _ = sym.search_symop([[rep2, -1], [isym, 1], [rep1, 1]])
                if sym.t_rev[isym_h] != 0:
                    raise ValueError(
                        "k-points {} and {} of the full mesh are related by a spatial "
                        "operation but reached from the irreducible mesh with and without "
                        "time reversal; the dmn file cannot represent an antiunitary "
                        "operation".format(ik1+1, ik2+1))

                # band side: d(R,k) = repmat[k, h]
                d = sym.repmat[iks, isym_h, :, :]

                # Wannier side: the rotations and phases that symmetrize_expand applies,
                # built from the same decomposition. Taking the rotation of R directly
                # would be wrong whenever R k leaves the first Brillouin zone: the
                # rotation matrices pick up a phase from the fractional translation when
                # the umklapp vector is not zero, and that phase only comes out right if
                # the factors are kept apart.
                mat = np.conj(m1).T @ self._rot(Rmat, Rshift, isym_h, ks) \
                                   @ self._rot(Rmat, Rshift, rep2, ks)

                # where symmetrize_expand conjugated Amn, it did so at both k-points,
                # so the relation stays linear and only the matrices are conjugated
                if sym.t_rev[rep1] == 1:
                    d = np.conj(d)
                    mat = np.conj(mat)

                self.d_matrix_band[i, ir, :, :] = d
                self.d_matrix_wann[i, ir, :, :] = np.conj(mat).T

    @staticmethod
    def _rot(Rmat, Rshift, isym, k):
        """Rotation of the projections for one operation, with the phase at k.

        The same matrix that Amn.symmetrize_expand contracts Amn with.
        """
        phase2 = np.einsum("a,na->n", k, Rshift[isym, :, :], optimize=True)
        phase = np.exp(-1j * 2*np.pi * phase2)
        return np.einsum("ln,n->ln", Rmat[isym,:,:], phase, optimize=True)

    def check(self, amn_full, thr=1e-4):
        """Check the dmn against Amn in the full BZ.

        Verifies the relation wannier90 relies on,

            A(Rk) = d(R,k) A(k) D^dagger(R,k),

        for every irreducible k and every symmetry operation. Returns the largest
        deviation found.

        The relation can only hold as well as the Amn written by pw2wannier90 is
        itself symmetric under the little group, which sets the deviation of a sound
        calculation: 5e-14 (H), 1e-08 (diamond), 2e-07 (Sn), 1e-06 (GaAs, whose Amn
        is only symmetric to 6e-07 itself). The default threshold sits above those
        and well below the 1e-02 of a case where a degenerate multiplet is cut by
        num_bands, which makes repmat non-unitary and genuinely breaks the relation.
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
