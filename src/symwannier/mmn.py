#!/usr/bin/env python

import numpy as np
import itertools
import logging

from symwannier.io_utils import open_text_or_gz

class Mmn:
    """Reader for mmn files (overlap matrices between neighboring k-points)."""

    def __init__(self, file_mmn, nnkp, sym=None, log=None):
        """Load Mmn data and set up k-space neighbor mappings.

        Parameters
        ----------
        file_mmn : str
            Path to mmn file (optionally gzipped).
        nnkp : Nnkp
            Parsed nnkp object.
        sym : Sym, optional
            Symmetry data for IBZ handling.
        log : logging.Logger, optional
            Logger instance.
        """
        self.log = log or logging.getLogger(__name__)
        if not self.log.handlers:
            logging.basicConfig(level=logging.INFO, format="%(message)s")
        self.nnkp = nnkp
        self.sym = sym

        # Determine if immn based on file extension
        ibz = file_mmn.endswith(".immn")

        fp, used_path = open_text_or_gz(file_mmn, desc="mmn file")
        self.log.debug(f"Reading mmn from {used_path}")
        with fp:
            self._read_mmn(fp, ibz)

        self._mmn_full_klist()

    def check_bands(self, eig=None, thr=1e-4):
        """Find the bands the symmetry expansion cannot reproduce.

        <u_mk|u_n,k+b> is by definition the adjoint of <u_n,k+b|u_mk>, so the
        expanded Mmn has to satisfy M(k,b) = M(k+b,-b)^dagger. A band whose
        degenerate multiplet is cut by num_bands breaks it: the symmetry mixes that
        band with one that is not in the file. The identity holds for the bands
        below the first such band and fails from there on, so the number of leading
        bands that satisfy it is the number of bands the expansion reproduces.

        Which bands those are changes from k-point to k-point, because a multiplet
        that is cut at one k-point need not be degenerate at another, so num_bands
        cannot be chosen to avoid the problem. An energy can: every band below the
        returned ceiling is reproduced at every k-point, which is what the outer
        disentanglement window (dis_win_max) has to stay below.

        Parameters
        ----------
        eig : ndarray, optional
            Eigenvalues of the full BZ, [nk, num_bands], used for the energy ceiling.
        thr : float, optional
            Tolerance on the identity.

        Returns
        -------
        n_bands : ndarray[int]
            Number of leading bands reproduced at each k-point.
        e_max : float or None
            Highest energy up to which every band is reproduced at every k-point;
            None when no band is affected or when eig is not given.
        """
        bvec = self.nnkp.bvec_crys
        opp = - np.ones([len(bvec)], dtype=int)
        for i in range(len(bvec)):
            for j in range(len(bvec)):
                if np.allclose(bvec[i], -bvec[j]):
                    opp[i] = j
                    break
        if np.any(opp < 0):
            self.log.info("the b-vector shell has no -b for every b; skipping the band check")
            return np.full([self.nk], self.num_bands), None

        n_bands = np.full([self.nk], self.num_bands, dtype=int)
        for ik in range(self.nk):
            dev = np.max(np.abs(self.mmn[ik,:,:,:]
                                - np.conj(self.mmn[self.kb2k[ik,:], opp, :, :]).transpose(0,2,1)), axis=0)
            run = 0.0
            for n in range(self.num_bands):
                run = max(run, dev[n,:n+1].max(), dev[:n+1,n].max())
                if run > thr:
                    n_bands[ik] = n
                    break

        affected = np.flatnonzero(n_bands < self.num_bands)
        if len(affected) == 0:
            self.log.info("Mmn: the expansion reproduces all %d bands", self.num_bands)
            return n_bands, None

        e_max = None
        if eig is not None:
            e_max = min( eig[ik, n_bands[ik]-1] for ik in affected if n_bands[ik] > 0 )
        self.log.warning(
            "  Warning: the symmetry expansion does not reproduce the highest bands at "
            "%d of %d k-points (as few as %d of %d bands); a degenerate multiplet is cut "
            "by num_bands there", len(affected), self.nk, n_bands.min(), self.num_bands)
        if e_max is not None:
            self.log.warning(
                "           every band below %.4f eV is reproduced at every k-point; keep "
                "the outer window dis_win_max below that value", e_max)
        return n_bands, e_max

    def write_mmn(self, file_mmn):
        """Write overlap matrices to file in wannier90 format."""
        with open(file_mmn, "w") as fp:
            fp.write("Mmn created by mmn.py\n")
            fp.write("{} {} {}\n".format(self.num_bands, self.nk, self.nb))
            for ik, ib in itertools.product(range(self.nk), range(self.nb)):
                k = self.sym.full_kpoints[ik]
                b = self.nnkp.bvec_crys[ib]
                ikb = self.kb2k[ik,ib]
                g = k + b - self.sym.full_kpoints[ikb]
                fp.write("{0}  {1}  {2[0]}  {2[1]}  {2[2]}\n".format(ik+1, ikb+1, np.round(g).astype("int")))
                # loop order: m (column) outer, n (row) inner  ->  mmn[ik,ib].T.ravel()
                mmn_kb = self.mmn[ik,ib,:,:].T.ravel()
                np.savetxt(fp, np.column_stack([mmn_kb.real, mmn_kb.imag]), fmt="%18.12f  %18.12f")

    def _read_mmn(self, fp, ibz):
        """Parse mmn file content and populate overlap matrices.

        Parameters
        ----------
        fp : file object
            Open file pointer to mmn file.
        ibz : bool
            True if file is an IBZ mmn (immn), False otherwise.

        Sets the following attributes:
        
        num_bands : int
            Number of bands.
        nks : int
            Number of irreducible k-points.
        nb : int
            Number of b-vectors.
        mmn : ndarray
            Overlap matrices M^k,b_mn.
        kb2k : ndarray
            Index of k+b.
        kpb_info : ndarray
            Information about k+b.
        """
        first_line = fp.readline()
        if not first_line:
            raise ValueError("Empty mmn file")

        if ibz:
            # immn requires symmetry information
            if self.sym is None:
                raise Exception("IBZ Mmn requires symmetry information.")
            self.log.info("Reading IBZ mmn file")
        else:
            # Regular mmn can be used with or without symmetry
            self.log.info("Reading mmn file")

        header = fp.readline()
        if not header:
            raise ValueError("mmn file missing header line")
        self.num_bands, self.nks, self.nb = [ int(x) for x in header.split() ]

        self.mmn = np.zeros([self.nks, self.nb, self.num_bands, self.num_bands], dtype=complex)
        self.kpb_info = np.zeros([self.nks, self.nb, 5], dtype=int)

        block_size = self.num_bands * self.num_bands
        for ik, ib in itertools.product(range(self.nks), range(self.nb)):
            head = fp.readline()
            if not head:
                raise ValueError("mmn file ended unexpectedly while reading header")
            d = [ int(x) for x in head.split() ]
            assert ik == d[0]-1, "{} {}".format(ik, d[0])
            self.kpb_info[ik,ib,:] = d

            block_lines = list(itertools.islice(fp, block_size))
            if len(block_lines) != block_size:
                raise ValueError("mmn file ended unexpectedly while reading data block")
            flat = np.fromstring(" ".join(block_lines), sep=" ")
            data = flat.reshape(self.num_bands, self.num_bands, 2)
            self.mmn[ik,ib,:,:] = data[:,:,0].T + 1j * data[:,:,1].T

    def _mmn_full_klist(self):
        if self.sym is None:
            assert self.nks == self.nnkp.nk
        else:
            assert self.sym.nkf == self.nnkp.nk
        self.nk = self.nnkp.nk   # nk for full k list
        assert self.nb == self.nnkp.nb

        self.kb2k = - np.ones([self.nk, self.nb], dtype=int)  # -1 becomes non-negative when defined
        ####### simple case (without symmetry) #######
        if self.sym is None:
            # Reorder one k-point at a time to avoid duplicating the full
            # MMN array, which can be large for many bands and k-points.
            for ik in range(self.nk):
                mmn_ik = self.mmn[ik, :, :, :].copy()
                kpb_info_ik = self.kpb_info[ik, :, :].copy()
                found = np.zeros(self.nb, dtype=bool)

                for ib in range(self.nb):
                    info = kpb_info_ik[ib, :]
                    bvec = self.nnkp.calc_bvec(info)
                    matches = np.flatnonzero(np.all(np.isclose(self.nnkp.bvec, bvec), axis=1))
                    if len(matches) != 1:
                        raise ValueError(
                            f"MMN block ({ik + 1}, {ib + 1}) matches "
                            f"{len(matches)} nnkp b-vectors"
                        )

                    ibt = matches[0]
                    if found[ibt]:
                        raise ValueError(
                            f"duplicate MMN b-vector {ibt + 1} "
                            f"at k-point {ik + 1}"
                        )

                    ikb = info[1] - 1
                    if not 0 <= ikb < self.nk:
                        raise ValueError(
                            f"MMN block ({ik + 1}, {ib + 1}) has invalid "
                            f"neighbor k-point {info[1]}"
                        )

                    found[ibt] = True
                    self.mmn[ik, ibt, :, :] = mmn_ik[ib, :, :]
                    self.kpb_info[ik, ibt, :] = info
                    self.kb2k[ik, ibt] = ikb

                if not np.all(found):
                    missing = np.flatnonzero(~found)[0]
                    raise ValueError(
                        f"missing MMN b-vector {missing + 1} "
                        f"at k-point {ik + 1}"
                    )

        ####### symmetrized case #######
        else:
            mmn = np.zeros([self.nk, self.nb, self.num_bands, self.num_bands], dtype=complex)

            bvec_equiv = self.sym.kpoint_equiv_info(self.nnkp.bvec_crys)
            for ikf, kf in enumerate(self.sym.full_kpoints):
                # s[isym] . irr_k[iki] = full_k[ikf]
                iki = self.sym.equiv[ikf]
                isym1 = self.sym.equiv_sym[ikf]
                ki = self.sym.irr_kpoints[iki]
                for ibf, bf in enumerate(self.nnkp.bvec_crys):
                    ## s[isym] . bvec[ibi] ~= bvec[ bvec_equiv[ibf,isym1] ]  i.e. S.bi = bf
                    ibi = bvec_equiv[ibf,isym1]
                    bi = self.nnkp.bvec_crys[ibi]
                    ikbi = self.sym.search_ik_full(ki + bi)
                    ikbf = self.sym.search_ik_full(kf + bf)
                    isym2 = self.sym.equiv_sym[ikbi]
                    isym3 = self.sym.equiv_sym[ikbf]
                    # g[isym2]^-1 g[isym1]^-1 g[isym3] = s[isym,:,:], ft[isym,:,:] + tdiff
                    isym, factor, tdiff = self.sym.search_symop( [[isym2, -1], [isym1, -1], [isym3, 1]] )

                    ikbi_eq = self.sym.equiv[ ikbi ]

                    if np.sum(np.abs(self.sym.repmat[ikbi_eq,isym,:,:])) < 1e-5:
                        print(ikbi)
                        print(ikbi_eq)

                    self.kb2k[ikf,ibf] = ikbf
                    if self.sym.t_rev[isym2] == 1:
                        mmn[ikf,ibf,:,:] = np.einsum("mn,np->mp", self.mmn[iki,ibi,:,:], np.conj(self.sym.repmat[ikbi_eq,isym,:,:]), optimize=True) * factor
                    else:
                        mmn[ikf,ibf,:,:] = np.einsum("mn,np->mp", self.mmn[iki,ibi,:,:], self.sym.repmat[ikbi_eq,isym,:,:], optimize=True) * factor
                    if self.sym.t_rev[isym1] == 1:
                        mmn[ikf,ibf,:,:] = np.conj(mmn[ikf,ibf,:,:])

                    # e^{-i b_i T} 
                    #   from g0^-1 e^{-i b_f r} => e^{i b_i r} g0^-1 x e^{-i b_i T}
                    # e^{i (k_i+b_i) Tdiff}
                    #   from h = g0^-1(ki+bi) g0^-1(kf) g0(kf+bf)
                    kbi_eq = self.sym.irr_kpoints[ikbi_eq]
                    arg1 = - np.dot(bi, self.sym.ft[isym1,:])
                    arg2 = - np.dot(kbi_eq, tdiff)
                    if self.sym.t_rev[isym2] == 1:
                        arg2 *= -1
                    if self.sym.t_rev[isym1] == 1:
                        arg1 *= -1
                        arg2 *= -1
                    if self.sym.t_rev[isym] == 1:
                        arg2 *= -1
                    phase = np.exp( 2j*np.pi* (arg1 + arg2) )
                    mmn[ikf,ibf,:,:] *= phase

            self.mmn = mmn

        assert np.all( self.kb2k[:,:] >= 0 )  # check all kb2k is defined

