#!/usr/bin/env python

import numpy as np
import itertools
import scipy.linalg
import logging
import os


class Nnkp:
    """Reader for nnkp files (neighbor k-point information from wannier90)."""
    def __init__(self, file_nnkp, log=None, snap_kp=True):
        """Parse nnkp file and compute b-vector list.

        Parameters
        ----------
        file_nnkp : str
            Path to nnkp file.
        log : logging.Logger, optional
            Logger instance.
        snap_kp : bool, optional
            Snap the k-points to the exact i/mp_grid rationals of the
            Monkhorst-Pack grid (default: True). See ``_snap_kpoints``.
        """
        self.log = log or logging.getLogger(__name__)
        if not self.log.handlers:
            logging.basicConfig(level=logging.INFO, format="%(message)s")

        self.read_file(file_nnkp)

        if snap_kp:
            self._snap_kpoints()

        self.calc_bvec_list()

    def _snap_kpoints(self):
        """Snap the k-points to the exact rationals of a Gamma-centered grid.

        wannier90 writes the .nnkp k-points at a fixed 8-decimal format, so a
        value such as 1/12 is stored as 0.08333333 (~3e-9 below the exact
        double). On a (possibly half-shifted) Monkhorst-Pack grid the exact
        k-points are (2*i + s)/(2*N); we recover them by rounding k*2N to the
        nearest integer. This restores the bit-exact double used internally by
        the DFT code and avoids a ~1e-5 error in the H(k)<->H(R) Fourier
        transforms (write_hr/write_tb). Non-standard grids are left unchanged.
        """
        k = self.kpoints
        # infer the grid size per direction from the number of distinct values
        n = np.array([len(np.unique(np.round(k[:, d], 6))) for d in range(3)])
        two_n = 2 * n
        m = np.round(k * two_n)                      # = 2*i + shift (integers)
        resid = np.max(np.abs(k * two_n - m))        # distance from the integers
        if resid < 1e-10:
            return                                   # already exact (e.g. dyadic grid)
        if resid > 1e-4:
            self.log.info("  Note: k-points are not on a standard Monkhorst-Pack "
                          "grid; leaving them unchanged (--snap-kp).")
            return
        self.log.info("  Note: snapping .nnkp k-points to exact i/mp_grid rationals "
                      "(finite-digit output, residual %.1e). Disable with "
                      "--snap-kp false.", resid)
        self.kpoints = m / two_n


    def read_file(self, file_nnkp):
        """Parse nnkp file and extract lattice, kpoints, neighbors, and projections.

        Parameters
        ----------
        file_nnkp : str
            Path to wannier90 nnkp file.

        Raises
        ------
        FileNotFoundError
            If nnkp file cannot be found.
        """
        if not os.path.exists(file_nnkp):
            raise FileNotFoundError(f"nnkp file not found: {file_nnkp}")
        try:
            with open(file_nnkp) as fp:
                lines = fp.readlines()
                for i, line in enumerate(lines):
                    if "begin real_lattice" in line:
                        self.a = np.genfromtxt(lines[i+1:i+4], dtype=float)

                    if "begin recip_lattice" in line:
                        self.b = np.genfromtxt(lines[i+1:i+4], dtype=float)

                    if "begin kpoints" in line:
                        self.nk = int(lines[i+1])
                        self.kpoints = np.genfromtxt(lines[i+2:i+2+self.nk], dtype=float)

                    if "begin nnkpts" in line:
                        self.nb = int(lines[i+1])
                        dat = np.genfromtxt(lines[i+2:i+2+self.nk*self.nb], dtype=int)
                        self.nnkpts = dat.reshape(self.nk, self.nb, 5)

                    if "begin projections" in line or "begin spinor_projections" in line:
                        spinors = ("begin spinor_projections" in line)
                        self.num_wann = int(lines[i + 1])
                        self.proj = []
                        atom_orb_strlist = []
                        atom_pos_strlist = []
                        atom_pos_r_strlist = []
                        self.nw2n = np.zeros([self.num_wann], dtype=int)
                        self.nw2l = np.zeros([self.num_wann], dtype=int)
                        self.nw2m = np.zeros([self.num_wann], dtype=int)
                        self.nw2r = np.zeros([self.num_wann, 3], dtype=float)
                        # read projections
                        for j in range(self.num_wann):
                            if (spinors):
                                proj_str = lines[i + 2 + 3 * j]
                            else:
                                proj_str = lines[i + 2 + 2 * j]
                            proj_dat = proj_str.split()
                            self.nw2l[j] = int(proj_dat[3])
                            self.nw2m[j] = int(proj_dat[4])
                            self.nw2r[j, :] = [float(x) for x in proj_dat[0:3]]
                            atom_orb_strlist.append(proj_str[0:40])
                            atom_pos_strlist.append(proj_str[0:35])
                        # set atom_pos_r, atom_pos, atom_orb
                        #   for example, Fe case
                        #   atom_pos_r: [[0.0, 0.0, 0.0]]
                        #   atom_pos: [[0, 1, 2]]
                        #   atom_orb: [[0, 1], [2, 3, 4, 5, 6, 7], [8, 9, 10, 11, 12, 13, 14, 15, 16, 17]]
                        atom_orb_uniq = sorted(set(atom_orb_strlist), key=atom_orb_strlist.index)
                        atom_pos_uniq = sorted(set(atom_pos_strlist), key=atom_pos_strlist.index)
                        self.atom_orb = []
                        for orb_str in atom_orb_uniq:
                            indexes = [j for j, x in enumerate(atom_orb_strlist) if x == orb_str]
                            self.atom_orb.append(indexes)
                        self.atom_pos = []
                        self.atom_pos_r = []
                        for pos_str in atom_pos_uniq:
                            indexes = [j for j, x in enumerate(atom_orb_uniq) if pos_str in x]
                            self.atom_pos.append(indexes)
                            self.atom_pos_r.append([float(x) for x in pos_str.split()[0:3]])
                        # print ("atom_pos_r: " + str(self.atom_pos_r))
                        # print ("atom_pos: " + str( self.atom_pos))
                        # print ("atom_orb: " + str(self.atom_orb))
                        self.natom = len(self.atom_pos_r)
                        for i, pos in enumerate(self.atom_pos):
                            for p in pos:
                                for n in self.atom_orb[p]:
                                    self.nw2n[n] = i
                        # for j in range(self.num_wann):
                        #    print("nw {:3d} : n = {:3d}, l = {:3d}, m = {:3d}".format(j, self.nw2n[j], self.nw2l[j], self.nw2m[j]))

        except Exception as e:
            print("failed to read: " + file_nnkp)
            print('type:' + str(type(e)))
            print('args:' + str(e.args))
            print(str(e))

    def calc_bvec(self, d):
        """
        d : array of integer [5]
        return : bvec in cartesian coordinates
        """
        k = self.kpoints[ d[0]-1 ]
        kb = self.kpoints[ d[1]-1 ]
        bvec = kb - k + np.array( d[2:5] )
        return np.matmul(bvec, self.b)

    def bvec_num(self, b):
        for ib in range(self.nb):
            if np.allclose(self.bvec[ib,:], b):
                return ib
        assert False, b

    def calc_bvec_list(self):
        self.bvec = np.zeros([self.nb, 3])
        self.bvec_crys = np.zeros([self.nb, 3])
        self.wb = np.zeros([self.nb])
        for i in range(self.nb):
            self.bvec[i,:] = self.calc_bvec( self.nnkpts[0,i,:] )
            self.bvec_crys[i,:] = self.k_cart2crys( self.bvec[i,:] )
            #k = self.nnkp.kpoints[self.nnkp.nnkpts[0,i,0]-1]
            #kb = self.nnkp.kpoints[self.nnkp.nnkpts[0,i,1]-1]
            #self.bvec[i,:] = kb - k + np.array(self.nnkp.nnkpts[0,i,2:5])

        #self.bvec = np.matmul(self.bvec, self.nnkp.b)  # crystal coordinate -> cartesian
        #for i in range(self.nb):
        #    self.wb[i] = 3/self.nb/np.linalg.norm(self.bvec[i,:])**2
        bmat = np.zeros([self.nb, 9])
        for i in range(self.nb):
            bmat[i,:] = [ self.bvec[i,a]*self.bvec[i,b] for a,b in itertools.product(range(3), range(3)) ]
        delta_ab =  np.array([ a==b for a,b in itertools.product(range(3), range(3)) ]).astype(int)
        self.wb = np.matmul(delta_ab, scipy.linalg.pinv(bmat))

    def k_crys2cart(self, k):
        return np.matmul(k, self.b)

    def k_cart2crys(self, k):
        return np.matmul(self.a, k)/(2*np.pi)

