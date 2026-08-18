#!/bin/sh
#
# Carbon diamond, 6x6x6 k-mesh. Regenerates the IBZ inputs (diamond_mp6.i*)
# from scratch and then constructs symmetry-adapted Wannier functions.
#
# Edit the tool paths below and place C.pz-vbc.UPF in this directory.

MPIRUN=""
QE_PWX="path/to/pw.x"
QE_PW2WANNIER="path/to/pw2wannier90.x"
WANNIER="path/to/wannier90.x"
SYMWANNIER="../../src"

# 1. SCF and NSCF
$MPIRUN $QE_PWX < scf.in  > scf.out
$MPIRUN $QE_PWX < nscf.in > nscf.out

# 2. .nnkp (also produces the .win kpoints -> 8-decimal k such as 0.16666667)
$WANNIER -pp diamond_mp6

# 3. IBZ Mmn/Amn/Eig and symmetry info (irr_bz = .true.)
$MPIRUN $QE_PW2WANNIER < pw2wan.in > pw2wan.out

export PYTHONPATH=$SYMWANNIER

# 4. symmetry-adapted Wannier functions
#    --snap-kp (default true) snaps the .nnkp k-points to the exact i/mp_grid
#    rationals; use --snap-kp false to keep the finite-digit values.
python -m symwannier.wannierize -S -H diamond_mp6

# 5. quick check that the snapping matters on this non-dyadic grid
python verify.py
