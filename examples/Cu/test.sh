#!/bin/sh

MPIRUN="mpirun -n 16"
QE_PWX="path/to/pw.x"
QE_PW2WANNIER="path/to/pw2wannier90.x"
WANNIER="path/to/wannier90.x"
SYMWANNIER="../../symWannier/src"

$MPIRUN $QE_PWX < scf.in > scf.out
$MPIRUN $QE_PWX < nscf.in > nscf.out

$WANNIER -pp pwscf

$MPIRUN $QE_PW2WANNIER < pw2wan.in > pw2wan.out

export PYTHONPATH=$SYMWANNIER

# restore mmn, amn
python $SYMWANNIER/symwannier/write_full_data.py pwscf
cp pwscf.win w90/

cd w90
$MPIRUN $WANNIER pwscf

cd ..

# symmetry adapted Wannier
python $SYMWANNIER/symwannier/wannierize.py -S pwscf
