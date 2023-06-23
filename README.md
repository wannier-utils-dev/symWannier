# symWannier

A collection of codes to construct Wannier functions using crystal symmetry.

## Installation

First, pw2wannier90 should be modified.

Copy src/pw2wannier90/qe7.0/pw2wannier90.f90 to your QE directory and compile it.
Alternatively, you can also use the develop version of QE in Gitlab (https://gitlab.com/QEF/q-e).

## How to use

SCF calculation and NSCF calculation to obtain wavefunctions in the irreducible Brillouin zone (IBZ).
```
pw.x < scf.in
pw.x < nscf.in
```
In the NSCF calculation, k points should be automatically generated as
```
K_POINTS {automatic}
8 8 8 0 0 0
```

Prepare prefix.win and generate prefix.nnkp.
```
wannier90.x -pp prefix
```
Run modified pw2wannier90.x with ```irr_bz = .true.``` to compute prefix.immn, prefix.iamn, prefix.ieig. Symmetry infomation is stored in prefix.isym.
```
pw2wannier90.x < pw2wan.in
```
Calculate Mmn, Amn and Eig in the full BZ using ```write_full_data.py```.
```
SYMWANNIERDIR=path/to/symWannier/src
export PYTHONPATH=$SYMWANNIERDIR
python $SYMWANNIERDIR/symwannier/write_full_data.py 
```
Run wannier90 as usual.
```
wannier90.x prefix
```

There are sample input and script files, in examples/Cu and examples/Fe.

## Paper

For more information, please see
[T. Koretsune, Comp. Phys. Comm. 285 108645 (2023).](https://doi.org/10.1016/j.cpc.2022.108645)

We hope that you cite this reference when you publish the results using this code.

