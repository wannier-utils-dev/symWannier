# symWannier

A collection of codes to construct Wannier functions using crystal symmetry.
By using this package, user only have to calculate the wave functions, overlap matrices (Mmn) and projection matrices (Amn) in the irreducible Brillouin zone (IBZ).
Mmn and Amn files in the full BZ for wannier90 package (http://www.wannier.org/) can be generated using the symmetry information.
Users can also use the Python scripts included in this package to construct maximally-localized Wannier functions.
There is an option to calculate symmetry-adapted Wannier functions.

## Installation

First, please download the scripts:
```
git clone git@github.com:wannier-utils-dev/symWannier.git
```

In addition, pw2wannier90 should be modified.

Copy src/pw2wannier90/qe7.0/pw2wannier90.f90 to your QE directory and compile it.
Alternatively, you can also use the develop version of QE in Gitlab (https://gitlab.com/QEF/q-e).

## IBZ calculation using wannier90

Run SCF calculation and NSCF calculation to obtain wavefunctions in the irreducible Brillouin zone (IBZ).
```
pw.x < scf.in
pw.x < nscf.in
```
In the NSCF calculation, k points should be automatically generated as
```
K_POINTS {automatic}
8 8 8 0 0 0
```

Prepare prefix.win and generate prefix.nnkp. You can use the same prefix.win that is used in the original wannier90 calculation.
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

python $SYMWANNIERDIR/symwannier/write_full_data.py prefix
```
Run wannier90 as usual.
```
wannier90.x prefix
```

There are sample input and script files, in examples/Cu and examples/Fe.


## Symmetry adapted Wannier functions

Once prefix.immn, prefix.iamn, prefix.ieig, and prefix.isym are obtained, you can construct Wannier functions as follows.
```
SYMWANNIERDIR=path/to/symWannier/src
export PYTHONPATH=$SYMWANNIERDIR

python $SYMWANNIERDIR/symwannier/wannierize.py -s prefix
```
To construct symmetry-adapted Wannier functions, use "-S" option instead of "-s".
```
python $SYMWANNIERDIR/symwannier/wannierize.py -S prefix
```


## Paper

For more information, please see
[T. Koretsune, Comp. Phys. Comm. 285 108645 (2023).](https://doi.org/10.1016/j.cpc.2022.108645)

We hope that you cite this reference when you publish the results using this code.

