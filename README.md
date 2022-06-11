# symWannier

A collection of codes to construct Wannier functions using crystal symmetry.

## Installation

We need to modify pw2wannier90.

Copy src/pw2wannier90/qe7.0/pw2wannier90.f90 to your QE directory and compile it.

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
Prepare prefix.win and generate prefix.nnkp. Uniform k points in prefix.win should be consistent with NSCF calc (-0.5 <= kx, ky, kz < 0.5). 
```
wannier90.x -pp prefix
```
Run modified pw2wannier90.x with ```irr_bz = .true.``` to compute prefix_ibz.mmn, prefix_ibz.amn, prefix_ibz.eig. Symmetry infomation is stored in prefix_sym.dat.
```
pw2wannier90.x < pw2wan.in
```
Calculate Mmn, Amn and Eig in the full BZ.
```
python write_full_data.py 
```
Run wannier90 as usual.
```
wannier90.x prefix
```

There are sample input files and a script file, test.sh, in examples/Cu and examples/Fe.


