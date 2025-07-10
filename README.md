# symWannier

A collection of codes to construct Wannier functions using crystal symmetry.
By using this package, users only need to calculate the wave functions, overlap matrices (Mmn), and projection matrices (Amn) in the irreducible Brillouin zone (IBZ).
Mmn and Amn files for the full Brillouin zone, required by the [wannier90](http://www.wannier.org/) package, can be generated using the symmetry information.
Users can also use the Python scripts included in this package to construct maximally-localized Wannier functions.
There is an option to calculate symmetry-adapted Wannier functions.

## Installation

You can install this package using pip:

```
pip install symWannier
```

## IBZ calculation using wannier90

Run SCF and NSCF calculations to obtain wave functions in the irreducible Brillouin zone (IBZ).
```
pw.x < scf.in
pw.x < nscf.in
```

In the NSCF calculation, k-points should be generated automatically as
```
K_POINTS {automatic}
8 8 8 0 0 0
```

Prepare prefix.win and generate prefix.nnkp. You can use the same prefix.win that is used in the original wannier90 calculation.
```
wannier90.x -pp prefix
```

Run pw2wannier90.x with ```irr_bz = .true.``` to compute prefix.immn, prefix.iamn, prefix.ieig. Symmetry infomation is stored in prefix.isym.

Note: Quantum ESPRESSO version 7.3 or later is required for this step.
```
pw2wannier90.x < pw2wan.in
```

Calculate Mmn, Amn and Eig in the full BZ using ```expand_wannier_inputs.py```.
```
symwannier expand prefix
```

Run wannier90 as usual.
```
wannier90.x prefix
```

There are sample input and script files, in examples/Cu and examples/Fe.


## Symmetry adapted Wannier functions

Once prefix.immn, prefix.iamn, prefix.ieig, and prefix.isym are obtained, you can construct Wannier functions as follows.
```
symwannier wannierize -s prefix
```
To construct symmetry-adapted Wannier functions, use the "-S" option instead of "-s".
```
symwannier wannierize -S prefix
```


## Paper

For more information, please see
[T. Koretsune, Comp. Phys. Comm. 285 108645 (2023).](https://doi.org/10.1016/j.cpc.2022.108645)

We hope that you cite this reference when you publish the results using this code.

