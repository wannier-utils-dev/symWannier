# symWannier

A collection of codes to construct Wannier functions using crystal symmetry.
By using this package, users only need to calculate the wave functions, overlap matrices (Mmn), and projection matrices (Amn) in the irreducible Brillouin zone (IBZ).
Mmn and Amn files for the full Brillouin zone, required by the [wannier90](http://www.wannier.org/) package, can be generated using the symmetry information.
Users can also use the Python scripts included in this package to construct maximally-localized Wannier functions.
There is an option to calculate symmetry-adapted Wannier functions.

## Installation

If you are working from a local checkout of this repository, install it with:

```bash
pip install -e .
```

For development, including the test dependencies, use:

```bash
pip install -e '.[dev]'
```

This installs the Python package and the `symwannier` CLI from the current source tree.

If you only want the published Python package and do not need this repository checkout,
you can instead use:

```bash
pip install symWannier
```

The Python install is enough for:

- `symwannier expand prefix`
- `symwannier wannierize ... prefix`

It is not enough for the full Quantum ESPRESSO workflow starting from NSCF output.
In particular, `pip install symWannier` does not install or build:

- `wannier90.x`
- Quantum ESPRESSO
- the modified `pw2wannier90.x` source bundled in this repository

If you want to generate `prefix.isym`, `prefix.immn`, `prefix.iamn`, and `prefix.ieig`
yourself, you still need a Quantum ESPRESSO build that includes the modified
`pw2wannier90.x`.

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

Run pw2wannier90.x with `irr_bz = .true.` to compute `prefix.immn`, `prefix.iamn`,
and `prefix.ieig`. Symmetry infomation is stored in `prefix.isym`.

If you want to use projectability-based disentanglement from atomic projectors,
`atom_proj = .true.` can be combined with `irr_bz = .true.` in the modified
`pw2wannier90.x` source bundled in this repository. In that case, the generated
`prefix.iamn` can still be expanded and used by the Python-side `wannierize`
workflow.

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

Once `prefix.immn`, `prefix.iamn`, `prefix.ieig`, and `prefix.isym` are obtained,
you can construct Wannier functions as follows.
```
symwannier wannierize -s prefix
```
To construct symmetry-adapted Wannier functions, use the `-S` option instead of `-s`.
```
symwannier wannierize -S prefix
```

If you want to define disentanglement windows from AMN projectability rather than
from the energy window in `prefix.win`, add the `-P` option.

Typical combinations are:

```bash
symwannier wannierize -P -s prefix
symwannier wannierize -P -S prefix
```


## Paper

For more information, please see
[T. Koretsune, Comp. Phys. Comm. 285 108645 (2023).](https://doi.org/10.1016/j.cpc.2022.108645)

We hope that you cite this reference when you publish the results using this code.
