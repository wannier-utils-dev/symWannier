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

Calculate Mmn, Amn and Eig in the full BZ using ```symwannier expand```.
```
symwannier expand prefix
```

Two further files can be produced from the same IBZ data, neither of which
pw2wannier90.x writes together with ```irr_bz```.

```
symwannier expand --dmn prefix
```
also writes prefix.dmn, the symmetry file wannier90 reads with
```site_symmetry = .true.```. The file cannot represent an antiunitary
operation, so it is not written for every system; where that happens the command
says so rather than writing something wrong.

```
symwannier expand --spn prefix
```
also expands prefix.ispn to prefix.spn, the spin matrix elements postw90 uses.
Run pw2wannier90.x with ```write_spn = .true.``` and ```spn_formatted = .true.```
to obtain prefix.ispn.

Run wannier90 as usual.
```
wannier90.x prefix
```

There are sample input and script files in `examples/Cu`, `examples/Fe`,
`examples/Fe_atom_proj`, `examples/Ni_atom_proj`, and
`examples/Ni_atom_proj_ext`.


## Symmetry adapted Wannier functions

Once prefix.immn, prefix.iamn, prefix.ieig, and prefix.isym are obtained, you can construct Wannier functions as follows.
```
symwannier wannierize -s prefix
```
To construct symmetry-adapted Wannier functions, use the "-S" option instead of "-s".
```
symwannier wannierize -S prefix
```

The Wannier center involves ```Im ln M_nn```, which is multivalued. The branch is
taken closest to the projection centers, so that a projection written at a periodic
image far from the origin still localizes. wannier90 calls the same idea
```guiding_centres``` and has it off by default; here it is on unless the win file
says ```guiding_centres = .false.```, which restores the principal branch.

The spread is minimized by conjugate gradients with wannier90's parabolic line
search, as wannier90 does. ```--minimizer lbfgs``` uses L-BFGS instead, which
needs the fewest iterations, and ```--minimizer sd``` steepest descent, which is
slow where the spread is badly conditioned. With the symmetry information of
```-s``` and ```-S```, every step of either is kept to rotations that respect the
crystal symmetry. Without it the minimization is unconstrained, and where the
symmetric Wannier functions are a saddle point of the spread - the atom-centred
s, p, d orbitals of bcc Fe are one - L-BFGS leaves them for symmetry-broken
orbitals of lower spread, which wannier90 does not reach.

### Running on several MPI ranks

```symwannier wannierize``` distributes the k-points of the full mesh over MPI
ranks when started under a launcher:
```
pip install symWannier[mpi]        # mpi4py, built against the MPI of your mpirun
mpirun -n 8 symwannier wannierize -s prefix
```
Under a scheduler's launcher, give it the MPI plugin of the MPI that mpi4py is
built against, for example ```srun --mpi=pmix```. Started without one, the
processes know nothing of each other, symwannier cannot tell them from separate
serial runs, and each of them computes everything and writes the same files.

Each rank holds the overlap matrices of its own block of k-points only, which is
where the memory goes, so the run can be spread over nodes for memory as well as
for time. The results are those of a serial run up to the order in which the sums
over k are taken, which moves them by about 1e-12. Only rank 0 prints and writes
files. ```symwannier expand``` is not distributed and refuses to start on more
than one rank.

Prefer ranks to threads: start one rank per core with ```OMP_NUM_THREADS=1```.
The overlap matrices are small matrices in large batches, which BLAS threads
speed up little - for a 1000 k-point calculation, 8 threads made the run 4%
faster, while 8 ranks made it 4.5 times faster.


## Bundled pw2wannier90 interfaces

Patched `pw2wannier90.f90` sources are provided under `src/pw2wannier90`.
The `qe7.6` version supports the `atom_proj` workflow, including projector
exclusion, UPF projectors, external projectors, and IBZ expansion for
SymWannier. The 4 x 4 x 4 regression inputs are in `examples/Fe_atom_proj`,
`examples/Ni_atom_proj`, and `examples/Ni_atom_proj_ext`.

## Paper

For more information, please see
[T. Koretsune, Comp. Phys. Comm. 285 108645 (2023).](https://doi.org/10.1016/j.cpc.2022.108645)

We hope that you cite this reference when you publish the results using this code.

