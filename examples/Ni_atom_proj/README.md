# Ni atom_proj

This example tests the ordinary QE 7.6 `atom_proj` interface together with
`irr_bz=.true.`. The NSCF and Wannier meshes are both 4 x 4 x 4.

Place `Ni.pbe-n-rrkjus_psl.0.1.UPF` in this directory, set the executable
paths at the top of `test.sh`, and run:

```sh
./test.sh
```

The UPF supplies the complete s, p, and d multiplets used for nine Wannier
functions. The script follows the same command sequence as
`examples/Fe/test.sh`, then additionally runs SymWannier projectability
disentanglement with `-P -S`. Both the ordinary symmetry-aware
disentanglement and projectability-disentanglement routes passed in the
reference calculation.

## Verifying the unitarity of -S

`verify.py` reruns `-S` on the IBZ inputs bundled under `tests/inputs` (no QE or
wannier90 needed) with and without the nearest-unitary projection that follows
the final symmetrization, and with the outer window as given and cut below the
multiplet that `num_bands = 25` splits at k = (0,0,1/4):

```sh
python verify.py
```

It prints the unitarity of the final U, the point-group residual of
H(k) = U^dag diag(eps) U on the full mesh, the same residual of H(k)
rebuilt from the written hr.dat, and two other effects of U^dag U != 1: the
eigenvalues of H(k) changing with the energy zero, and their error from
ignoring the Wannier-function overlap U^dag U.
