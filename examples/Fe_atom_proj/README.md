# Fe atom_proj

This example tests the QE 7.6 `atom_proj` interface with
`atom_proj_exclude`, `irr_bz=.true.`, noncollinear magnetism, and no SOC.
The NSCF and Wannier meshes are both 4 x 4 x 4.

Place `Fe.pbe-spn-rrkjus_psl.0.2.1.UPF` in this directory, set the
executable paths at the top of `test.sh`, and run:

```sh
./test.sh
```

The input excludes eight of the 26 UPF atomic projectors and retains 18
Wannier functions. The script follows the same command sequence as
`examples/Fe/test.sh`, then additionally runs SymWannier projectability
disentanglement with `-P -S`. Both the ordinary symmetry-aware
disentanglement and projectability-disentanglement routes passed in the
reference calculation.
