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
