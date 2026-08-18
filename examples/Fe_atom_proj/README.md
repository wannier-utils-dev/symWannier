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

The projectability route uses `dis_proj_min=0.01` and `dis_proj_max=0.90`
together with the original energy windows. As in Wannier90, the frozen
subspace is the union of states in the inner energy window and states above
`dis_proj_max`, restricted to the outer energy window. On the reference data,
Wannier90 and SymWannier give total spreads of 17.333924328 and
17.333923686 Angstrom^2, respectively. The bundled pytest also checks the
frozen-window size (12--16 states), outer-window size (18--28 states), and
the SymWannier spread. Since `num_iter=0`, these values test disentanglement
and the initial gauge rather than a subsequent MLWF minimization.
