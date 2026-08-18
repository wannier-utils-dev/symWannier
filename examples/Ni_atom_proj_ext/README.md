# Ni atom_proj_ext + atom_proj_exclude

This QE 7.6 regression example uses a 4 x 4 x 4 NSCF/Wannier mesh and tests
`atom_proj_ext=.true.`, `atom_proj_exclude`, and `irr_bz=.true.`
together.

Place `Ni.pbe-n-rrkjus_psl.0.1.UPF` in this directory. Alternatively,
`NI_UPF=/path/to/Ni.UPF` may be used by `make_atom_proj.py`, but the QE
inputs still expect the pseudopotential in the current directory.

The generator reads the s, p, and d `PP_CHI` radial wavefunctions and writes
`atom_proj/Ni.dat`. Each radial channel is duplicated, producing 18
external projectors:

```text
s1, s2, p1(3), p2(3), d1(5), d2(5)
```

`pw2wan.in` excludes projector indices 2, 6--8, and 14--18, retaining the
complete `s1 + p1(3) + d1(5)` multiplets (nine Wannier functions). This
specifically tests external-projector metadata compaction with
`proj_excl_map`.

Set the executable paths at the top of `test.sh`, as in the Fe and Cu
examples, and run:

```sh
./test.sh
```

The tested QE source is
`src/pw2wannier90/qe7.6/pw2wannier90.f90`. The reference calculation
completes AMN and IBZ MMN and reaches `JOB DONE`. SymWannier expansion gives
`symmetrize Gk diff1 = 1.61420e-3`; the ordinary UPF Ni case gives
`4.41654e-3`. Wannier90 and SymWannier energy disentanglement give total
spreads 6.673898109 and 6.67399836 Angstrom^2, respectively.

Projectability disentanglement is intentionally omitted from this script.
For this deliberately duplicated external radial set, its automatically
selected subspace is undercomplete; that optional failure is independent of
the external-projector metadata fix.
