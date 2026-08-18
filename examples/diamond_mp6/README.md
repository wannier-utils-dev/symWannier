# diamond_mp6 — k-point snapping demo (carbon diamond, 6×6×6)

A small carbon-diamond example on a **6×6×6** k-mesh. Unlike the 8×8×8 grids in
`examples/Cu` and `examples/Fe`, `6` has the factor 3, so the mesh k-points
`k = i/6` are **not** exactly representable in the 8-decimal format that
wannier90 writes into the `.nnkp` file. This makes the example show the effect
of the `--snap-kp` option.

## Quick check (symWannier only, no DFT tools needed)

The IBZ inputs (`diamond_mp6.i*`, `diamond_mp6.nnkp`, `diamond_mp6.win`) are
included, so you can verify the fix with only symWannier:

```
python verify.py
```

Expected output:

```
System: 6x6x6 mesh = 216 k-points, 4 bands -> 4 Wannier functions (no disentanglement)

(A) DFT bands vs Wannier bands   max | eps_Wannier(k) - eps_DFT(k) |
    reference: QE eigenvalues from .ieig (external)
    max over 216 mesh k-points and 4 bands, in eV:
      --snap-kp false  (H(R) from 8-digit .nnkp k) : 2.34e-06 eV
      --snap-kp true   (H(R) from exact i/mp_grid) : 1.98e-10 eV
    floor of this metric, with NO Fourier transform at all:
      max| eig H(k) - eps_DFT | = 1.98e-10 eV   (MLWF gauge U, ||U^dag U - 1|| = 6.9e-12)

(B) Fourier self-consistency   max | H(R)->H(k) - H(k) |
    reference: symWannier's own H(k) (internal; the gauge floor cancels)
    max over 216 mesh k-points and 4x4 matrix elements, in eV:
      --snap-kp false  (H(R) from 8-digit .nnkp k) : 1.76e-06 eV
      --snap-kp true   (H(R) from exact i/mp_grid) : 3.74e-14 eV

Underlying cause: the written H(R) coefficients themselves differ by
  max| H(R)_false - H(R)_true | = 2.07e-07 eV
```

(the last digits vary a little with the BLAS/platform; the point is ~1e-6
before the fix versus the respective floors after it)

**What the two numbers mean.** They use different references and answer
different questions:

- **(A) is the error a user of the model sees**: the interpolated Wannier bands
  against the DFT eigenvalues from QE (an *external* reference). Its floor is
  set not by the Fourier transform but by the unitarity of the MLWF gauge
  matrix `U` — `H(k)` is built as `U†·diag(ε_DFT)·U`, so a `U` that is unitary
  only to ~7e-12 already shifts the eigenvalues by ~2e-10 eV before any
  transform. The snapped value lands exactly on that pre-existing floor, so
  nothing of the truncation is left.
- **(B) isolates the bug**: comparing against symWannier's own `H(k)` (an
  *internal* reference) puts the same `H` on both sides, so the gauge floor
  cancels and only the Fourier-transform error remains. This is a
  self-consistency check, *not* a comparison against DFT.

**Where the error comes from.** symWannier builds `H(R)` (written to
`hr.dat`/`tb.dat`) as `H(R) = (1/Nk) Σ_k exp(-2πi k·R) H(k)`. Without
`--snap-kp` those phases use the 8-digit `.nnkp` k (`1/6` stored as
`0.16666667`, ~3e-9 low), which corrupts the `H(R)` coefficients by ~2e-7.
The *backward* transform (band interpolation `H(R)→H(k)`) always evaluates at
double-precision k and never re-reads the `.nnkp`, so it is not where precision
is lost — the verification above therefore builds `H(R)` both ways but always
evaluates back at the exact mesh k.

Equivalently, run the Wannierization both ways and compare the Hamiltonian:

```
python -m symwannier.wannierize -S -H --snap-kp true  diamond_mp6   # exact rationals (prints a note)
python -m symwannier.wannierize -S -H --snap-kp false diamond_mp6   # finite-digit .nnkp k
```

`--snap-kp` defaults to `true`; pass `--snap-kp false` to keep the raw
`.nnkp` values.

## Reproducing the inputs from scratch (needs QE + wannier90)

`test.sh` runs the full pipeline (`pw.x` → `wannier90.x -pp` →
`pw2wannier90.x` with `irr_bz=.true.`) to regenerate the `.i*` files. Edit the
tool paths at the top and provide the `C.pz-vbc.UPF` pseudopotential, then:

```
sh test.sh
```
