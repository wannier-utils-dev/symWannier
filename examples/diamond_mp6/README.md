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
Fourier round-trip error  max| H(R)->H(k) - H(k) |  at the mesh points:
  --snap-kp false  (finite-digit .nnkp k) : ~2e-06
  --snap-kp true   (exact i/mp_grid k)    : ~5e-14

(1) deliverable H(R) coefficients (written to hr.dat/tb.dat), no interpolation:
      max| H(R)_false - H(R)_true | = ~2e-07

(2) interpolate H(R) back at the EXACT double-precision mesh k, vs true H(k):
      --snap-kp false  (H(R) from 8-digit .nnkp k) : ~2e-06
      --snap-kp true   (H(R) from exact i/mp_grid ) : ~4e-14
```

(the last digits vary a little with the BLAS/platform; the point is ~1e-6 vs
machine precision)

**Why this is the right test.** A Wannier tight-binding model reproduces the
ab-initio `H(k)` exactly at the k-points of the Wannierization mesh. So
transforming the model's `H(k)` to `H(R)` and back must recover `H(k)` at those
points to machine precision. With the finite-digit `.nnkp` k-points (`1/6`
stored as `0.16666667`, ~3e-9 low) the phases `2π·k·R` summed over the cell
amplify the error to ~1e-6; snapping k to the exact rational `i/6` (default)
restores machine precision (~1e-13).

**A note on what this measures (lines (1)/(2)).** The round-trip above uses the
*same* k for `H(k)→H(R)` and `H(R)→H(k)`, which slightly overstates the role of
the backward transform. The error actually lives entirely in the forward build
of `H(R)`: the coefficients written to `hr.dat`/`tb.dat` are themselves wrong by
~2e-7 (line (1)), independent of any interpolation. Wannier band interpolation
(`H(R)→H(k)`) always evaluates at double-precision k and never re-reads the
`.nnkp`, so it does not depend on `.nnkp` precision — yet interpolating the
corrupted `H(R)` even at the *exact* mesh k still misses the true `H(k)` by ~1e-6
(line (2), `--snap-kp false`). Snapping fixes the `H(R)` construction, so both
metrics reach machine precision.

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
