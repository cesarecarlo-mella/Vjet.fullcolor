# generators_deq -- closed DE and numerical evaluation of the generators and of all 84 finite remainders

The generators T[1..579] of `../../Zgqq_finite_remainders/Generators.wl` are pure combinations of iterated integrals
II[a1..an] over the 20-letter alphabet, weights 0-6. They are everything the finite remainders need: each of the 84
remainders (`alpha1-3` = HA1-3, `gamma1-3` = HA4-6, 14 colour layers each) is a sum of
q * RationalFunctionBasis[r] * T[g] * TranscendentalMultipliers[m].

- One multiplier, #13, is I/r with r = +Sqrt[x y z] (decided 2026-10-28).
- It multiplies exactly the 254 rows whose generators are odd under r -> -r (an odd number of the letters 7, 8). They
  occur only in the layers Nc^1, Nc^-1 and Nc^-3 with Nf^0 NfZ^0.
- The old denominator-bracket files had plain I there. That sector has no independent numerical check.
- `remainder_eval` evaluates any multiplier of the form constant * r^k.

- T[1..60] are the H[1..60] of `../Generators.wl` (identical word by word), and T[61] = 1.
- `Alphabet20.wl` (variables t = y, u = z = 1 - x - y) is letter by letter the alphabet used here, in (x, y)
  (`data/trident_definition/alphabet.wl`: letter 1 = x, 2 = y, 3 = 1 - x - y).
- Only 14 letters appear: 1-8 and 11-16. Letters 9, 10, 17-20 never occur.
- d II[a1..an] = dlog W_an II[a1..a(n-1)].
- Base point c1, shuffle regularised: II[] = 1 and every non-empty II -> 0.
- All 579 are integrable. 80 of them have first entries 1 - x or 1 - y (letters 4, 5).

## System

Close the generators under d: at each weight, the generators plus the strips (last letter removed) of the weight
above.

| weight | 0 | 1 | 2 | 3 | 4 | 5 | 6 | total |
|---|---|---|---|---|---|---|---|---|
| generators | 1 | 5 | 24 | 75 | 203 | 247 | 24 | 579 |
| linearly independent | 1 | 5 | 24 | 75 | 113 | 245 | 24 | 487 |
| extra functions to close | | | | | | 16 | | 16 |
| **system** | 1 | 5 | 24 | 75 | 113 | 261 | 24 | **503** |

- **Linear dependencies.** 92 generators are exact rational combinations of the others of the same weight: 90 at
  weight 4, and T[570], T[571] at weight 5. The rank is the same modulo two different primes. They need no slot of
  their own: they are read off at the end, like all the others. `out/dependent.wl` gives them (and T[61] = One)
  as exact combinations of the basis.
- **Closure.** Weights 1-4 are already closed: the derivatives of the generators stay in their span. Only weight 5
  needs 16 extra functions, from the derivatives of the weight-6 generators.
- The old 302-function system of the 60 generators lies inside this one, at every weight.

**Exact DE:** dF = Sum_k dlog W_k M[k] . F, with F the 503 functions (`exact_deq.py`):
- `out/basis.wl`: the functions F in II[...], weight 6 -> 0. At each weight the independent generators come first
  (with their names), then the extra functions F[5, j]; One = 1.
- `out/deq.wl`: labels and M[1..20] as SparseArrays, 130914 non-zero entries. Closed and checked word by word.
- The rationals are large (up to about 120 digits): the strips are expressed in the 113 independent weight-4
  generators, whose integer coefficients are large.

**Numerical gauge:** an orthonormal basis of each weight in word space. The letter blocks are then the strip maps,
with every entry |A| <= 1.

## Numerical evaluation (method of `INTEGRALS/4P1M/3L1M/trident_deq`, Czakon-Tancredi)

`cpp/tridents_bs.cpp` is the trident_deq evaluator, unchanged:
- Bulirsch-Stoer on the complex-deformed straight path from P0 = (0.3, 0.3), deformation (0.1, 0.2);
- r = Sqrt[xyz] is integrated along;
- the generators are read off at the end point.

The data are in `out/generators_bs.bin` (TRIDBS01 format, NT = 579, generator i in the slot of its weight).

```sh
make                                              # builds cpp/tridents_bs
python3 gen_eval.py 0.25 0.35                     # the 579 generators at one point
python3 remainder_eval.py 0.25 0.35               # the 84 finite remainders at one point
```
```python
from remainder_eval import evaluate_remainders, KEYS, LABELS
R, info = evaluate_remainders([(0.3, 0.2), (0.1, 0.7)])    # R.shape == (2, 84)
# column j: KEYS[j] = ('alpha1', 'Nc^0*Nf^0*NfZ^1'), LABELS[j] = 'C0F0Z1_HA1' (the old file names)
```

Points are (x, y) with letter 1 = x, letter 2 = y. The remainder files are in t = y, u = 1 - x - y.

**Speed:** about 0.06 s per point, for all 579 generators and the 84 remainders; 0.12 s with the error estimate.

**Error estimate** (`evaluate_remainders(..., error=True)`). Every point is evaluated a second time. That
evaluation starts from independent start values and runs on a different path; for an averaged point it also uses half
the stencil step.
- **Independent start values.** The start values at P0 come from the y-first regularisation instead of the x-first
  one. The second data file differs from the first only in its first 72588 bytes, the start values. The system blocks
  are identical, so only `out/start_y.head` is stored. `gen_eval.second_data()` builds the full file once, next to
  the first one or in the temp directory.
- **Why both changes are needed.** The start values are an error shared by every evaluation from one file. Near the
  edges the rational functions amplify it: it was the main error there. A second path from the same file missed it,
  by up to a factor 10.
- **The estimate.** The reported error is `ERR_FACTOR` = 10 times the difference of the two evaluations. Against
  exact values the actual error was below it for 100% of the values (Nc^3 layer, 406 grid points) and for 99% (30
  points near the corner x, y -> 0, against the exact order-18 series of reproduceLC). The raw difference is
  typically 1-3 times smaller than the actual error.

**Near the edges** the rational functions cancel against large logarithms, and double precision loses digits there.
The poles go up to 1/x^3, 1/y^3 and 1/z^3, and up to (x+y)^-5 etc. at the corners. Tightening the tolerance does not
help: the limit is roundoff. Typical worst relative error over the 84 remainders, by distance d from the edge:

| d | x -> 0 | y -> 0 or z -> 0 | corners |
|---|---|---|---|
| 0.1 | 1e-10 | 1e-9 | 1e-9 |
| 0.03 | 1e-9 | 1e-7 | 1e-5 |
| 0.01 | 1e-9 | 1e-5 | 1e-3 |
| 0.003 | 1e-8 | 1e-4 | O(1) |

These regions need the corner series.

**Accuracy** over the triangle at the default tolerance 1e-10 (`gen_eval_check.py`):
- relative error median 1e-13 (worst 5e-9, weight 6, at random points; 7e-9 near the edges);
- **precision floor about 1e-10:** some generators have integer coefficients up to about 5e4 spread over hundreds of
  words, and an O(1) value then comes from large cancellations. That loses about 5 digits in double precision. It
  shows up as the 1.5e-10 difference between the x-first and y-first start values at P0, at any series order.

**Spurious curves.** Letters 11-16 vanish on curves inside the triangle (e.g. letter 11: x = y - y^2). The
functions are smooth there, but the DE has poles. Two different things can happen:
- **The path goes through a zero.** The complex zeros of a letter can lie exactly on the deformed path; for
  example, for (7/30, 19/30) letter 11 vanishes on the default path at t = 3/4, and the stepper fails. Along the
  path every singular polynomial (the letters in use, and x z + y, x y + z for the odd letters) has degree <= 4 in
  t. So `gen_eval.choose_deformation` computes its roots before integrating and keeps the default deformation only
  if they stay at least 1e-3 away from 0 < t < 1; otherwise it takes the candidate with the largest margin. This
  switches about 1% of random points, and 2 of the 406 grid points; it costs about 3 ms per point. Different
  deformations agree to 3e-11. `info[p]['delta']` records the deformation used.
- **The end point is on a zero curve** (the phase-space point itself, e.g. (0.21, 0.7) on letter 11, or
  (0.24, 0.36) on letter 14). No deformation helps: the direct integration fails at several such points. Close to
  the curves of letters 15 and 16, the rational functions of the remainders also have poles. They cancel in the
  remainder, but the direct evaluation loses digits: 5e-10 at distance 1e-6, 7e-8 at 1e-8. These letters are
  spurious and the remainders are smooth there. So a point close to a curve of letters 11-16 is replaced by the
  symmetric average across the curve. All lengths are relative to rho = min(x, y, z), because all spurious curves
  run into the corners:
  - **When.** The point is averaged if its distance to the curve is below `NEAR` * rho.
    - `NEAR` = 1e-8 for letters 11-14. These are poles of the DE only, and the direct integration works until
      almost exact hits.
    - `NEAR` = 1e-5 for letters 15, 16. These are also poles of the rational functions.
  - **The stencil.** The points P +- h n, P +- 2h n, P +- 3h n, with n the unit normal of the curve (the gradient of
    the letter). The step is h = min(`H` = 1e-3, rho / 100).
  - **The combination.** f(P) = (15 S1 - 6 S2 + S3) / 10 + O(h^6), with S_k = (f(P + k h n) + f(P - k h n)) / 2.
    The O(h^6) term is about 36 (h / rho)^6 <= 4e-11.
  - **A mistake that the absolute thresholds made.** A first version used an absolute distance of 1e-5 and
    h = rho / 10. It averaged diagonal points near the corner x, y -> 0, which lie within 1e-5 of the letter-11
    curve. With h = rho / 10 the O(h^6) error there was 1e-5.

  The stencil points are ordinary points about 1e-3 away from the curve. Each one is integrated on its own
  complexified path, with the path check above. The normal is used, not the path direction, because the path
  direction can be nearly tangent to the curve. The plain 2-point average would only be good to (h/rho)^2 ~ 1e-4.
  - `gen_eval.evaluate` averages the generators.
  - `remainder_eval.evaluate_remainders` averages the assembled remainders. This covers the poles of the DE and
    those of the rational functions at once.
  - `info[p]['averaged']` records the letter, the distance and h. An averaged point costs 6 integrations.
  - Random points essentially never need it. Grids with simple rational points can: none of the 406 grid points do.

Start values (`iisys.start_values`):
- log-power series from c1, along x = 0 and then in Sqrt[x] at y = 1/2, giving the values at (0.02, 0.5);
- carried once to P0 with a tight tolerance (`gen_run.py`).
- At y = 1/2, 60 terms are enough (the same values as 140 terms). At y = 0.7 they are not (6e-7); 100 are.

## Checks

| check | result |
|---|---|
| **84 remainders vs exact values**: leading-colour layer C3F0Z0, HA1-6, at the 406 grid points of `reproduceLC/HA_series_Nc3_vs_exact_digits` (GPL form, >= 15 digits) (`remainder_check.py`, tolerance 1e-12) | relative difference median 1-5e-13, worst 7e-10 (corner points next to the edges, \|HA\| ~ 3e3) |
| **error estimate vs actual error** (`remainder_check.py`; 30 points near the corner against the exact order-18 series) | actual error below 10 x the raw difference for 100% (grid) and 99% (corner) of the values |
| **end points on the curves of letters 11-16**, 9 points (`spurious_check.py`): averaged remainders, step H vs H/2 | 4e-11 to 3.5e-10 (8e-10 at (0.21, 0.7), near the hypotenuse) |
| same 9 points moved 1e-4 off the curve: direct remainders vs the average forced there | 3e-11 to 4e-10 (2.6e-9 at (0.21, 0.7)) |
| tolerance 1e-10 vs 1e-13, and the two sides of the path deformation, 60 random + 9 near-edge points | median 1e-13, worst 7e-9 |
| x -> 0 first vs y -> 0 first regularisation at c1, and both deformation sides, at P0 | basis functions agree to 7e-10 absolute, generators to 1.5e-10 relative |
| direct c1 series (100 terms) near the edges vs integration from P0 | 1.3e-11, 4.2e-10, 2.2e-9 |
| 24 tridents through the same pipeline (`trident_check.py`, independent of the generators) vs the exact corner series, 21 points near c1, c2, c3 | worst 7e-13 (eps^1), 5e-11 (eps^5), 3e-8 (eps^6), as good as trident_deq |

The representation itself (`Zgqq_finite_remainders`) was checked exactly against the old denominator-bracket files:
all 84 layers, every (II word, constant) coefficient, at random rational points. That was before multiplier #13
became I/r.

The shippable package is assembled from these files by `../../package/build_dist.py`.

## Files

| file | content |
|---|---|
| `iisys.py` | closure, orthonormal system, 1-d log-power series start values, TRIDBS01 export (works for any set of pure II combinations) |
| `gen_run.py` | build: system, both regularisations, carry to P0; writes `out/generators_bs.bin`, `out/start_y.head`, `out/meta.json`, `out/system.pkl` |
| `gen_eval.py` | python interface for the generators: choice of path deformation, average across spurious curves (`plan`, `collapse`) |
| `remainder_eval.py` | the 84 finite remainders from the generators (averaged as a whole near spurious curves) |
| `exact_deq.py` | the exact DE (`out/basis.wl`, `out/deq.wl`, `out/dependent.wl`) |
| `gen_eval_check.py`, `remainder_check.py`, `spurious_check.py`, `trident_check.py` | the checks above (`trident_check.py` reads `2610.anc/trident_basis`) |
| `out/meta.json` | names, weights and system dimensions (all that `gen_eval.py` needs) |
| `out/start_y.head` | the head (start values at P0) of the data file built from the y-first regularisation, for the error estimate |
| `out/generators_bs_y.bin` | that full file, built on first use by `gen_eval.second_data()` (not in git) |
| `out/system.pkl` | the numerical system: not in git (60 MB), rebuilt by `gen_run.py` in 20 s; needed by `exact_deq.py` and the checks |

To rebuild everything:

```sh
python3 gen_run.py ../../Zgqq_finite_remainders/Generators.wl out        # about 30 s
python3 exact_deq.py out/system.pkl ../../Zgqq_finite_remainders/Generators.wl out
```
