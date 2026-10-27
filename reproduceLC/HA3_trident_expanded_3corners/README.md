# C3F0Z0 HA3: full expansion around the three corners (Trident series + rational functions)

The starting point is `zj.fullcolor/C3F0Z0_HA3.wl`:

    HA3 = Σ R_k(t,u) Trident[k] + I·Pi·Σ R'_k(t,u) Trident[k]

Here t = y, u = z = 1−x−y. Each Trident function is replaced by its corner series
(`trident_selected_cK_ord20.m`). The rational prefactors R_k are **also expanded**, so the results
are pure series in the corner variables.

## Contents of `series/`

There are three objects, one per corner, plus their Bernoulli forms. All coefficients are exact.

| corner | variables (p, q) | native file | Bernoulli file (bp = −Log[1−p]) | truncation | lowest powers |
|---|---|---|---|---|---|
| c1 | x, y | `C3F0Z0_HA3_expanded_c1.m` → `HA3exp["C3F0Z0",3,"c1"]` | `C3F0Z0_HA3_bernoulli_c1.m` → `HA3bern[…,"c1"]` (bx, by) | total degree ≤ 14 | y^−1 |
| c2 | y, z | `…_expanded_c2.m` | `…_bernoulli_c2.m` (by, bz) | total degree ≤ 14 | y^−1, z^−1 |
| c3 | x, z | `…_expanded_c3.m` | `…_bernoulli_c3.m` (bx, bz) | total degree ≤ 17 | z^−1 |

- Logarithms appear as Log[p] and Log[q]; they are kept unexpanded in the Bernoulli form.
- Constants appear as Pi, Zeta[3], Zeta[5] and I.
- Each file is written as `name = ( … );` so that `Get` reads the whole expression.

Example:

```
Get["series/C3F0Z0_HA3_expanded_c2.m"];
HA3exp["C3F0Z0", 3, "c2"]
```

## Why the truncation is in total degree, and why it is 14, 14 and 17

### The spurious factor

At each corner, the common denominator D of the R_k factorises as D = M·B·U:

- M is a monomial;
- U is a unit (it does not vanish at the corner);
- B is a homogeneous spurious factor that is not a monomial.

| corner | M | B | deg M + deg B |
|---|---|---|---|
| c1 | x·y² | (x+y)³ | 6 |
| c2 | y²·z | (y+z)³ | 6 |
| c3 | x·z | (x+z) | 3 |

1/B has no double expansion around the corner, so the product cannot be expanded term by term.

### How the expansion is built (`expand_corner.py`, exact arithmetic with FLINT)

1. Compute Q = Σ (D·R_k) T_k. Here D·R_k is a polynomial and T_k is the order-20 trident series, so
   Q is exact up to total degree 20.
2. Divide by U as a series.
3. Divide by B **degree by degree**. Each homogeneous part must be exactly divisible by B. It is, at
   every corner and for every log/constant structure: this shows the spurious pole cancels.
4. Shift by 1/M.

### What this costs

The result is exact up to total degree 20 − deg M − deg B, which gives 14, 14 and 17 for c1, c2 and
c3. Higher terms are **not fixed** by the order-20 trident series.

The reason is that multiplication by B is not invertible on truncated series. For example, at c1 a
rectangle of order 18 in (x, y) would need trident series up to total degree 42. With trident
series of order 24 (or 21 at c3), the expansion would reach total degree 18.

## Checks

- **c1 against the GPL series:** the expansion agrees with `HAser["C3F0Z0",3]` (the GPL-based
  series in `HA_series_Nc3_xy0_order18`) restricted to total degree ≤ 14. All 5487 coefficients are
  identical.
- **Divisibility:** there were no failures in step 3 at any corner. See `results/run.log`.

## Digit maps against the exact GPLs (`figs/Nc3_HA3_expanded_3corners_digits.*`)

- Grid: the same 406 points as before.
- Digits: min(16, −log10|S−E|/|E|), computed separately for Re and Im.
- Exact reference E: HAG["C3F0Z0",3] evaluated with `gpl_eval.py`.

Median digits:

| | c1 | c2 | c3 | best of 3 |
|---|---|---|---|---|
| native Re / Im | 1.9 / 2.6 | 4.5 / 5.6 | 2.2 / 2.7 | 6.7 / 7.9 |
| Bernoulli Re / Im | 3.1 / 4.0 | 6.0 / 7.0 | 3.6 / 4.3 | **8.9 / 9.9** |

Taking the best of the three corners with the Bernoulli form, the worst point has 4.8 digits (Re)
and 5.2 digits (Im). Full table: `results/digit_stats.txt`.

## Scripts (`scripts/`)

| script | what it does |
|---|---|
| `expand_corner.py <cK> <outdir>` | builds the native series |
| `bern_eval.py <cK> <dir> <exact.json> [HAser.m]` | builds the Bernoulli form, runs the c1 check, and evaluates on the grid |
| `plot_expanded_digits.py <dir> <outdir>` | draws the figure |

Inputs needed: `C3F0Z0_HA3.wl`, `sel_cK.m` (the trident series) and the exact-value json.
