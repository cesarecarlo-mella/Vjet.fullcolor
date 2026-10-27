# HA coefficients as GPLs in (x, y)

This is the procedure of `../convertLC.nb` applied to every file in
`../../remainder_denominator_bracket/` that is built only from letters 1–6. Each output file
`<layer>_HA<i>_G.m` defines `HAG["<layer>", i]`. Nothing here is compared with the existing
results.

## Conventions (identical to the notebook)

- **Letters:** `all.ints/a_tildes/letters.dict.m`. In the files used here only letters 1–6
  appear: l1 = log x, l2 = log y, l3 = log(1-x-y), l4 = log(1-x), l5 = log(1-y),
  l6 = log(x+y).
- **Ordering of the II:** dII[a1,...,an] = dlog l[an] II[a1,...,a(n-1)], so a1 is the innermost
  letter, as in the notebook's recursion.
- **Path and base point:** the path is (0,0) -> (0,y) -> (x,y), starting at (0,0) with tangential
  regularisation, i.e. log x and log y are dropped at the base point (G(...;0) = 0). This is what
  the notebook's two GIntegrate calls do: first along y at x = 0, then along x at fixed y.
- **Variables:** rational prefactors use u -> 1-x-y and t -> y, as in the notebook.
- **G definition:** `G[a1,...,an, z] = Int_0^z dt/(t-a1) G[a2,...,an, t]`, the PolyLogTools
  convention.
  - G of y has indices {0, 1}.
  - G of x has indices {0, 1, 1-y, -y}.
  - The result is directly in the fibration basis {x, y}.
- **Constants** are kept as they are in the source: Pi, I, Zeta[3], Zeta[5].

## How it is computed

`ii2g.py` applies Chen's formula for the composed path:

    II[a1..an] = sum_k G(cx(an),...,cx(a(k+1)); x) * G(cy(ak),...,cy(a1); y)

Here cy and cx are the GPL indices of dlog l restricted to the two legs; a letter that vanishes
on a leg (l1 and l4 at x = 0, l2 and l5 at fixed y) kills that term. The result is the same
path-ordered integral that the notebook builds recursively. The notebook only sets
G[0,0] -> 0; here every G(...;0) is set to 0, which matters from weight 3 on.

The coefficients are collected per (constant, G(x), G(y)) monomial, then put over a common
denominator and factored.

## Check

`check_chen.py` compares the G form with a direct numerical path-ordered integral along the same
path for 18 random words of weight 1–3 at (x, y) = (0.31, 0.22). The maximum difference is
9e-15.

## Coverage

| layer  | HA | letters | done |
|--------|----|---------|------|
| C1F2Z0 (the notebook case) | 1–6 | 1–6 | yes |
| C9F2Z0, C1F1Z1, C9F1Z1, C2F1Z0, C3F0Z0, C8F1Z0 | 1–6 | 1–6 | yes |
| C0F1Z0, C8F0Z1 | 1–6 | + 15, 16 | no |
| C0F0Z1, C2F0Z1 | 1–6 | + 11–16 | no |
| C1F0Z0, C7F0Z0, C9F0Z0 | 1–6 | + 7, 8 (square roots) | no |

The remaining layers are not done, for these reasons:
- **Quadratic letters 11–16.** These give G indices that are roots in y, such as
  (1 ± Sqrt[1-4y])/2, 1 ± Sqrt[y] and -y ± Sqrt[y]. Some of these lie on the integration path
  inside the region (for example x² - x + y changes sign), so an i0 prescription has to be
  chosen.
- **Letters 7 and 8.** These involve Sqrt[x y (1-x-y)] and are not GPLs in x at fixed y without
  rationalising first.

## Files

- `run_all.log` (laptop) and `run_cloud.log` (cloud): timings and the number of G monomials per
  file.
- `run_all.sh` / `run_rest.sh`: the runners used.
