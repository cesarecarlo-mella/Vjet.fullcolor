# C3F0Z0 HA3: Trident functions replaced by their corner series

## Inputs

- **Expression:** `zj.fullcolor/C3F0Z0_HA3.wl`, which has the form
  HA3 = Σ R_k(t,u) Trident[k] + I·Pi·Σ R'_k(t,u) Trident[k].
  - There are 1079 + 215 terms; every coefficient R is rational in t, u, with t = y and u = 1−x−y.
  - The Trident functions are real (phase-normalised).
- **Series:** `mathematica-results/trident/trident_selected_cK_ord20.m`.
  - They contain the 1079 Trident functions needed here, at corners c1 (x,y), c2 (y,z) and c3 (x,z).
  - Each series is rectangular, to order 20 in both corner variables.
- **Exact reference:** HAG["C3F0Z0",3], the GPL form, evaluated with `gpl_eval.py`.
  - It is accurate to about 1e-16.
  - The values are those of `HA_series_Nc3_vs_exact_digits/results/digits_order18.json`.

## 1. Exact check at corner 1, order 18 (`exact_check_c1.py`)

**Why the series cannot be expanded directly.** The coefficients R_k contain a spurious
denominator (x+y)^3, since 1−u = x+y. The sum Σ R_k T_k therefore cannot be expanded as a
double series directly, which is also why Series[] in the notebook aborts.

**What is compared instead.** Multiply by the common denominator
D = t²·u·(u−1)³·(t−1)·(t+u−1)·(t+u)³. This gives

    D·HA3 = Σ P_k T_k,   with P_k = D·R_k polynomial.

**Range of validity.**
- **Right-hand side, Σ P_k T_k:** the coefficient of x^a y^b only needs T coefficients with
  powers ≤ (a,b), so it is exact for a, b ≤ 20.
- **Left-hand side, D·HAser:** HAser["C3F0Z0",3] is the order-18 expansion of the GPL form, so
  this side is exact for a, b ≤ 18.

**How it is compared.** All coefficients with a, b ≤ 18 are compared, separately for every
monomial Log[x]^i Log[y]^j Pi^p Zeta3^q Zeta5^r I^s, modulo two 31-bit primes.

**Result: identical.** All 97 log/constant structures and 9248 nonzero coefficients agree, with no
mismatch, and every Log[x] term cancels. Log: `results/exact_check_c1.log`.

## 2. Numerical digit maps for c1, c2 and c3 (`eval_trident.py`, `plot_trident_digits.py`)

**Method.**
- Each Trident[k] is replaced by its cK series, truncated to orders a, b ≤ N, with N = 18 and N = 20.
- The prefactors R are evaluated exactly at each point.
- The sum is computed in long double on the same 406-point grid as before.
- Digits are min(16, −log10|S−E|/|E|), computed separately for Re and Im.

**Precision caveat.** The sum has large internal cancellations: Σ|R T| / |Σ R T| is about 2e5
(median) for Re, 1e3 for Im, and up to 3e7. Close to the corners this caps the real part at about
13–15 digits.

| N = 18, median digits (Re / Im) | c1 | c2 | c3 | best of 3 | HA3 series from GPLs (c1) |
|---|---|---|---|---|---|
| median | 4.1 / 4.8 | 7.3 / 8.1 | 4.3 / 5.1 | **10.6 / 11.4** | 3.1 / 4.2 |
| points ≥ 8 digits | 24% / 26% | 43% / 50% | 23% / 27% | **89% / 96%** | 19% / 24% |

With N = 20, the best of the three corners reaches a median of 11.6 / 12.4 digits, with 96% / 99% of
the points at ≥ 8 digits. The worst point has 7.0 digits (Re).

Full statistics: `results/digit_stats.txt`.

**Figures** (`figs/Nc3_HA3_trident_series_digits_N{18,20}`):
- columns: c1, c2, c3, best of the three corners, and (for N = 18 only) the HA3 series from GPLs at c1;
- rows: real part and imaginary part.

## Files

| file | purpose |
|---|---|
| `tparse.py` | parses `trident_selected_cK_ord20.m` into numpy arrays |
| `wparse.py` | parses the `.wl` expression |
| `eval_trident.py` | usage: `<corner> <N> <exact.json> <out.json>` |
| `exact_check_c1.py` | the exact check of section 1 |
| `plot_trident_digits.py` | the digit maps of section 2 |
