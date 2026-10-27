# Nc^3 layer: series around (0,0) vs exact GPLs — digit maps

## What is compared
- **exact**: `HAG["C3F0Z0", i]` (i = 1..6) from `HA_GPL/`, i.e. the GPL form built from the II's.
  The G's are evaluated with `gpl_eval.py`: Gauss–Legendre panels with spectral cumulative
  integration in long double.
  - x-words use the indices {0, 1, 1−y, −y}; y-words use {0, 1}.
  - Cross-checked against GiNaC (30 digits) and against the published 2024 ancillary at 4 points.
    The largest relative difference was between 1e-19 and 7e-16.
  - So "16 digits" in the maps means **at least about 15** digits.
- **series**: `HAser["C3F0Z0", i]` from `HA_series_Nc3_xy0_order{15,18}/`, evaluated in two ways:
  - *native*: as written, in powers of x and y (with Log[y]);
  - *Bernoulli*: re-expanded in bx = −log(1−x), by = −log(1−y) at the same truncation.

## Grid and digit definition
- Grid points: x = i/30, y = j/30 with i + j ≤ 29, giving 406 points inside the triangle.
- Digits are computed separately for Re and Im: `d = min(16, −log10(|s − e| / |e|))`.
  - A cell is chopped (blank) when |e| < 1e-4.
  - d = 16 when |s − e| < 1e-16 |e|.
- The summary figure takes, at each point, the **minimum over the 6 HA**.

## Results (minimum over HA)
| N | expansion | median Re | median Im | share of points with ≥8 digits (Re / Im) |
|---|---|---|---|---|
| 15 | native    | 2.0 | 2.8 | 9% / 13% |
| 15 | Bernoulli | 4.8 | 5.2 | 29% / 30% |
| 18 | native    | 2.6 | 3.4 | 14% / 18% |
| 18 | Bernoulli | 6.0 | 6.3 | 38% / 40% |

Every case reaches 16 digits near (0,0). Accuracy falls below 1 digit as x + y → 1, the hypotenuse.

## Files
- `series_vs_exact.py <series_dir> <N> <decay.m> <outdir> [K] [ncores]` writes `results/digits_order{N}.json`.
  Each row holds x, y, the exact values, and the native and Bernoulli digits as [Re, Im] for each HA.
- `plot_series_digits.py` produces the figures in `figs/` (in the style of `plot_w6_paper.py`):
  - `Nc3_HA_digits_{native,bernoulli}_N{15,18}`: Re row and Im row × 6 HA;
  - `Nc3_HA_digits_min_summary`: rows N = 15, 18; columns native Re/Im, Bernoulli Re/Im.
