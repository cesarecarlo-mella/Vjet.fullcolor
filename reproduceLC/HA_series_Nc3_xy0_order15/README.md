# Nc³ layer (C3F0Z0), HA1–6: expansion around (x, y) = (0, 0) to order 15

The Nc³ layer is the largest of the three leading-colour layers: up to 2227 G terms per HA,
against at most 736 for Nc² Nf. Each file `C3F0Z0_HA<i>_series15.m` defines
`HAser["C3F0Z0", i]`.

## What is expanded

The input is the GPL form `HA_GPL/C3F0Z0_HA<i>_G.m`, which was checked against the 2024
ancillary to 1e-36. Everything is expanded: the G functions **and** the rational prefactors.

The procedure is the one in `convertLC.nb`:

    ExpandPolyLogs[ExpandPolyLogs[expr, {x,0,15}], {y,0,15}]

That is, the expansion is iterated:
1. first in x to x^15 at fixed y (x << y);
2. then every coefficient in y to y^15.

All coefficients are exact rationals. The constants (Pi, I, Zeta[3], Zeta[5]) stay symbolic.

## Form of the result

    sum  c * const * x^m * y^n * Log[x]^kx * Log[y]^ky,   0 <= m <= 15,  n <= 15

The x-expansion of 1/(x+y)^k and of G(..,-y,..;x) produces powers of x/y, i.e. negative powers
of y, in the intermediate steps. They cancel in the sum: the lowest surviving power is y^-1
(y^0 for HA4). The 1/x and 1/x^2 prefactors also cancel, so no negative powers of x survive.

| HA | terms |
|----|------:|
| 1 | 5119 |
| 2 | 11424 |
| 3 | 11184 |
| 4 | 6702 |
| 5 | 11424 |
| 6 | 11382 |

## Check (`check.log`)

The series is compared with the exact G form evaluated with GiNaC at 40 digits. The table gives
the relative difference.

| (x, y) | 1/1000, 1/50 | 1/100, 1/20 | 1/50, 1/10 |
|---|---|---|---|
| HA1 | 1e-30 | 6e-24 | 8e-19 |
| HA2 | 1e-29 | 9e-23 | 2e-17 |
| HA3 | 2e-30 | 1e-23 | 2e-18 |
| HA4 | 7e-26 | 2e-19 | 1e-14 |
| HA5 | 1e-27 | 5e-21 | 5e-16 |
| HA6 | 8e-28 | 2e-21 | 2e-16 |

As expected, the truncation error is controlled by (x/y)^16 and y^16.

## Code

- `xyseries.py <in_G.m> <out.m> <name> [N]`: the expansion. It uses exact rationals, with
  python-flint for the y-polynomials, and takes about 20–70 s per HA.
- `check_series.py <in_G.m> <series.m> [points]`: the numerical check. It needs `ginsh`.
