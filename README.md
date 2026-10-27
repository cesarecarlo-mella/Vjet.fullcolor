# Vjet.fullcolor

Full-colour V+jet helicity amplitudes, and their series expansions at the three corners
(c1: x, y -> 0; c2: y, z -> 0; c3: x, z -> 0; x + y + z = 1).

## Layout

| folder | content |
|---|---|
| `trident_coefficients/` | Helicity coefficients in the trident basis: `C<c>F<f>Z<z>_HA<i>.wl` = sum of rational(t, u) * `Trident[i, w]`, with t = y, u = z. |
| `denominator_bracket/` | Helicity coefficients, denominator-bracket form (+ `letter_stats/`). |
| `remainder_denominator_bracket/` | Finite remainders, denominator-bracket form. `reproduceLC/HA_GPL` reads it as `../../remainder_denominator_bracket`, so keep it at the top level. |
| `reproduceLC/` | Leading-colour reproduction and HA3 series checks. |
| `blowup_series/` | Master series of all 22 families in blow-up variables (array/modular pipeline, cluster-ready). See its README. |
| `tools/` | Physical-variable series converters: JSON -> Mathematica lists (`convert_*`), extraction of a `Trident[i, w]` selection (`extract_trident_selection.py`, `trident_selection.m`). |
| `poles/` | Pole predictions and checks (the large `.m` files are not in git). |
| `inversion_reducible/` | Kira reduction targets. |
| `data/` | Large local data, **not in git**: `download/` (raw zgqq amplitudes), `series_trident_all_weights/` (trident series, JSON, physical variables, order 20), `mathematica-results/` (the same as `.m`, per family and `trident/`), `_xfer/` (transfer copies). |

## Series of the amplitudes (plan)

The rational coefficients have spurious poles at the corners (x + y at c1, y + z at c2, x + z at c3), so
the series cannot be multiplied in (x, y) with a rectangular truncation. Coefficients and trident
functions are both expanded in the blow-up variables of `blowup_series` (c1: x = t y^2, c2: z = t y^2,
c3: x = t z^2) and multiplied there. Only the assembled amplitude goes back to (x, y).

1. `blowup_series/run_all_families.py`: family masters in blow-up variables at 22/66 (enough for
   N = 20), in the three charts and their mirrors (needed for the crossed families).
2. `blowup_series/code/trident_rotation.py`: the family solutions -> `Trident[i, w]` at the three
   points, with the trident definition in `data/trident_definition/` (only 1318 rows of it, in
   `blowup_series/trident/rotation_rows.json`).
3. `blowup_series/code/assemble_amplitude.py`: expansion of the coefficients in blow-up variables,
   contraction M . Trident[i, w], cancellation check, output in physical and Bernoulli variables.

All three steps are validated at low order (5/15) against the published trident series and the GPL
expansion of HA3. See `blowup_series/README.md`. Production runs on the cluster.

## Tools usage

    python3 tools/convert_trident_to_mathematica.py data/series_trident_all_weights data/mathematica-results/trident
    python3 tools/extract_trident_selection.py tools/trident_selection.m data/mathematica-results/trident data/mathematica-results/trident

Python environment for `blowup_series`: `python3 -m venv .venv && .venv/bin/pip install -r blowup_series/requirements.txt`.
