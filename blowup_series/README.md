# blowup_series -- family and trident series in blow-up variables

Series of the master integrals of all 22 families (A B1 B2 C D E1 E2 F1 F2 G1 G2 G3 G4 H I1 I2 I3
L3 M3 N3 P3 R2) and of the 1318 trident functions, weights 0..6, at the three expansion points,
**left in the blow-up variables** (t, BIG). The helicity coefficients have spurious poles
(e.g. 1/(x+y)), so they are multiplied with the functions in these variables. Only the assembled
amplitude goes back to (x, y, z).

The recursion is the array/modular pipeline of `series.all/all_families_array` (python3 + numpy +
gmpy2, no Mathematica/FORM), with the cache-key fix. New here:

- `code/build_letters.py` builds the 20 letters at any order and in any chart. It replaces
  `00_build_letters.wl`, which limited the families to 15/45.
- Mirror charts (4-6) for the crossed families, see below.
- `code/select_trident_rows.py` + `code/trident_rotation.py`: families -> tridents.
- `run_all_families.py` reads this layout, writes no phys files and deletes the prime residues
  after each corner.

## Points and charts

There are three expansion points. The amplitude is expanded in the standard charts 1-3:

| point | chart | BIG | other variable | phys map of a term |
|---|---|---|---|---|
| x, y -> 0 | 1: x = t y^2 | y | z = 1 - y - t y^2 | t^e y^m -> x^e y^(m-2e) |
| y, z -> 0 | 2: z = t y^2 | y | x = 1 - y - t y^2 | t^e y^m -> z^e y^(m-2e) |
| x, z -> 0 | 3: x = t z^2 | z | y = 1 - z - t z^2 | t^e z^m -> x^e z^(m-2e) |

Log[t] -> Log[SMALL] - 2 Log[BIG]. Exponents can be half-integers (algebraic masters).
In the helicity-coefficient files, t = y and u = z (so x = 1 - t - u): their `t` is not the blow-up `t`.

**Mirror charts 4, 5, 6** are the same three points with SMALL and BIG exchanged:
4: y = t x^2, 5: y = t z^2, 6: z = t x^2. They are needed because the tridents are defined through
the crossed families. A crossing permutes x, y, z and can swap the two small variables: e.g. at
x, y -> 0, BASISLC[F, "x23", i](x, y) = F_i(y, x) is F at its own point but in the chart y = t x^2.
In (x, y) this is a relabelling, but the truncated region of a chart (x^a y^b with a <= ORDt,
b <= ORDy - 2a) is not symmetric, so the swapped function needs F's mirror-chart run. Without the
mirrors each point reaches rank 1310 of 1318, and HA3 uses the missing directions. A mirror run
uses the same atilde, sol0 and boundary constants (the log-free x^0 y^0 terms, the same in both charts).

## Order

A final (x, y) series to order N needs the functions up to t^(N+q), BIG^(3N+p). Here q and p are the
largest poles in t and BIG of the coefficients in the chart. For `C3F0Z0_HA3.wl`: q = 1 at every
point, p = 4 (x,y), 3 (y,z), 3 (x,z). For N = 20 that is (21, 64). **Production order: 22/66.**
`letters/` is built to 25/75 in all six charts, so any order up to that works.

## Run (cluster)

    python3 -m venv venv && venv/bin/pip install -r requirements.txt     # once
    venv/bin/python run_all_families.py --dry-run --ordt 22 --ordy 66     # show the order of families
    nohup venv/bin/python run_all_families.py -n 64 --ordt 22 --ordy 66 -w 6 > run_22_66.log 2>&1 &
    # then the tridents at the three points, weights 0..6
    venv/bin/python code/trident_rotation.py --dist dist_22_66 --out trident/dist_22_66 -n 64

    # single family / chart, by hand
    venv/bin/python code/run_array.py -w 6 -n 64 -c 1 --ordt 22 --ordy 66 \
        --workdir families/I3/dist_22_66 --inputs families/I3/inputs --letters letters --clean-residues

`run_all_families.py` runs all six charts by default. Everything can be resumed: rerun the same command,
and a chart whose weights are all reconstructed is skipped. Use `--keep-residues` only if you want to
extend to weight 7 later. Rough scale: H at 20/60 took about 35 min per chart on newton1. 22/66 has
about 1.3x the grid of 20/60 and needs somewhat more primes.

Output per family, in `families/<F>/dist_22_66/`:
- `exact_c<c>_w<w>.m`: one series per master, in t, BIG, Log[t], Log[BIG], Pi, Zeta, I (c = 1..6).
- `rec_c<c>_w<w>.pkl`: the same content as `{(tag, master, it, iy): mpq}`, faster to load.

Output of the rotation, in `trident/dist_22_66/`: `trident_c<c>_w<w>.pkl` and `.m` (c = 1..3, w = 0..6),
the 1318 trident series in the blow-up variables of point c.

## Tridents

`trident/rotation_rows.json` (12 MB) holds 1318 independent rows of the ancillary matrix m
(`BASISLC[F, crossing, i] = phase * m . Trident`), chosen cheapest first by `select_trident_rows.py`
from `data/trident_definition/trident_basis`. 931 tridents are single family elements, the rest
small combinations. `trident_rotation.py` reads, at each point, every selected crossed element from
the family chart it needs, solves m_sel . Trident = BASISLC_sel / phase modulo primes, reconstructs
and verifies. The 7 algebraic tridents (1, 2, 3, 15, 16, 17, 1283) come out as I * (real series).
They are divided by I to get the real, phase-normalised functions that the published series and the
amplitude coefficients use.

Crossings (x = s12, y = s13, z = s23): BASISLC[F, crossing, i](x, y, z) = F_i(X, Y, Z) with (X, Y, Z) =
"": (x,y,z), x12: (x,z,y), x13: (z,y,x), x23: (y,x,z), x123: (z,x,y), x132: (y,z,x).

## Amplitude series

    venv/bin/python code/assemble_amplitude.py ../trident_coefficients/*.wl \
        --tridents trident/dist_22_66 --out ../amplitude_series -n 16

For every coefficient file `<name>.wl` (= sum R(t,u) Trident[i, w] + I Pi sum R'(t,u) Trident[i, w]) and
every point, `assemble_amplitude.py`:
- expands the coefficients in the point's chart. It asserts that the chart resolves every
  denominator, and reports the poles q, p and the final order N = min(22 - q, (66 - p)/3), i.e. 20 for HA3;
- contracts with the trident series modulo primes (BLAS), then reconstructs and verifies;
- checks that no power below what the physical denominators allow survives (the spurious poles cancelled);
- writes two plain Mathematica expressions (`Get[...]`), each with a one-line comment header:
  - `<name>_c<c>_xyz_ord<N>.m`: physical variables (c1: x, y; c2: z, y; c3: x, z), rectangle of order N;
  - `<name>_c<c>_ber_ord<N>.m`: Bernoulli variables v -> 1 - Exp[-s_v] (sx, sy, sz), order N in each,
    with the LogX, LogY, LogZ placeholders of `05_bernoulli.wl`.

The trident residues are cached per prime in `trident/dist_22_66/residues/` and reused by every
coefficient file. At 22/66, one weight-6 residue set is a few GB, and each of the `-n` jobs holds one.

## Full cluster sequence

    venv/bin/python run_all_families.py -n 64 --ordt 22 --ordy 66 -w 6          # 22 families x 6 charts
    venv/bin/python code/trident_rotation.py --dist dist_22_66 --out trident/dist_22_66 -n 64
    venv/bin/python code/assemble_amplitude.py ../trident_coefficients/*.wl \
        --tridents trident/dist_22_66 --out ../amplitude_series -n 16

## Continued regions (A > 1)

The charts are in `code/charts.py`: R<A>_<b>v (s = 0) and R<A>_v<b> (s = 1) at the high-energy corner of each region.

- **Boundary constants.** `transport/` computes them for all families: the threshold corner (old CC1), the
  s = 0 charts (old CC2/CC3) and the s = 1 charts (old CC2e/CC3e). See `transport/README.md`.
- **Threshold corner series.** The threshold corner of region A has no chart of its own. Its series is the decay
  corner's series continued:

      venv/bin/python code/continue_threshold.py --workdir families/<F>/dist_22_66      # -> exact_cont_c<c>_w<w>.m

## Letters

    venv/bin/python code/build_letters.py 25 75                 # -> letters/, all six charts (about 15 s)
    venv/bin/python code/build_letters.py 15 45 --out /tmp/l --compare <old letters dir>

Each build checks itself numerically against the true letters (40 digits). It also asserts that the
chart resolves every letter (after pulling out t^a B^b the rest is nonzero at the point).

## Validated (exact, coefficient by coefficient)

- Letters at 15/45: all 20 letters, charts 1-3, identical to `series.all/all_families_array/letters`
  (00_build_letters.wl). With the principal root for l7, l8 in every chart.
- Letters at 25/75: the 14 H letters, charts 1-3, identical to `series.all/blow_up_array/pipeline/inputs`.
- Families with the new letters, weights 1-6, against the cluster phys results (`cluster_results.zip`)
  in the window the order determines: M3 at 6/18, I3 and H at 5/15, charts 1-3 **and mirrors 4-6**
  (a mirror run compared with the standard phys file of its point). All identical.
- Tridents: all 22 families at 5/15 in the six charts, rotated, against the published trident series
  (`data/series_trident_all_weights`, order 20) at all three points, weights 0-6, window 5 x 5
  (`code/validate_tridents.py`). All identical.
- Amplitude: `C3F0Z0_HA3.wl` from the 5/15 tridents. At all three points the spurious poles cancel.
  At x, y -> 0 all 780 terms of the window x, y <= 3 are identical to the independent expansion of
  the GPL form (`reproduceLC/HA_series_Nc3_xy0_order18`, HAser["C3F0Z0", 3]). The Bernoulli files are
  identical to Mathematica's `Series` of the physical ones, at all three points.
