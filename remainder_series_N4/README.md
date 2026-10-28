# Series of the 84 finite remainders at the three decay corners (order N = 4, local test run)

These are the three-loop finite remainders of `Zgqq_finite_remainders`: coefficients alpha1-3 and gamma1-3, each with
14 colour layers. Each is expanded around the three corners of the decay triangle (x, y, z = 1 - x - y).

    <coefficient>/<colour layer>/seriesc1.m      corner c1: x, y -> 0, in x and y
                                 seriesc2.m      corner c2: y, z -> 0, in y and z   (x = 1 - y - z)
                                 seriesc3.m      corner c3: x, z -> 0, in x and z   (y = 1 - x - z)
                                 seriesBerc1.m ... the same in Bernoulli variables

- **Colour layers.** `Nc3_Nf0_NfZ0` is Nc^3 Nf^0 NfZ^0, and so on; the remainder is the sum over the layers of
  Nc^c Nf^f NfZ^z times the series.
- **Files.** Every file is a plain Mathematica expression (`Get[...]`), with a one-line comment header.
- **Native variables.**
  - A rectangle of order N = 4 in the two small variables, with Log[...] of the small variables and the constants Pi,
    Zeta[3], Zeta[5], I.
  - Powers can be negative (physical poles, e.g. y^-1).
- **Bernoulli variables.**
  - v -> 1 - Exp[-s_v] for the two small variables, to order N in each.
  - Log[x], Log[y], Log[z] are left as the placeholders LogX, LogY, LogZ (substitute at the point).

## How they were made

1. **The functions.** The closed system of the 579 iterated-integral generators (586 functions, `iter int
   /generators_deq/out/deq.wl`) was expanded as a 23rd family, ITER, in the blow-up charts (blowup_series, run_array.py,
   order 8/24):
   - c1: chart 1, x = t y^2;
   - c2: chart 5, y = t z^2;
   - c3: chart 3, x = t z^2.
2. **The boundary constants.**
   - c1: zero (the regularised base point).
   - c2, c3: exact transport along the edges y = 0 and x = 0 in GPLs (PolyLogTools); they are MZVs.
3. **Assembly in the blow-up variables** (`blowup_series/code/assemble_remainders.py`): rational functions times
   generators times constants, with I/r and r = Sqrt[x y z].
4. **Back to the variables of the corner**, only at the end. The spurious poles of the rational functions cancelled
   in all 252 series: no term below the physical poles survives.

## Checks

- **Exact, leading colour.** The Nc^3 layer of all six coefficients at c1 is identical, term by term, to the
  independent exact expansion of the GPL form (`reproduceLC/HA_series_Nc3_xy0_order18`): 6195 terms in the window
  x, y <= 4.
- **Numerical, all 84 layers at all three corners.** At distance 0.05 from the corner the series agree with the
  numerical evaluator (`iter int /generators_deq/remainder_eval.py`) within the size of the order-4 terms; the worst
  difference is 0.14 of that size (`blowup_series/code/validate_remainder_series.py`).

**Order.** N = 4 comes from the ITER series at 8/24; this is a local test. The production order N = 20 needs the ITER
series at 24/72 (cluster).
