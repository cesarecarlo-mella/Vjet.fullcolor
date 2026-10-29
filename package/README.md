# Zgqq_evaluator: full-colour three-loop finite remainders for V → q̄ q g

This package evaluates the **three-loop finite remainders of the helicity amplitudes for V → q̄ q g in full colour**,
at any point of the decay region and, analytically continued, of the three production regions. It gives all six
helicity coefficients in all 14 colour layers, with an error estimate for every value. A point takes about 0.06 s, or 0.12 s with the error estimate. Loading the module once
takes another 0.5 s.

## What is computed

- **Process:** V_L → q̄_R(p1) + q_L(p2) + g(p3), for both gluon helicities, as in `decay.m` of the planar ancillary
  files [1].
- **Kinematics:**
  - x = s12/q², y = s13/q², z = s23/q² = 1 − x − y, with q² = M_V²;
  - all results are at q² = 1, as in [1];
  - the decay region is x > 0, y > 0, x + y < 1;
  - the production regions are described under "Production (scattering) regions" below.
- **Coefficients:** α1, α2, α3, γ1, γ2, γ3, the helicity amplitude coefficients of [1] in the same order.
- **Colour:** every coefficient is the sum of 14 colour layers Nc^c Nf^f NfZ^z:

  | Nc power | layers |
  |---|---|
  | Nc³ | Nf⁰ NfZ⁰ |
  | Nc² | Nf¹, NfZ¹ |
  | Nc¹ | Nf⁰, Nf¹ NfZ¹, Nf² |
  | Nc⁰ | Nf¹, NfZ¹ |
  | Nc⁻¹ | Nf⁰, Nf¹ NfZ¹, Nf² |
  | Nc⁻² | Nf¹, NfZ¹ |
  | Nc⁻³ | Nf⁰ |

  - NfZ multiplies the contributions in which the vector boson couples to a closed quark loop; see [1] for the
    precise definition.
  - The leading-colour layer (Nc³) is the planar result of [1]: it reproduces α[i,3] and γ[i,3] of `decay.m` to 1e-13.

## Install

You need a C++17 compiler and Python 3 with numpy (and mpmath for the quad option). Boost.Odeint (header-only) is
included.

```sh
tar xzf Zgqq_evaluator-1.2.0.tar.gz && cd Zgqq_evaluator-1.2.0
make                       # builds core/cpp/tridents_bs and core/cpp/iter_quad (a few seconds)
make test                  # installation test, about 30 s; PYTHON=/path/to/python3 if needed
make test-quad             # test of the quad-precision option, about 2 min (needs mpmath)
```

The test checks five things:
1. decay: the leading-colour layer at 45 points against the exact published values;
2. decay: all 84 values at 12 points against the reference values of this release;
3. production regions: all 84 values at 3 points of each region against the reference values of this release;
4. production region y: the planar layers (Nc³, Nc² Nf, Nc Nf²) against the published production results [2];
5. the command line.

Tested on macOS (arm64, Apple clang). On Linux, use `make CXX=g++`.

## Use

**Command line:**
```sh
python3 zgqq_eval.py 0.2 0.3                          # all 84 values with error estimates
python3 zgqq_eval.py --points pts.txt --out res.csv   # one point "x y" per line; .csv or .json output
python3 zgqq_eval.py 0.2 0.3 --Nc 3 --Nf 5 --NfZ 1    # also the six coefficients summed over the colour layers
python3 zgqq_eval.py 0.2 0.3 --no-error               # twice as fast, without the error estimate
python3 zgqq_eval.py 0.2 0.3 --quad                   # quad precision: all values printed with 32 digits
python3 zgqq_eval.py 0.2 0.3 --region z               # production region z > 1, (X, Y) of its triangle (see below)
python3 zgqq_eval.py --physical -- -0.5 -0.25         # physical (x, y) anywhere: decay or production, chosen per point
```

**Python:**
```python
import zgqq_eval as Z
res = Z.evaluate([(0.2, 0.3), (0.1, 0.5)])
res.values[p, i, l]       # complex: point p, coefficient Z.COEFFICIENTS[i], colour layer Z.LAYERS[l]
res.errors[p, i, l]       # absolute error estimate, same shape
v, e = Z.colour_sum(res, Nc=3, Nf=5, NfZ=1.0)          # shape (npts, 6): the summed coefficients and their errors

q = Z.evaluate([(0.2, 0.3)], precision='quad')         # quad precision (double-double, about 32 digits)
q.values + q.values_lo    # the full value: q.values is it rounded to double, q.values_lo the rest

r = Z.evaluate([(0.2, 0.3)], region='z')               # production region z > 1, (X, Y) of its triangle
p = Z.evaluate_physical([(-0.5, -0.25), (0.2, 0.3)])   # physical (x, y): the region is chosen per point (p.region)
Z.region_point(-0.5, -0.25)                            # -> ('z', X, Y);  Z.physical_point('z', X, Y) -> (x, y, z)
```

## Production (scattering) regions

The production regions have one invariant above q² and the other two negative. Each is parametrised by its own
triangle (X, Y), X > 0, Y > 0, X + Y < 1, with Z = 1 − X − Y the inverse of the large invariant (Z → 0 is the
high-energy limit):

| region | (x, y, z) | large invariant | incoming partons (decay labels p1 = q̄, p2 = q, p3 = g) |
|---|---|---|---|
| `'x'` | (1/Z, −X/Z, −Y/Z) | s12 > q² | p1, p2: the quark-antiquark channel |
| `'y'` | (−X/Z, 1/Z, −Y/Z) | s13 > q² | p1, p3: the channel of `lmL.m`/`lpL.m` in [2] (q̄ g → V q̄ there) |
| `'z'` | (−X/Z, −Y/Z, 1/Z) | s23 > q² | p2, p3: the other quark-gluon channel |

- **What is returned:** the decay coefficients α_i, γ_i analytically continued into the region, with the Feynman
  prescription (every invariant + i0). The root r = √(xyz) continues to −√|xyz|.
- **Not done here:** the crossing itself, i.e. relabelling momenta and helicities and the spinor phases of the
  incoming partons, is left to the user.
- **Cross-check:** in region y these are directly the coefficients of the published production results [2] (`lmL.m`:
  α[i,3], `lpL.m`: γ[i,3]): the planar layers agree to 1.5e-11 (test 4).
- **Precision:** both double and quad work in all regions. The start values at each region's hub (0.3, 0.3) come from
  a quad integration along the i0 path from the decay start values; two paths agree to 1e-27. They agree with an
  independent continuation through the corner to 1e-11.
- **Accuracy:** in double, the same as in the decay triangle: about 1e-12 in the bulk, degrading towards the edges of
  the region's triangle. Note that X, Y → 0 is a collinear/soft corner of the production kinematics and Z → 0 is high
  energy.

`Z.LAYERS` are strings like `'Nc^-1*Nf^1*NfZ^1'`, and `Z.POWERS` gives their exponents (c, f, z).

## Accuracy

**Error estimate.** Every point is evaluated twice:
- the two evaluations start from independent boundary values and integrate along different complexified paths;
- the reported error is 10 × their difference.

Against exact values, the actual error was below the estimate for 100% of the values at 406 points across the
triangle, and for 99% at 30 points near the corner x, y → 0. The estimate is typically 3–10 times larger than the
actual error.

**Where the evaluator is accurate.**
- In the bulk, the relative accuracy is about 1e-12, and at worst 1e-9.
- Close to the edges of the triangle, the rational functions (up to 1/x³, 1/y³, 1/z³, and (x+y)⁻⁵ etc. at the
  corners) cancel against large logarithms, and double precision loses digits.

Typical worst relative error over the 84 values, by distance d from the edge:

| d | x → 0 | y → 0 or z → 0 | corners |
|---|---|---|---|
| 0.1 | 1e-10 | 1e-9 | 1e-9 |
| 0.03 | 1e-9 | 1e-7 | 1e-5 |
| 0.01 | 1e-9 | 1e-5 | 1e-3 |
| 0.003 | 1e-8 | 1e-4 | O(1) |

Check the error estimate when you go close to the edges. The planned corner expansions will cover those regions.

**Quad precision** (`precision='quad'`, `--quad`). Everything is evaluated in double-double arithmetic (about 32
digits): the integration, the 579 functions, the rational functions and the sums.
- The start values at (0.3, 0.3) hold 27–28 digits. They come from an exact double series at the corner, transported
  in quad from two start points that agree to 3e-28.
- The error estimate comes from two paths. It is typically 1e-26 to 1e-23 relative in the bulk. Near the edges quad
  loses digits to the same cancellation as double, e.g. 1e-18 at d = 0.004.
- Cost: 1–10 s per point and path, against 0.07 s in double.

Mixed precision (double functions, quad rational functions) is not offered: it gains at most a factor of a few near
the edges, because the loss is in the functions themselves.

## How it works

**Representation.** Every colour layer is Σ q · R_r(y, z) · T_g · c_m:
- R_r: 129 rational functions;
- T_g: 579 pure functions, combinations of iterated integrals over a 20-letter alphabet (base point x = y = 0,
  shuffle regularised);
- c_m: 15 constants built from π, ζ3, ζ5 and i. One of them is i/r, with r = +√(xyz), which multiplies the functions
  that are odd under r → −r.

The analytic data are in `amplitudes/`, in Mathematica format. The row format is explained in the header of
`amplitudes/FiniteRemainders/*.wl`.

**Numerical evaluation.**
- **The system.** The 579 functions satisfy a closed system of differential equations: 503 independent functions,
  canonical form dF = Σ_k dlog W_k M_k F.
- **The integration.** It is solved numerically from boundary values at (0.3, 0.3), with a Bulirsch–Stoer method
  along a complex-deformed straight path, for every point (`core/cpp/tridents_bs.cpp`).
- **Spurious letters.** Six letters of the alphabet vanish on curves inside the triangle. The amplitudes are smooth
  there, but the differential equations (and some rational functions) have poles.
  - The path deformation is chosen so that the path keeps away from these zeros.
  - A point lying on such a curve is replaced by a symmetric average across it, with O(h⁶) accuracy.
  - Both happen automatically; `res.info` records them.
- **Quad precision.** The 579 functions are integer combinations of 586 functions with an integer system of
  differential equations. That system is integrated in double-double arithmetic (`core/cpp/iter_quad.cpp`). A point
  on a spurious curve is averaged with a 5-level stencil, O(h¹⁰).

## Validation

| check | result |
|---|---|
| Nc³ layer vs the exact published values [1], 406 points across the triangle | median relative difference 1e-13, worst 7e-10 |
| Nc³ layer vs its exact expansion around x, y → 0, 30 points | consistent with the error estimate |
| points exactly on the spurious curves: step h vs h/2, and vs direct evaluation 1e-4 away | 1e-10 to 3e-9 |
| the representation vs an independent form of the full-colour remainders | identical, coefficient by coefficient |

| production regions: start values at the hub, two quad i0 paths | 1e-27 |
| production regions: quad i0 continuation vs the continuation through the corner (independent) | 1e-11 |
| region y, planar layers vs the published production results [2], 3 points x 6 coefficients x 3 layers | worst 1.5e-11 |
| production regions: double vs quad, 6 points per region | median 2e-12; within the error estimate 504/504 per region |
| production regions, all 84 values at one point each: evaluator vs direct quad integration from the decay along an i0 path (no region system) | real and imaginary parts, every colour layer: 1e-28 to 5e-23 |
| the same vs the independent continuation through the corner (double) | 1e-13 to 2e-8 (the accuracy of that double continuation) |

In the decay, the three subleading layers that contain the i/r constant (Nc¹ Nf⁰, Nc⁻¹ Nf⁰, Nc⁻³ Nf⁰) have been
checked only against the quad evaluation of this package (the same analytic input), not against an independent
computation. In the production regions their analytic continuation is checked as above
(script `iter int /generators_deq/region_continuation_check.py` of the development repository).

## Limitations

- The crossing of the production regions (momentum and helicity relabelling, phases) is left to the user.
- Points must lie strictly inside the triangle (decay, or a region's (X, Y)), and the accuracy degrades near its edges
  (see above).
- The error estimate in double is conservative: it includes a term for the systematic error of the double system,
  amplified by the cancellations in the sum over terms, which overestimates the actual error ~1000 times in the bulk.

## Files

| path | content |
|---|---|
| `zgqq_eval.py` | the interface (Python and command line) |
| `core/` | the evaluator: `gen_eval.py` (the 579 functions), `remainder_eval.py` (the remainders), `cpp/` (the C++ integrators), `out/` (the differential equations and boundary values); quad: `quad_eval.py`, `remainder_quad.py`, `quad_data.py`, `ddnp.py`, `out/iter_quad_P0.bin`, `out/P0_values_quad.m` |
| `amplitudes/` | the analytic representation: `FiniteRemainders/{alpha,gamma}{1,2,3}.wl`, `Generators.wl`, `RationalFunctionBasis.wl`, `TranscendentalMultipliers.wl`, `Alphabet20.wl` |
| `tests/` | the installation test and its reference data |
| `MANIFEST.sha256` | checksums of all files |

## References

[1] T. Gehrmann, P. Jakubčík, C. C. Mella, N. Syrrakos, L. Tancredi, *Planar three-loop QCD helicity amplitudes for
V+jet production at hadron colliders*, and its ancillary files (`decay.m`).

[2] The same ancillary files, `FiniteRemainders/lmL.m` and `lpL.m`: the planar production results (q̄ g → V q̄).

Reference for the full-colour results: to be added.
