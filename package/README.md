# Zgqq_evaluator: full-colour three-loop finite remainders for V → q̄ q g

This package evaluates the **three-loop finite remainders of the helicity amplitudes for V → q̄ q g in full colour**,
at any point of the decay region. It gives all six helicity coefficients in all 14 colour layers, with an error
estimate for every value. A point takes about 0.06 s, or 0.12 s with the error estimate. Loading the module once
takes another 0.5 s.

## What is computed

- **Process:** V_L → q̄_R(p1) + q_L(p2) + g(p3), for both gluon helicities, as in `decay.m` of the planar ancillary
  files [1].
- **Kinematics:**
  - x = s12/q², y = s13/q², z = s23/q² = 1 − x − y, with q² = M_V²;
  - all results are at q² = 1, as in [1];
  - the decay region is x > 0, y > 0, x + y < 1.
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

You need a C++17 compiler and Python 3 with numpy. Boost.Odeint (header-only) is included.

```sh
tar xzf Zgqq_evaluator-1.0.0.tar.gz && cd Zgqq_evaluator-1.0.0
make                       # builds core/cpp/tridents_bs (a few seconds)
make test                  # installation test, about 20 s; PYTHON=/path/to/python3 if needed
```

The test checks three things:
1. the leading-colour layer at 45 points against the exact published values;
2. all 84 values at 12 points against the reference values of this release;
3. the command line.

Tested on macOS (arm64, Apple clang). On Linux, use `make CXX=g++`.

## Use

**Command line:**
```sh
python3 zgqq_eval.py 0.2 0.3                          # all 84 values with error estimates
python3 zgqq_eval.py --points pts.txt --out res.csv   # one point "x y" per line; .csv or .json output
python3 zgqq_eval.py 0.2 0.3 --Nc 3 --Nf 5 --NfZ 1    # also the six coefficients summed over the colour layers
python3 zgqq_eval.py 0.2 0.3 --no-error               # twice as fast, without the error estimate
```

**Python:**
```python
import zgqq_eval as Z
res = Z.evaluate([(0.2, 0.3), (0.1, 0.5)])
res.values[p, i, l]       # complex: point p, coefficient Z.COEFFICIENTS[i], colour layer Z.LAYERS[l]
res.errors[p, i, l]       # absolute error estimate, same shape
v, e = Z.colour_sum(res, Nc=3, Nf=5, NfZ=1.0)          # shape (npts, 6): the summed coefficients and their errors
```

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

## Validation

| check | result |
|---|---|
| Nc³ layer vs the exact published values [1], 406 points across the triangle | median relative difference 1e-13, worst 7e-10 |
| Nc³ layer vs its exact expansion around x, y → 0, 30 points | consistent with the error estimate |
| points exactly on the spurious curves: step h vs h/2, and vs direct evaluation 1e-4 away | 1e-10 to 3e-9 |
| the representation vs an independent form of the full-colour remainders | identical, coefficient by coefficient |

The three subleading layers that contain the i/r constant (Nc¹ Nf⁰, Nc⁻¹ Nf⁰, Nc⁻³ Nf⁰) have no independent
numerical check yet.

## Limitations

- Decay kinematics only. The production channels are not included.
- Points must lie strictly inside the triangle, and the accuracy degrades near its edges (see above).

## Files

| path | content |
|---|---|
| `zgqq_eval.py` | the interface (Python and command line) |
| `core/` | the evaluator: `gen_eval.py` (the 579 functions), `remainder_eval.py` (the remainders), `cpp/` (the C++ integrator), `out/` (the differential equations and boundary values) |
| `amplitudes/` | the analytic representation: `FiniteRemainders/{alpha,gamma}{1,2,3}.wl`, `Generators.wl`, `RationalFunctionBasis.wl`, `TranscendentalMultipliers.wl`, `Alphabet20.wl` |
| `tests/` | the installation test and its reference data |
| `MANIFEST.sha256` | checksums of all files |

## References

[1] T. Gehrmann, P. Jakubčík, C. C. Mella, N. Syrrakos, L. Tancredi, *Planar three-loop QCD helicity amplitudes for
V+jet production at hadron colliders*, and its ancillary files.

Reference for the full-colour results: to be added.
