# generators_deq -- closed DE and numerical evaluation of the 60 generators

The generators H[1..60] of `../Generators.wl` are pure combinations of iterated integrals II[a1..an] over the
20-letter alphabet (`data/trident_definition/alphabet.wl`), weights 2-6. Conventions:

- d II[a1..an] = dlog W_an II[a1..a(n-1)].
- Base point c1, shuffle regularised: II[] = 1 and every non-empty II -> 0.

## System

Close the generators under d: at each weight, the generators plus the strips (last letter removed) of the weight
above. This gives

| weight | 6 | 5 | 4 | 3 | 2 | 1 | 0 |
|---|---|---|---|---|---|---|---|
| functions (generators) | 16 (16) | 88 (15) | 97 (16) | 72 (9) | 23 (4) | 5 | 1 |

302 functions in total.

- **Exact DE:** dF = Sum_k dlog W_k M[k] . F, with the generators first at each weight. It is in `out/deq.wl`
  (labels and M[1..20] as SparseArrays), and the functions F in II[...] are in `out/basis.wl` (`exact_deq.py`).
  - Closed and checked word by word.
  - The rationals are large (up to about 54 digits) because of the choice of extra functions.
- **Numerical gauge:** an orthonormal basis of each weight in word space. The letter blocks are then the strip maps,
  with every entry |A| <= 1.

## Numerical evaluation (method of `INTEGRALS/4P1M/3L1M/trident_deq`, Czakon-Tancredi)

`cpp/tridents_bs.cpp` is the trident_deq evaluator, unchanged:
- Bulirsch-Stoer on the complex-deformed straight path from P0 = (0.3, 0.3);
- r = Sqrt[xyz] is integrated along;
- the generators are read off at the end point.

The data are in `out/generators_bs.bin` (TRIDBS01 format, NT = 60, generator i in the slot of its weight).

```sh
make                                              # builds cpp/tridents_bs
python3 gen_eval.py 0.25 0.35                     # the 60 values at one point
python3 -c "from gen_eval import evaluate; H, info = evaluate([(0.3, 0.2), (0.1, 0.7)])"   # H.shape == (2, 60)
```

Accuracy and speed over the triangle (`gen_eval_check.py`), at the default tolerance 1e-10:
- relative error median 5e-13, worst 6e-9 (weight 6);
- 0.03 s per point.

Start values (`iisys.start_values`):
- log-power series from c1, along x = 0 and then in Sqrt[x] at y = 1/2, giving the values at (0.02, 0.5);
- carried once to P0 with a tight tolerance (`gen_run.py`).

## Checks

| check | result |
|---|---|
| 24 tridents through the same pipeline (`trident_check.py`) vs the exact corner series, 21 points near c1, c2, c3 | worst 7e-13 (eps^1), 5e-11 (eps^5), 3e-8 (eps^6), as good as trident_deq |
| x -> 0 first vs y -> 0 first regularisation at c1, and the two sides of the path deformation, all 302 functions at P0 | agree to 5e-12 (single-valued in the triangle) |
| direct c1 series near the edges vs integration from P0 | 4e-12 to 3e-11 |

## Files

| file | content |
|---|---|
| `iisys.py` | closure, orthonormal system, 1-d log-power series start values, TRIDBS01 export (works for any set of pure II combinations) |
| `gen_run.py` | build for Generators.wl: system, both regularisations, carry to P0, write `out/generators_bs.bin` |
| `gen_eval.py` | python interface |
| `gen_eval_check.py`, `trident_check.py` | the checks above (`trident_check.py` reads `2610.anc/trident_basis`) |
| `exact_deq.py` | the exact DE (`out/basis.wl`, `out/deq.wl`) |
| `out/system.pkl` | the numerical system (needed by `gen_eval.py` for the weights) |

To rebuild everything:

```sh
python3 gen_run.py ../Generators.wl out
python3 exact_deq.py out/system.pkl out
```
