# ITER -- the iterated-integral system as a 23rd family

The closed system of the 579 generators of `Zgqq_finite_remainders/Generators.wl`, in the basis of
`iter int /generators_deq/out/basis.wl`: integral, and redundant at weights 6 and 5. It has 586 functions; the generators in
terms of them are in `generators.wl` there. Every entry of the connection is an integer, at most 3.5e5; the 6 -> 5 block
holds only scale factors <= 6. The weight plays the role of the eps order.

Inputs (`inputs/`), from `iter int /generators_deq/family_inputs.py` and `transport/iter_corner_transport.wl`:

| file | content |
|---|---|
| `atilde.m` | sum_k M[k] l[k] (586 x 586) |
| `labels.m` | the functions in order |
| `sol0.m` | weight 0: B[0,1] = 1/96 (One = 96 B[0,1]) |
| `bdy_c1_w*.m` | 0: the shuffle-regularised base point, chart 1 (x = t y^2) |
| `bdy_c3_w*.m` | exact transport along x = 0 (GPLs, PolyLogTools), chart 3 (x = t z^2) |
| `bdy_c2_w*.m` | exact transport along y = 0, in the regularisation of the **mirror chart 5** (y = t z^2) |
| `edge_x0.m`, `edge_y0.m` | the connection on the two edges, {A0, A1} with dF/ds = (A0/s + A1/(s-1)) F |

**Charts: use 1, 3 and 5 only.**
- The constants at c2 depend on the chart. The letters (y+z)^2 - z and (y+z)^2 - y vanish non-transversally there, so the
  log-free terms differ between chart 5 (y << z) and chart 2 (z << y).
- The pipeline reads `bdy_c2` for both charts. Chart 2 with these constants is wrong: O(1) differences from the
  numerical evaluator. Chart 5 resolves all the rational functions at c2, so chart 2 is not needed.
- The mirrors 4 and 6 have not been checked either.

    python3 code/run_array.py -w 6 -n 2 -c 1,3,5 --ordt 5 --ordy 15 --workdir families/ITER/dist_5_15 \
        --inputs families/ITER/inputs --letters letters --clean-residues
    python3 code/validate_iter.py families/ITER/dist_5_15 families/ITER/inputs "../iter int /generators_deq/out/generators.wl" 1,3,5

**Validation at 5/15** (`code/validate_iter.py`): the series of the 579 generators against the numerical evaluator, at
t = 0.01, BIG = 0.04 and at t = BIG = 0.1:

| chart | median relative difference | worst |
|---|---|---|
| 1 | 3e-13 | 2e-8 |
| 3 | 2e-13 | 3e-10 |
| 5 | 3e-13 | 3e-9 |

The worst values are at the point closest to the edge, where the numerical side is the less precise one.
