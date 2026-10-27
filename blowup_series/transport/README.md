# transport -- boundary constants of the continued regions

Boundary constants for the series in the three continued regions A > 1 (A = x, y, z; the other two variables
negative), for all 22 families, weights 0-6. They are computed exactly with PolyLogTools. The charts and their
names are in `../code/charts.py`:

- Coordinates: A = 1/v and b = -bh/v for the other two variables, with v + bh + ch = 1.
- Each region has three corners:

| corner | where | series | constants | old name (A = z) |
|---|---|---|---|---|
| threshold | bh, ch -> 0 (v -> 1) | the decay corner continued (`../code/continue_threshold.py`) | `bdy_R<A>_th` (`threshold.py`) | CC1 |
| high energy, s = 0 | bh, v -> 0 with bh << v | chart `R<A>_<b>v` (bh = t v^2) | `bdy_R<A>_<b>v` (`edge_transport.wl`) | CC2 (b = y), CC3 (b = x) |
| high energy, s = 1 | bh, v -> 0 with v << bh | chart `R<A>_v<b>` (v = t bh^2) | `bdy_R<A>_v<b>` (`divisor_transport.wl`) | CC2e, CC3e |

Region A touches the decay triangle at the decay corner where the other two variables vanish:

- z > 1 at c1 (x, y -> 0)
- x > 1 at c2 (y, z -> 0)
- y > 1 at c3 (x, z -> 0)

Everything is transported out of that corner.

**Continuation convention (as for CC).** Log b -> log|b| + I Pi for both small variables, and Sqrt[x y z] -> -Sqrt|x y z|.

## The three steps

1. **`threshold.py`: threshold constants.** No DE is solved.
   - `low`: the x^0 y^0 coefficient of the decay series, a polynomial in Log[t] and Log[BIG], for every master.
     It reads a decay run at any order (`rec_c<c>_w<w>.pkl`) and writes `families/<F>/inputs/low_c<c>_w<w>.m`.
     These files are kept in the repo, so the old runs are no longer needed.
   - `consts`: the threshold constant is that polynomial at Log[t] -> -I Pi, Log[BIG] -> I Pi. Writes
     `bdy_R<A>_th_w<w>.m`; the w = 0 file is `sol0`.
2. **`edge_transport.wl`: the s = 0 charts.** Transport along the edge bh = 0 from the threshold corner to the
   high-energy corner.
   - On the edge the even letters are monomials in s and 1 - s, and the odd letters are constant.
   - GIntegrate, weight by weight, seeded by the threshold constants.
   - At s -> 1: s -> 1 - t, change of fibration (`fibrtto1.m`), drop G[0..0, t], then t -> 0.
3. **`divisor_transport.wl`: the s = 1 charts.** Transport along the exceptional divisor of the high-energy corner,
   with bh = s p, v = (1 - s) p, ch = 1 - p.
   - On the divisor the letters are rational in q, with s = q^2/(1 + q^2). The odd letter becomes
     (1 - I rho)/(1 + I rho) with rho rational in q. Its sign comes from r = -Sqrt|xyz| and the sign of a in the region.
   - G functions in q, from q = 0 (the s = 0 constants) to the regularised limit q -> oo, by inversion.
   - Divergent G(1, ...; 1) are shuffle-regularised; their coefficients must cancel, and the script checks this.
   - Whatever is left as G(...; 1), Catalan or Log[2] is fitted on Q[i Pi, MZV]. The values come from GiNaC at
     75 digits, and a fit is accepted only below 1e-60.

The scripts stop with `!!` on any failed check:
- a letter with a pole off {0, 1} on an edge;
- an odd letter that is not rational in q;
- a divergence that does not cancel;
- a constant that is not in Q[i Pi, MZV].

## Run

```sh
python3 threshold.py low --runs <dir with <F>/rec_c*_w*.pkl>      # once (done: low_c*_w*.m are in families/)
python3 threshold.py consts                                       # -> bdy_R{Z,X,Y}_th_w{0..6}.m (done)
python3 run_transport.py --dry-run                                # what is left
nohup python3 run_transport.py -n 16 > logs/run_transport.log 2>&1 &
```

- `run_transport.py` runs 66 (family, region) chains. Each chain is one edge job followed by two divisor jobs.
- Each job uses one Mathematica kernel, and `-n` jobs run at a time.
- It is resumable: a job whose last-weight file exists is skipped. Logs go to `logs/<F>_<A>_<step>.log`.
- PolyLogTools must be installed, and `ginsh` (GiNaC) must be on the PATH for the fits.
- Timing for H (371 masters): edge (both s = 0 charts) 1 min; one divisor job about 17 min.

## Validation (exact, all components)

| test | result |
|---|---|
| `bdy_RZ_th` vs `transport.boundaries/H/dist/bdy_CC1`, w 0-6 | identical |
| `bdy_RZ_yv`, `bdy_RZ_xv` vs `bdy_CC2`, `bdy_CC3`, w 0-6 | identical |
| `bdy_RZ_vy` vs `analytic_continuation/results/bdy_CC2e_H` (divisor transport + MZV fit), w 1-6 | identical (after PolyGamma[3,1/4] + PolyGamma[3,3/4] = 16 Pi^4) |
| edge jobs, all 22 families x 3 regions | 66/66 ok, every constant in Q[i Pi, Zeta[3], Zeta[5]] |
| `continue_threshold.py` (chart variables) vs direct continuation of the physical series, H, 6 charts, w 1-3 | identical |

`compare.py` (python) and `compare_exact.wls` (Mathematica, reduces the PolyGamma/Hurwitz identities that old files
leave unreduced) do the exact comparisons.

Note on the fit: the odd masters carry an extra I, so real and imaginary parts are each fitted on the full
weight-w basis of Pi^k x MZV, at 100 digits, accepting only relations with integers below 10^13.
